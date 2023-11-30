#include <iostream>

#include "cgleaf.H"
#include "cgpairtree.H"
#include "cgslist.H"

REAL CLeaf::m_distances[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE];

// Switch the rotamer coordinates held by this leaf to different 
// precomputed values.
void CLeaf::changeRotamer(int rotIndex)
{
  assert(rotIndex < SIDECHAIN::m_aalist[getType()]->m_nRotamers);
  
  m_undoRotIndex = m_rotIndex;
  m_rotIndex = rotIndex;

  const ROTAMER & rot = 
    SIDECHAIN::m_aalist[getType()]->getRotamer(rotIndex, m_rotType);

  // Set the new coordinates as well as their BV and energy, which were
  // precomputed.
  m_positions = rot.m_positions;
  m_bv = rot.m_bv;
  m_energy = rot.m_energy;
  m_prob = rot.m_prob;
  m_chi1mean = rot.m_chi1mean;
  m_chi2mean = rot.m_chi2mean;
  m_chi3mean = rot.m_chi3mean;
  m_chi4mean = rot.m_chi4mean;
  m_chi1std = rot.m_chi1std;
  m_chi2std = rot.m_chi2std;
  m_chi3std = rot.m_chi3std;
  m_chi4std = rot.m_chi4std;
}

// Compute all pairs of distances between the atoms of this leaf and
// the given leaf.
void CLeaf::computeDistances(CLeaf * pLeaf, const REAL rot[3][3],
			     const REAL trans[3])
{
  int size1 = SIDECHAIN::m_aalist[getType()]->m_size;
  int size2 = SIDECHAIN::m_aalist[pLeaf->getType()]->m_size;
  
  REAL cen[3], dist[3], vec1[3], vec2[3];

  // If this is a PRO backbone, we need to add the Cd atom because 
  // it is part of a group with the Ca atom for Electrostatic purposes
  if (getType() == BBP)
    {
      assert(getNext());
      assert(getNext()->getType() == PRO);
      MxVpV(vec1, m_rotate, getNext()->getPositions()[2], m_translate);
    }

  // If this is a PRO backbone, we need to add the Cd atom because 
  // it is part of a group with the Ca atom for Electrostatic purposes
  if (pLeaf->getType() == BBP)
    {
      assert(pLeaf->getNext());
      assert(pLeaf->getNext()->getType() == PRO);

      REAL temp[3];
      MxVpV(temp, pLeaf->m_rotate, pLeaf->getNext()->getPositions()[2], 
	    pLeaf->m_translate);
      MxVpV(vec2, rot, temp, trans); 
    }
  
  // Compute te distances between all pairs of atoms.
  for(int j = 0; j < size2; j++)
    {
      MxVpV(cen, rot, pLeaf->getPositions()[j], trans); 

      for (int i = 0; i < size1; i++)
	{	
	  VmV(dist, cen, getPositions()[i]);
	  m_distances[i][j] = Vlength2(dist);

	  // Add distances to the Cd of the second node (if type is BBP)
	  if (pLeaf->getType() == BBP)
	    {
	      VmV(dist, vec2, getPositions()[i]);
	      m_distances[i][size2] = Vlength2(dist);
	    }
	}

      // Add distances to the Cd of the first node (if type is BBP)
      if (getType() == BBP)
	{
	  VmV(dist, cen, vec1);
	  m_distances[size1][j] = Vlength2(dist);
	}
    }

  // Add distance between the Cd of the first and second nodes (both BBPs)
  if (getType() == BBP && pLeaf->getType() == BBP)
    {
      VmV(dist, vec2, vec1);
      m_distances[size1][size2] = Vlength2(dist);
    }
}
/*
void CLeaf::computePairEnergy(CNode * pNode, const REAL rot[3][3],
			      const REAL trans[3], CTerm * term,
			      bool bSeparated)
{
  assert(pNode->isLeaf());
  CLeaf * pLeaf = (CLeaf*) pNode;

  // GLY node does not have any atoms. 
  // There is no interaction with it.
  if (getType() == GLY || pLeaf->getType() == GLY)
    return;

  // If the BVs are too far away, no need to do anything
  if (getBV()->computeDistance(pLeaf->getBV(), rot, trans) > CUTOFF_DISTANCE)
    {
      term->reset();
      return;
    }

  // Fill the pairwise distances matrix.
  computeDistances(pLeaf, rot, trans);

  REAL sum = 0.0;

  int diff = pLeaf->getIndex() - getIndex();

  // Compute all vdW terms.
  sum += CTerm::computeVdW(getType(), pLeaf->getType(),
  			    m_distances, diff);

  // Compute all elctrostatic terms
  sum += CTerm::computeElectrostatics(getType(), pLeaf->getType(),
  			      m_distances, diff);

  // Compute all Solvation terms.
  sum += CTerm::computeSolvation(getType(), pLeaf->getType(),
  				 m_distances, diff);

  // Store the energy sum at the corresponding leaf of the energytree.
  assert(term);
  term->set(sum);

  // Save a pointer to this leaf in case we need to undo the last move.
  m_undoPairs.push_back(term);
  
  return;
}
*/
/*
void CLeaf::computeSelfEnergy(CTerm * term)
{
  // Insert the precomputed energy of interaction between atoms inside
  // this leaf.
  term->set(m_energy);
  m_undoPairs.push_back(term);
  
  return;
}
*/
bool CLeaf::findPairClash(CNode * pNode, const REAL rot[3][3],
			  const REAL trans[3], bool bSeparated)  
{
 assert(pNode->isLeaf());
 CLeaf * pLeaf = (CLeaf*) pNode;

 // GLY node does not have any atoms. 
 // There cannot be a clash with it.
 if (getType() == GLY || pLeaf->getType() == GLY)
   return false;
 
 int size1 = SIDECHAIN::m_aalist[getType()]->m_size;
 int size2 = SIDECHAIN::m_aalist[pLeaf->getType()]->m_size;
  
 REAL cen[3], dist[3];
 
 // Check all pairs of atoms for possible clash
 for(int j = 0; j < size2; j++)
   {
     MxVpV(cen, rot, pLeaf->getPositions()[j], trans); 
     int bb = SIDECHAIN::m_aalist[pLeaf->getType()]->m_aTypes[j];

     for (int i = 0; i < size1; i++)
       {	
	 VmV(dist, cen, getPositions()[i]);
	 int aa = SIDECHAIN::m_aalist[getType()]->m_aTypes[i];
	 
	 // Check exclusion list to see if this pair of atoms 
	 // should not be checked.
	 int ex = isExcluded(pLeaf->getIndex() - getIndex(), getType(), i, 
			     pLeaf->getType(), j);
	 if (ex != EXCLUDED &&
	     isStericClash(aa, bb, Vlength2(dist), ex == PAIR1_4))
	   {
	     return true;
	   }
       }
   }

 return false;
}

// Undo the changes to this leaf caused by the latest move.
void CLeaf::undo()
{
  if (m_nc & CNode::TRANSFORM)
    {
      m_angle = m_undoAngle;
      m_torsionE = m_undoTorsionE;
      McM(m_rotate, m_undoRotate);
    }

  if (m_nc & CNode::BOX)
    changeRotamer(m_undoRotIndex);
}
 
// Rotate the protein around the rotatable bond associated with 
// this leaf.
void CLeaf::rotate(REAL angle)
{ 
  int tmp1, tmp2;
  // Store undo information.
  //cout << "m_angle:  " << m_angle << endl;
  m_undoAngle = m_angle;
  McM(m_undoRotate, m_rotate);
  m_undoTorsionE = m_torsionE;

  // Change the angle and compute a new rotation matrix.
  m_angle += angle;
  //cout << "m_angle:  " << m_angle << endl;
  compute_rotation(m_rotate, getJoint(), m_angle);
  //cout << "m_joint " << getJoint()[0] <<  getJoint()[1] << getJoint()[2] << endl;
  /*
  for (tmp1 = 0; tmp1 < 3; tmp1++)
  {
    for (tmp2 = 0; tmp2 < 3; tmp2++)
    {
      cout << m_rotate[tmp1][tmp2] << endl;
    }
  }
  */
  // Recompute the torsion energy for this torsion angle.
  //m_torsionE = compute_dihedral(getAngle(), getIndex() % 2); //What is getIndex()????????
}

void CLeaf::translate(REAL transl[3], REAL trans[3], REAL trans_prev[3], REAL rota[3][3])
{

  REAL tmp_trans[3];
  REAL tmp_trans2[3];
  REAL tmp_trans3[3];
  REAL tmp_rot[3][3];
  REAL tmp_rot2[3][3];

  VmV(tmp_trans, trans, trans_prev);
  McM(tmp_rot, rota);
  Mqinverse(tmp_rot2, tmp_rot);
  //cout << "tmp_trans " << tmp_trans[0] << endl;
  MxV(tmp_trans3, tmp_rot2, tmp_trans);
  //cout << "tmp_trans3 " << tmp_trans3[0] << endl;
  //VcV(m_translate, tmp_trans3);
  VcV(transl, tmp_trans3);
  //cout << "m_translate " << m_translate[0] << "   " << m_translate[1] << "   " << m_translate[2] << endl;
  //MxV(m_translate, tmp_rot2, tmp_trans);
}
