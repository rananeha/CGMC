#ifndef _CGLEAF_H
#define _CGLEAF_H

#include "cgeef1.H"
#include "bv.H"
#include "cgnode.H"
#include "cgrss.H"
// A structure describing one angle change.
struct ANGLE_CHANGE {
  int m_index;
  REAL m_angle;
};

struct ANGLE_CHANGE_COMP {
  bool operator() (const ANGLE_CHANGE & aa, const ANGLE_CHANGE & bb)
  {
    return (aa.m_index < bb.m_index);
  }
};
struct CLOSE_NEIGHBOR {
  int m_index;
};

// A structure describing one rotamer change.
struct ROTAMER_CHANGE {
  int m_index;
  int m_rotIndex;
};

struct ROTAMER_CHANGE_COMP {
  bool operator() (const ROTAMER_CHANGE & aa, const ROTAMER_CHANGE & bb)
  {
    return (aa.m_index < bb.m_index);
  }
};

struct LIGAND_CHANGE {
  int m_index;
  int m_confIndex;
};
 
//class CTerm;

class CLeaf : public CNode {
public:
  // Create a leaf with the given rotamer value and type.
  CLeaf(ROTAMER::ROT_TYPE rType, int rotIndex, REAL * joint, 
//    CLeaf(int rotIndex, REAL * joint,
	int type, int index) : 	CNode(index), m_joint(joint), m_type(type), 
				m_rotType(rType), m_angle(0.0), m_energy(0.0), 
				//m_angle(0.0), m_energy(0.0),
				m_torsionE(0.0), m_size(0) 
  { changeRotamer(rotIndex); }

  // Create a leaf given a set of coordinates for the atoms.
  CLeaf(int index, REAL (*pos)[3], int type) : CNode(index), m_positions(pos),
					   m_type(type), m_angle(0.0), m_energy(0.0), 
					   m_torsionE(0.0), m_joint(NULL) {}
  
  void rotate(REAL angle);
  void translate(REAL transl[3], REAL trans[3], REAL trans_prev[3], REAL rota[3][3]);
  //virtual void computePairEnergy(CNode * pNode, const REAL rot[3][3],
	//		      const REAL trans[3], CTerm * term,
	//		      bool bSeparated);
  //virtual void computeSelfEnergy(CTerm * term);
  virtual bool findPairClash(CNode * pNode, const REAL rot[3][3],
			       const REAL trans[3], bool bSeparated);

  // A leaf is never in clash with itself.
  virtual bool findSelfClash()
  { return false; }
  virtual void undo();

  void changeRotamer(int rotIndex);

  // Accessors
  int getType() const
  { return m_type; }
  int getRotIndex() const
  { return m_rotIndex; }
  virtual bool isLeaf() const
  { return true; }
  CLeaf * getNext() const
  { return (CLeaf*) m_next; }
  
  COORDS  getPositions() const 
  { return m_positions; }
  const REAL * getJoint() const
  { return m_joint; }
  REAL getAngle() const
  { return m_angle; }
  //REAL getTorsionE() const
  //{ return m_torsionE; }
  //REAL getTorsionChange() const
  //{ return m_torsionE - m_undoTorsionE; }
  double * getCenter()
  { return m_sc; }
  int getLeafSize() const
  { return m_size; }

  //static REAL m_distances[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE];
  static REAL m_distances[13][13];

private:
  void computeDistances(CLeaf * pLeaf, const REAL rot[3][3],
			const REAL trans[3]);

  const REAL * m_joint;
  REAL m_angle;

  int m_type;

  REAL (*m_positions)[3];
  int m_rotIndex;
  ROTAMER::ROT_TYPE m_rotType;
  int m_nPos;
  REAL m_energy;
  REAL m_torsionE;
  REAL m_prob;
  REAL m_chi1mean;
  REAL m_chi2mean;
  REAL m_chi3mean;
  REAL m_chi4mean;
  REAL m_chi1std;
  REAL m_chi2std;
  REAL m_chi3std;
  REAL m_chi4std;
  double m_sc[3];
  int m_size;
  // undo variables
  REAL m_undoAngle;
  REAL m_undoTorsionE;
  int m_undoRotIndex;
};


#endif
