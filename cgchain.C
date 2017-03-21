#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <functional>
#include <string>
#include <vector>
#include <cfloat>
#include <fstream>
#include <iostream>
#include <cstring>

#include "cgpdb.H"
#include "cgchain.H"
#include "cgleaf.H"
#include "cgpairtree.H"
#include "cgeef1.H"
#include "cgrss.H"
#include "cgspheres.H"

// Create a chain from a file describing the protein's structure.
CChain::CChain(const char * fname, const char * fnameAX, const char * fname3D, const char * fn_quadB, const char * fn_quadE, const char * fn_quadI, const char * fn_tripB, const char * fn_tripE, const char * fn_tripI, const char * fn_sing, const char * fl_quadB, const char * fl_quadE, const char * fl_quadI, const char * fl_tripB, const char * fl_tripE, const char * fl_tripI, const char * fl_sing, int type)
{
  vector<int> AAtypes;
  vector<REAL> phis, psis;
  vector<int> rots;
  //vector<REAL> coorsx, coorsy, coorsz;
  vector<REAL> axis_x, axis_y, axis_z;
  
  if (type == ANGS_FILE)
    load_angs(fname, AAtypes, phis, psis, rots);
  else
    {
      cout << "Bad type of input file!" << endl;
      exit(0);
    }
  load_axes(fnameAX, axis_x, axis_y, axis_z);
  load_quad_buried_prob(fn_quadB, quadB_prob);
  load_quad_exposed_prob(fn_quadE, quadE_prob);
  load_quad_intermediate_prob(fn_quadI, quadI_prob);
  load_trip_buried_prob(fn_tripB, tripB_prob);
  load_trip_exposed_prob(fn_tripE, tripE_prob);
  load_trip_intermediate_prob(fn_tripI, tripI_prob);
  load_sing_prob(fn_sing, sing_prob);

  load_quadlig_buried_prob(fl_quadB, quadligB_prob);
  load_quadlig_exposed_prob(fl_quadE, quadligE_prob);
  load_quadlig_intermediate_prob(fl_quadI, quadligI_prob);
  load_triplig_buried_prob(fl_tripB, tripligB_prob);
  load_triplig_exposed_prob(fl_tripE, tripligE_prob);
  load_triplig_intermediate_prob(fl_tripI, tripligI_prob);
  load_singlig_prob(fl_sing, singlig_prob);

  create(AAtypes, phis, psis, rots, axis_x, axis_y, axis_z);

  //Compute the reference solvent energy and the torsion energy for this chain.
  //computeTorsionEnergy();
  //m_solventE = computeSolventRefEnergy();
}

// Load the chain from a description that specifies phi/psi angles and rotamer indices.
void CChain::load_angs(const char * fname, vector<int> & AAtypes, vector<REAL> & phis, vector<REAL> & psis, vector<int> & rots)
{
  ifstream fin(fname);
  if (!fin.is_open())
    {
      cout << "Could not open input file " << fname << endl;
      exit(0);
    }

  //string buf;
  char buf[100];

  float phi, psi;
  int rotIndex;
  char aa[5];
  while (!fin.eof())
    {
      fin.getline(buf, 99);
      
      if (strlen(buf) == 0)
	continue;
      int res = sscanf(buf, "%s %f %f %d", aa, &phi, &psi, &rotIndex);
      if (res < 1)
	{
	  cout << "File formatting error for .angs file: " << buf << endl;
	  exit(0);
	}

      AAtypes.push_back(getAA(aa));

      if (res >= 2)
       {
	phis.push_back((REAL) (phi*M_PI/180));
       }
      else
	phis.push_back(0.0);

      if (res >= 3)
       {
	psis.push_back((REAL) (psi*M_PI/180));
       }
      else
	psis.push_back(0.0);
      
      if (res == 4)
	{
	  rots.push_back(rotIndex);
	  if (rotIndex >= SIDECHAIN::m_aalist[AAtypes[AAtypes.size()-1]]->m_nRotamers)
	    {
	      cout << "Rotamer index " << rotIndex << " is too big for " << aa 
		   << rots.size() << " which has only " 
		   << SIDECHAIN::m_aalist[AAtypes[AAtypes.size()-1]]->m_nRotamers 
		   << " rotamers (the index is zero-based)!!!" << endl;
	      exit(0);
	    }
	}
      else
	rots.push_back(0);

      int i = rots.size()-1;
    }
  // The first phi and last psi are always 0
  phis[0] = 0.0; psis[psis.size()-1] = 0.0;

  cout << "loaded " << AAtypes.size() << " AAs" << endl;
  fin.close();
}
// Load translation and axes of rotation
void CChain::load_axes(const char * fnameAX, vector<REAL> & axes_x, vector<REAL> & axes_y, vector<REAL> & axes_z)
{
  ifstream fin(fnameAX);
  if (!fin.is_open())
    {
      cout << "Could not open coordinate file " << fnameAX << endl;
      exit(0);
    }

  //string buf;
  char buf[100];

  float cx, cy, cz;
  while (!fin.eof())
    {
       fin.getline(buf, 99);

       if (strlen(buf) == 0)
           continue;
       int res = sscanf(buf, "%f %f %f", &cx, &cy, &cz);
       if (res < 1)
         {
          cout << "File formatting error for .coords file: " << buf << endl;
          exit(0);
         }


        axes_x.push_back(cx);
        axes_y.push_back(cy);
        axes_z.push_back(cz);
      }
   fin.close();
   cout << "loaded  axes" << endl;
}
void CChain::load_quad_buried_prob(const char * fn_quadB, vector<REAL> & quadB_prob)
{
  cout << fn_quadB << endl;
  cout << "Start loading probabilities" << endl;
  ifstream fin(fn_quadB);

  if (!fin.is_open())
    {
      cout << "Could not open quadlig_prob_buried file " << fn_quadB << endl;
      exit(0);
    }

  //cout << "quadB to be started" << endl; 

  char buf[200];

  float c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20;

  while (!fin.eof())
    {
      fin.getline(buf, 199);
      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8, &c9, &c10, &c11, &c12, &c13, &c14, &c15, &c16, &c17, &c18, &c19, &c20);
      quadB_prob.push_back(c1);
      quadB_prob.push_back(c2);
      quadB_prob.push_back(c3);
      quadB_prob.push_back(c4);
      quadB_prob.push_back(c5);
      quadB_prob.push_back(c6);
      quadB_prob.push_back(c7);
      quadB_prob.push_back(c8);
      quadB_prob.push_back(c9);
      quadB_prob.push_back(c10);
      quadB_prob.push_back(c11);
      quadB_prob.push_back(c12);
      quadB_prob.push_back(c13);
      quadB_prob.push_back(c14);
      quadB_prob.push_back(c15);
      quadB_prob.push_back(c16);
      quadB_prob.push_back(c17);
      quadB_prob.push_back(c18);
      quadB_prob.push_back(c19);
      quadB_prob.push_back(c20);
   }
  fin.close();   
 
  cout << "loaded quadB" << endl;
}


void CChain::load_quad_exposed_prob(const char * fn_quadE, vector<REAL> & quadE_prob)
{
  cout << fn_quadE << endl;
  ifstream fin(fn_quadE);
  if (!fin.is_open())
    {
      cout << "Could not open quad_prob_exposed file " << fn_quadE << endl;
      exit(0);
    }

  char buf[200];

  float c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20;

  while (!fin.eof())
    {
      fin.getline(buf, 199);
      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8, &c9, &c10, &c11, &c12, &c13, &c14, &c15, &c16, &c17, &c18, &c19, &c20);

      quadE_prob.push_back(c1);
      quadE_prob.push_back(c2);
      quadE_prob.push_back(c3);
      quadE_prob.push_back(c4);
      quadE_prob.push_back(c5);
      quadE_prob.push_back(c6);
      quadE_prob.push_back(c7);
      quadE_prob.push_back(c8);
      quadE_prob.push_back(c9);
      quadE_prob.push_back(c10);
      quadE_prob.push_back(c11);
      quadE_prob.push_back(c12);
      quadE_prob.push_back(c13);
      quadE_prob.push_back(c14);
      quadE_prob.push_back(c15);
      quadE_prob.push_back(c16);
      quadE_prob.push_back(c17);
      quadE_prob.push_back(c18);
      quadE_prob.push_back(c19);
      quadE_prob.push_back(c20);
   }
  fin.close();

  cout << "loaded quadE" << endl;
}

void CChain::load_quad_intermediate_prob(const char * fn_quadI, vector<REAL> & quadI_prob)
{
  cout << fn_quadI << endl;
  ifstream fin(fn_quadI);
  if (!fin.is_open())
    {
      cout << "Could not open quad_prob_intermediate file " << fn_quadI << endl;
      exit(0);
    }

  char buf[200];

  float c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20;

  while (!fin.eof())
    {
      fin.getline(buf, 199);

      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8, &c9, &c10, &c11, &c12, &c13, &c14, &c15, &c16, &c17, &c18, &c19, &c20);

      quadI_prob.push_back(c1);
      quadI_prob.push_back(c2);
      quadI_prob.push_back(c3);
      quadI_prob.push_back(c4);
      quadI_prob.push_back(c5);
      quadI_prob.push_back(c6);
      quadI_prob.push_back(c7);
      quadI_prob.push_back(c8);
      quadI_prob.push_back(c9);
      quadI_prob.push_back(c10);
      quadI_prob.push_back(c11);
      quadI_prob.push_back(c12);
      quadI_prob.push_back(c13);
      quadI_prob.push_back(c14);
      quadI_prob.push_back(c15);
      quadI_prob.push_back(c16);
      quadI_prob.push_back(c17);
      quadI_prob.push_back(c18);
      quadI_prob.push_back(c19);
      quadI_prob.push_back(c20);
   }
  fin.close();

  cout << "loaded quadI" << endl;
}

void CChain::load_trip_buried_prob(const char * fn_tripB, vector<REAL> & tripB_prob)
{
  ifstream fin(fn_tripB);
  if (!fin.is_open())
    {
      cout << "Could not open trip_prob_buried file " << fn_tripB << endl;
      exit(0);
    }

  char buf[200];

  float c1;

  while (!fin.eof())
    {
      fin.getline(buf, 199);

      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f", &c1);

      tripB_prob.push_back(c1);
    }
  fin.close();

  cout << "loaded tripB" << endl;
}

void CChain::load_trip_exposed_prob(const char * fn_tripE, vector<REAL> & tripE_prob)
{
  ifstream fin(fn_tripE);
  if (!fin.is_open())
    {
      cout << "Could not open trip_prob_exposed file " << fn_tripE << endl;
      exit(0);
    }

  char buf[200];

  float c1;

  while (!fin.eof())
    {
      fin.getline(buf, 199);

      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f", &c1);

      tripE_prob.push_back(c1);
    }
  fin.close();

  cout << "loaded tripE" << endl;
}

void CChain::load_trip_intermediate_prob(const char * fn_tripI, vector<REAL> & tripI_prob)
{
  ifstream fin(fn_tripI);
  if (!fin.is_open())
    {
      cout << "Could not open trip_prob_exposed file " << fn_tripI << endl;
      exit(0);
    }

  char buf[200];

  float c1;

  while (!fin.eof())
    {
      fin.getline(buf, 199);

      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f", &c1);

      tripI_prob.push_back(c1);
    }
  fin.close();

  cout << "loaded tripI" << endl;
}

void CChain::load_sing_prob(const char * fn_sing, vector<REAL> & sing_prob)
{
  ifstream fin(fn_sing);
  if (!fin.is_open())
    {
      cout << "Could not open sing_prob file " << fn_sing << endl;
      exit(0);
    }

  char buf[200];

  float c1;

  while (!fin.eof())
    {
      fin.getline(buf, 199);

      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f", &c1);

      sing_prob.push_back(c1);
    }
  fin.close();

  cout << "loaded sing" << endl;
}

//##################################################################
void CChain::load_quadlig_buried_prob(const char * fl_quadB, vector<REAL> & quadligB_prob)
{

  ifstream fin(fl_quadB);
  if (!fin.is_open())
    {
      cout << "Could not open quadlig_prob_buried file " << fl_quadB << endl;
      exit(0);
    }

//  cout << "quadB to be started" << endl; 

  char buf[200];

  float c1, c2, c3, c4;

  while (!fin.eof())
    {
      fin.getline(buf, 199);
      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f %f %f %f", &c1, &c2, &c3, &c4);
      quadligB_prob.push_back(c1);
      quadligB_prob.push_back(c2);
      quadligB_prob.push_back(c3);
      quadligB_prob.push_back(c4);
   }
  fin.close();   

  cout << "loaded quadligB" << endl;
}


void CChain::load_quadlig_exposed_prob(const char * fl_quadE, vector<REAL> & quadligE_prob)
{
  ifstream fin(fl_quadE);
  if (!fin.is_open())
    {
      cout << "Could not open quad_prob_exposed file " << fl_quadE << endl;
      exit(0);
    }

  char buf[200];

  float c1, c2, c3, c4;

  while (!fin.eof())
    {
      fin.getline(buf, 199);

      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f %f %f %f", &c1, &c2, &c3, &c4);

      quadligE_prob.push_back(c1);
      quadligE_prob.push_back(c2);
      quadligE_prob.push_back(c3);
      quadligE_prob.push_back(c4);
   }
  fin.close();

  cout << "loaded quadligE" << endl;
}

void CChain::load_quadlig_intermediate_prob(const char * fl_quadI, vector<REAL> & quadligI_prob)
{
  ifstream fin(fl_quadI);
  if (!fin.is_open())
    {
      cout << "Could not open quad_prob_intermediate file " << fl_quadI << endl;
      exit(0);
    }

  char buf[200];

  float c1, c2, c3, c4;

  while (!fin.eof())
    {
      fin.getline(buf, 199);

      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f %f %f %f", &c1, &c2, &c3, &c4);

      quadligI_prob.push_back(c1);
      quadligI_prob.push_back(c2);
      quadligI_prob.push_back(c3);
      quadligI_prob.push_back(c4);
   }
  fin.close();

    cout << "loaded quadligI" << endl;
}

void CChain::load_triplig_buried_prob(const char * fl_tripB, vector<REAL> & tripligB_prob)
{
  ifstream fin(fl_tripB);
  if (!fin.is_open())
    {
      cout << "Could not open trip_prob_buried file " << fl_tripB << endl;
      exit(0);
    }

  char buf[200];

  float c1;

  while (!fin.eof())
    {
      fin.getline(buf, 199);

      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f", &c1);

      tripligB_prob.push_back(c1);
    }
  fin.close();

  cout << "loaded tripligB" << endl;
}

void CChain::load_triplig_exposed_prob(const char * fl_tripE, vector<REAL> & tripligE_prob)
{
  ifstream fin(fl_tripE);
  if (!fin.is_open())
    {
      cout << "Could not open trip_prob_exposed file " << fl_tripE << endl;
      exit(0);
    }

  char buf[200];

  float c1;

  while (!fin.eof())
    {
      fin.getline(buf, 199);

      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f", &c1);

      tripligE_prob.push_back(c1);
    }
  fin.close();

  cout << "loaded tripligE" << endl;
}

void CChain::load_triplig_intermediate_prob(const char * fl_tripI, vector<REAL> & tripligI_prob)
{
  ifstream fin(fl_tripI);
  if (!fin.is_open())
    {
      cout << "Could not open trip_prob_exposed file " << fl_tripI << endl;
      exit(0);
    }

  char buf[200];

  float c1;

  while (!fin.eof())
    {
      fin.getline(buf, 199);

      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f", &c1);

      tripligI_prob.push_back(c1);
    }
  fin.close();

  cout << "loaded tripligI" << endl;
}

void CChain::load_singlig_prob(const char * fl_sing, vector<REAL> & singlig_prob)
{
  ifstream fin(fl_sing);
  if (!fin.is_open())
    {
      cout << "Could not open sing_prob file " << fl_sing << endl;
      exit(0);
    }

  char buf[200];

  float c1;

  while (!fin.eof())
    {
      fin.getline(buf, 199);

      if (strlen(buf) == 0)
         continue;
      int res = sscanf(buf, "%f", &c1);

      singlig_prob.push_back(c1);
    }
  fin.close();

  cout << "loaded singlig" << endl;
}

//############################################################

void CChain::load_coors(const char * fname3D, vector<REAL> & coorsx, vector<REAL> & coorsy, vector<REAL> & coorsz)
{
  ifstream fin(fname3D);
  if (!fin.is_open())
    {
      cout << "Could not open coordinate file " << fname3D << endl;
      exit(0);
    }

  //string buf;
  char buf[100];

  float cx, cy, cz;
  char aa[5];
  while (!fin.eof())
    {
      fin.getline(buf, 99);

      if (strlen(buf) == 0)
        continue;
      int res = sscanf(buf, "%f %f %f", &cx, &cy, &cz);
      if (res < 1)
        {
          cout << "File formatting error for .coords file: " << buf << endl;
          exit(0);
        }


      coorsx.push_back(cx);
      coorsy.push_back(cy);
      coorsz.push_back(cz);
     }
  fin.close();
}

// Create the kinematic chain representaion of the protein from internal 
// coordinates (torsion angles).
void CChain::create(const vector<int> AAtypes, const vector<REAL> & phis, const vector<REAL> & psis, const vector<int> & rots, const vector<REAL> & axes_x, const vector<REAL> & axes_y, const vector<REAL> & axes_z)
{
  assert(AAtypes.size() > 1);
  const ROTAMER * rot;
  CLeaf * pLeaf1, * pLeaf2;

  //pLeaf1 = new CLeaf(ROTAMER::UP, 0, &(D_PHI_AXIS[0]), NTR, 0);
  double NTR_PHI_AXIS[3];
  double N_C_NTR_TRANS[3];
  double PSI_AXIS[3];
  double CA_C_TRANS[3];
  double PHI_AXIS[3];
  double C_CA_TRANS[3];

  NTR_PHI_AXIS[0] = (axes_x[1]-axes_x[0])/sqrt((axes_x[1]-axes_x[0])*(axes_x[1]-axes_x[0])+(axes_y[1]-axes_y[0])*(axes_y[1]-axes_y[0])+(axes_z[1]-axes_z[0])*(axes_z[1]-axes_z[0]));
  NTR_PHI_AXIS[1] = (axes_y[1]-axes_y[0])/sqrt((axes_x[1]-axes_x[0])*(axes_x[1]-axes_x[0])+(axes_y[1]-axes_y[0])*(axes_y[1]-axes_y[0])+(axes_z[1]-axes_z[0])*(axes_z[1]-axes_z[0]));
  NTR_PHI_AXIS[2] = (axes_z[1]-axes_z[0])/sqrt((axes_x[1]-axes_x[0])*(axes_x[1]-axes_x[0])+(axes_y[1]-axes_y[0])*(axes_y[1]-axes_y[0])+(axes_z[1]-axes_z[0])*(axes_z[1]-axes_z[0]));

  cout << "axes " << axes_x[1] << " " << axes_x[2] << " " << axes_x[2]-axes_x[1] << endl;
  N_C_NTR_TRANS[0] = axes_x[1]-axes_x[0];
  N_C_NTR_TRANS[1] = axes_y[1]-axes_y[0];
  N_C_NTR_TRANS[2] = axes_z[1]-axes_z[0];
  //N_C_NTR_TRANS[0] = axes_x[0];
  //N_C_NTR_TRANS[1] = axes_y[0];
  //N_C_NTR_TRANS[2] = axes_z[0];

  pLeaf1 = new CLeaf(ROTAMER::UP, 0, NTR_PHI_AXIS, NTR, 0);
  //pLeaf1 = new CLeaf(0, NTR_PHI_AXIS, NTR, 0);
  VcV(pLeaf1->m_translate, N_C_NTR_TRANS); 
  pLeaf1->setPrev(NULL);
  pLeaf1->rotate(phis[0]);
//pLeaf1->m_size = SIDECHAIN::m_aalist[m_links[j]->getType()]->m_size;
  m_links.push_back(pLeaf1);

  cout << "pLeaf1->m_translate" << pLeaf1->m_translate[0] << "   " << pLeaf1->m_translate[1] << "   " << pLeaf1->m_translate[2] << endl;
  
  //cout << NTR_PHI_AXIS[0] << "   " << NTR_PHI_AXIS[1] << "   " << NTR_PHI_AXIS[2] << endl;
  /* 
  for (int tmp1 = 0; tmp1 < 3; tmp1++)
  {
   for (int tmp2 = 0; tmp2 < 3; tmp2++)
   {
    cout << "pLeaf1->m_rotate" << pLeaf1->m_rotate[tmp1][tmp2] << endl;
   }
   //cout << "\n" << endl;
  }
  */
  //exit(0);
  CLeaf * last = pLeaf1;
  cout << AAtypes.size() << endl;
  for (int i = 0; i < AAtypes.size(); i++)
    {
      int bbn_type;
      if (i == AAtypes.size() - 1)
	bbn_type = CTR;
      else if (AAtypes[i+1] == PRO)
	bbn_type = BBP;
      //else if (AAtypes[i+1] == GLY)
	//bbn_type = BBG;
      else
	bbn_type = BBN;
      //exit(0);
      //cout << "i " << "  " << psis[i] << "  " << bbn_type << endl;
 
      if (i%2 == 0)
	{
          //cout << axes_x[2*i+2] << axes_y[2*i+2] << axes_z[2*i+2] << endl;
          //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
          //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
          PSI_AXIS[0] = (axes_x[2*i+2]-axes_x[2*i+1])/sqrt(pow((axes_x[2*i+2]-axes_x[2*i+1]),2)+pow((axes_y[2*i+2]-axes_y[2*i+1]),2)+pow((axes_z[2*i+2]-axes_z[2*i+1]),2));
          PSI_AXIS[1] = (axes_y[2*i+2]-axes_y[2*i+1])/sqrt(pow((axes_x[2*i+2]-axes_x[2*i+1]),2)+pow((axes_y[2*i+2]-axes_y[2*i+1]),2)+pow((axes_z[2*i+2]-axes_z[2*i+1]),2));
          PSI_AXIS[2] = (axes_z[2*i+2]-axes_z[2*i+1])/sqrt(pow((axes_x[2*i+2]-axes_x[2*i+1]),2)+pow((axes_y[2*i+2]-axes_y[2*i+1]),2)+pow((axes_z[2*i+2]-axes_z[2*i+1]),2));
          //cout << "PSI_AXIS[0]    " << PSI_AXIS[0] << " " << PSI_AXIS[1]  << " " << PSI_AXIS[2] << endl;
          CA_C_TRANS[0] = axes_x[2*i+2]-axes_x[2*i+1];
          CA_C_TRANS[1] = axes_y[2*i+2]-axes_y[2*i+1];
          CA_C_TRANS[2] = axes_z[2*i+2]-axes_z[2*i+1];
          if (CA_C_TRANS[0] == 0.0)
              CA_C_TRANS[0] = 0.001;
          if (CA_C_TRANS[1] == 0.0)
              CA_C_TRANS[1] = 0.001;
          if (CA_C_TRANS[2] == 0.0)
              CA_C_TRANS[2] = 0.001;

          pLeaf1 = new CLeaf(ROTAMER::UP, rots[i], PSI_AXIS, AAtypes[i], 2*i+1);
          //cout << "Link_m_bv " << pLeaf1->m_bv->getBVvertices()[1][0] << endl;
          //cout << "Link_m_pos " << pLeaf1->getPositions()[0][0] << endl;
          //exit(0);
	  VcV(pLeaf1->m_translate, CA_C_TRANS);
	  pLeaf1->rotate(psis[i]);
          cout << "pLeaf1->m_translate" << pLeaf1->m_translate[0] << "   " << pLeaf1->m_translate[1] << "   " << pLeaf1->m_translate[2] << endl;
          //cout << "psi_angle   " << psis[i] << endl;
          /*
          for (int tmp1 = 0; tmp1 < 3; tmp1++)
          {
            for (int tmp2 = 0; tmp2 < 3; tmp2++)
            {
               cout << "pLeaf1->m_rotate  " << pLeaf1->m_rotate[tmp1][tmp2] << endl;
            }
          }
          */
          //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
          //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
          if (bbn_type == CTR)
          { 
            pLeaf2 = NULL;
            continue;
          } 
          PHI_AXIS[0] = (axes_x[2*i+3]-axes_x[2*i+2])/sqrt(pow((axes_x[2*i+3]-axes_x[2*i+2]),2)+pow((axes_y[2*i+3]-axes_y[2*i+2]),2)+pow((axes_z[2*i+3]-axes_z[2*i+2]),2));
          PHI_AXIS[1] = (axes_y[2*i+3]-axes_y[2*i+2])/sqrt(pow((axes_x[2*i+3]-axes_x[2*i+2]),2)+pow((axes_y[2*i+3]-axes_y[2*i+2]),2)+pow((axes_z[2*i+3]-axes_z[2*i+2]),2));
          PHI_AXIS[2] = (axes_z[2*i+3]-axes_z[2*i+2])/sqrt(pow((axes_x[2*i+3]-axes_x[2*i+2]),2)+pow((axes_y[2*i+3]-axes_y[2*i+2]),2)+pow((axes_z[2*i+3]-axes_z[2*i+2]),2));
          //cout << "PHI_AXIS[0]    " << PHI_AXIS[0] << " " << PHI_AXIS[1]  << " " << PHI_AXIS[2] << endl;
          C_CA_TRANS[0] = axes_x[2*i+3]-axes_x[2*i+2];
          C_CA_TRANS[1] = axes_y[2*i+3]-axes_y[2*i+2];
          C_CA_TRANS[2] = axes_z[2*i+3]-axes_z[2*i+2];
          if (C_CA_TRANS[0] == 0.0)
              C_CA_TRANS[0] = 0.001;
          if (C_CA_TRANS[1] == 0.0)
              C_CA_TRANS[1] = 0.001;
          if (C_CA_TRANS[2] == 0.0)
              C_CA_TRANS[2] = 0.001;
          /*
          if (bbn_type == CTR)
          {
             pLeaf2 = NULL;
             continue;
          }*/
          pLeaf2 = new CLeaf(ROTAMER::UP, 0, PHI_AXIS, bbn_type, 2*i+2);
	  VcV(pLeaf2->m_translate, C_CA_TRANS);
	  //if (bbn_type == BBP)
	    //pLeaf2->rotate(-M_PI/3.0);
	  //else
	  pLeaf2->rotate(phis[i+1]);
          //cout << "U_PHI_AXIS" << U_PHI_AXIS[0] << "   " << U_PHI_AXIS[1] << "   " << U_PHI_AXIS[2] << endl;
          cout << "pLeaf2->m_translate" << pLeaf2->m_translate[0] << "   " << pLeaf2->m_translate[1] << "   " << pLeaf2->m_translate[2] << endl;
          /*
          for (int tmp1 = 0; tmp1 < 3; tmp1++)
          {
            for (int tmp2 = 0; tmp2 < 3; tmp2++)
            {
               cout << "pLeaf2->m_rotate" << pLeaf2->m_rotate[tmp1][tmp2] << endl;
            }
          }
          */
          //exit(0);
	}
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
      else
	{
          PSI_AXIS[0] = (axes_x[2*i+2]-axes_x[2*i+1])/sqrt(pow((axes_x[2*i+2]-axes_x[2*i+1]),2)+pow((axes_y[2*i+2]-axes_y[2*i+1]),2)+pow((axes_z[2*i+2]-axes_z[2*i+1]),2));
          PSI_AXIS[1] = (axes_y[2*i+2]-axes_y[2*i+1])/sqrt(pow((axes_x[2*i+2]-axes_x[2*i+1]),2)+pow((axes_y[2*i+2]-axes_y[2*i+1]),2)+pow((axes_z[2*i+2]-axes_z[2*i+1]),2));
          PSI_AXIS[2] = (axes_z[2*i+2]-axes_z[2*i+1])/sqrt(pow((axes_x[2*i+2]-axes_x[2*i+1]),2)+pow((axes_y[2*i+2]-axes_y[2*i+1]),2)+pow((axes_z[2*i+2]-axes_z[2*i+1]),2));
          //cout << "PSI_AXIS[0]    " << PSI_AXIS[0] << " " << PSI_AXIS[1]  << " " << PSI_AXIS[2] << endl;
          CA_C_TRANS[0] = axes_x[2*i+2]-axes_x[2*i+1];
          CA_C_TRANS[1] = axes_y[2*i+2]-axes_y[2*i+1];
          CA_C_TRANS[2] = axes_z[2*i+2]-axes_z[2*i+1];
          if (CA_C_TRANS[0] == 0.0)
              CA_C_TRANS[0] = 0.001;
          if (CA_C_TRANS[1] == 0.0)
              CA_C_TRANS[1] = 0.001;
          if (CA_C_TRANS[2] == 0.0)
              CA_C_TRANS[2] = 0.001;

          pLeaf1 = new CLeaf(ROTAMER::DOWN, rots[i], PSI_AXIS, AAtypes[i], 2*i+1);
	  VcV(pLeaf1->m_translate, CA_C_TRANS);
	  pLeaf1->rotate(psis[i]);
          //cout << "psi_angle   " << psis[i] << endl;
          cout << "pLeaf1->m_translate" << pLeaf1->m_translate[0] << "   " << pLeaf1->m_translate[1] << "   " << pLeaf1->m_translate[2] << endl;
          //cout << "U_PSI_AXIS" << U_PSI_AXIS[0] << "   " << U_PSI_AXIS[1] << "   " << U_PSI_AXIS[2] << endl;
          /*
          for (int tmp1 = 0; tmp1 < 3; tmp1++)
          {
            for (int tmp2 = 0; tmp2 < 3; tmp2++)
            {
               cout << "pLeaf1->m_rotate" << pLeaf1->m_rotate[tmp1][tmp2] << endl;
            }
          }
          */
          //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
          //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
          
          //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
          //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
          if (bbn_type == CTR)
           {
             pLeaf2 = NULL;
             continue;
           }
          PHI_AXIS[0] = (axes_x[2*i+3]-axes_x[2*i+2])/sqrt(pow((axes_x[2*i+3]-axes_x[2*i+2]),2)+pow((axes_y[2*i+3]-axes_y[2*i+2]),2)+pow((axes_z[2*i+3]-axes_z[2*i+2]),2));
          PHI_AXIS[1] = (axes_y[2*i+3]-axes_y[2*i+2])/sqrt(pow((axes_x[2*i+3]-axes_x[2*i+2]),2)+pow((axes_y[2*i+3]-axes_y[2*i+2]),2)+pow((axes_z[2*i+3]-axes_z[2*i+2]),2));
          PHI_AXIS[2] = (axes_z[2*i+3]-axes_z[2*i+2])/sqrt(pow((axes_x[2*i+3]-axes_x[2*i+2]),2)+pow((axes_y[2*i+3]-axes_y[2*i+2]),2)+pow((axes_z[2*i+3]-axes_z[2*i+2]),2));
          //cout << "PHI_AXIS[0]:  " << (axes_x[2*i+2]-axes_x[2*i+1])/sqrt((axes_x[2*i+2]-axes_x[2*i+1])*(axes_x[2*i+2]-axes_x[2*i+1])+(axes_y[2*i+2]-axes_y[2*i+1])*(axes_y[2*i+2]-axes_y[2*i+1])+(axes_z[2*i+2]-axes_z[2*i+1])*(axes_z[2*i+2]-axes_z[2*i+1])) << endl;
          //cout << "phi   " << phis[i+1] << endl;
          C_CA_TRANS[0] = axes_x[2*i+3]-axes_x[2*i+2];
          C_CA_TRANS[1] = axes_y[2*i+3]-axes_y[2*i+2];
          C_CA_TRANS[2] = axes_z[2*i+3]-axes_z[2*i+2];
          if (C_CA_TRANS[0] == 0.0)
              C_CA_TRANS[0] = 0.001;
          if (C_CA_TRANS[1] == 0.0)
              C_CA_TRANS[1] = 0.001;
          if (C_CA_TRANS[2] == 0.0)
              C_CA_TRANS[2] = 0.001;
          /*
          if (bbn_type == CTR)
           {
             pLeaf2 = NULL;
             continue;
           }*/
          pLeaf2 = new CLeaf(ROTAMER::DOWN, 0, PHI_AXIS, bbn_type, 2*i+2);
	  VcV(pLeaf2->m_translate, C_CA_TRANS);
	  pLeaf2->rotate(phis[i+1]);
          //cout << "phi_angle   " << phis[i+1] << endl;
          cout << "pLeaf2->m_translate" << pLeaf2->m_translate[0] << "   " << pLeaf2->m_translate[1] << "   " << pLeaf2->m_translate[2] << endl;
          //cout << "D_PHI_AXIS" << D_PHI_AXIS[0] << "   " << D_PHI_AXIS[1] << "   " << D_PHI_AXIS[2] << endl;
          /*
          for (int tmp1 = 0; tmp1 < 3; tmp1++)
          {
            for (int tmp2 = 0; tmp2 < 3; tmp2++)
            {
               cout << "pLeaf2->m_rotate" << pLeaf2->m_rotate[tmp1][tmp2] << endl;
            }
          }
          */
          //exit(0);
          //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
          //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
          
	}
      
      m_links.push_back(pLeaf1);
      m_links.push_back(pLeaf2);

      last->setNext(pLeaf1);
      pLeaf1->setPrev(last);
      pLeaf1->setNext(pLeaf2);
      pLeaf2->setPrev(pLeaf1);
      last = pLeaf2;
  }
  //last->setNext(NULL);
  /*
  REAL trans_prev[3];
  REAL trans[3];
  REAL rota[3][3];
  REAL Rtemp[3][3];

  Midentity(rota);
  CLeaf * pLeaf = m_links[0];
  McM(rota, pLeaf->m_rotate);

  trans_prev[0] = coorsx[0];
  trans_prev[1] = coorsy[0];
  trans_prev[2] = coorsz[0];
  trans[0] = coorsx[1];
  trans[1] = coorsy[1];
  trans[2] = coorsz[1];
  REAL transl[3];
  VmV(pLeaf->m_translate, trans, trans_prev);

  for (int i = 1; i < getLength()-1; i++)
  {
    CLeaf * pLeaf = m_links[i];

    trans[0] = coorsx[i+1];
    trans[1] = coorsy[i+1];
    trans[2] = coorsz[i+1];
    trans_prev[0] = coorsx[i];
    trans_prev[1] = coorsy[i];
    trans_prev[2] = coorsz[i];
    //exit(0);
    pLeaf->translate(transl, trans, trans_prev, rota);
    pLeaf->m_translate[0] = transl[0];
    pLeaf->m_translate[1] = transl[1];
    pLeaf->m_translate[2] = transl[2];
    //cout << "********" << pLeaf->m_translate[0] << "   " << pLeaf->m_translate[1] << "   " << pLeaf->m_translate[2] << endl;
    //exit(0);
    McM(Rtemp, rota);
    MxM(rota, Rtemp, pLeaf->m_rotate);
  }
  */
  //exit(0);
  m_frames.resize(getLength());
  //computeTorsionEnergy();
}
/*
REAL CChain::computeSolventRefEnergy()
{
  REAL sum = 0.0;

 for (int i = 0; i < getLength(); i++)
   {
     CLeaf * pLeaf = m_links[i];
     int type = pLeaf->getType();
     if (type == GLY)
	    continue;
     
     int size = SIDECHAIN::m_aalist[type]->m_size;
     for (int j = 0; j < size; j++)
       sum += deltaG_ref[SIDECHAIN::m_aalist[type]->m_aTypes[j]];
    }

 return sum;
}
*/
void CChain::readjustLinkCoords()
{
  vector<FRAME>::iterator it = m_frames.begin();
  //Videntity(it->first.v);
  Midentity(it->second.m);

  int l = 0;
  REAL rot[3][3];
  REAL Rtemp[3][3];
  REAL trans[3] = {0.291,0.005,1.427};
  //REAL tmp_trans[3] = {-0.023, 0.031, -0.008};
  const CLeaf * pLeaf = getLink(0);

  McM(rot, pLeaf->m_rotate);
  VcV(it->first.v, trans);
  Vprint(trans);
  while(pLeaf->getNext() != NULL)
    {
      MxVpV(trans, rot, pLeaf->m_translate, trans);
      //McM(Rtemp, rot);
      //MxM(rot, Rtemp,  pLeaf->m_rotate);
      pLeaf = pLeaf->getNext();
      it++;
      VcV(it->first.v, trans);
      McM(it->second.m, rot); // Always identity matrix
      Vprint(trans);
      //pLeaf = pLeaf->getNext();
      //MxVpV(trans, rot, pLeaf->m_translate, trans);
      //VcV(trans, pLeaf->m_translate);
      //McM(Rtemp, rot);
      //MxM(rot, Rtemp,  pLeaf->m_rotate);
      //if (l == 4)
      //   exit(0);
      //l++;
    }
  //exit(0);
  /*
  McM(rot, pLeaf->m_rotate);
  //VcV(trans, pLeaf->m_translate);
  VcV(trans, tmp_trans);
  while(pLeaf->getNext() != NULL)
    {
      it++;

      pLeaf = pLeaf->getNext();
      MxVpV(trans, rot, pLeaf->m_translate, trans);
      VcV(it->first.v, trans);
      McM(it->second.m, rot);
      //MxM(rot, Rtemp,  pLeaf->m_rotate);
      //if (l == 4)
      //   exit(0);
      //l++;     
    }
  //exit(0);
  */
}

void CChain::refillLinkCoords()
{
  vector<FRAME>::iterator it = m_frames.begin();
  REAL Rot_inv[3][3];
  REAL Rtemp[3][3];

  const CLeaf * pLeaf = getLink(0);
  McM(Rtemp, it->second.m);
  cout << " m_rotate3 " << Rtemp[0][0] << " " << Rtemp[0][1] << " " << Rtemp[0][2] << endl;
  while(pLeaf->getNext() != NULL)
  {
     it++;
     cout << " m_rotate3 " << Rtemp[0][0] << " " << Rtemp[0][1] << " " << Rtemp[0][2] << endl;
     pLeaf = pLeaf->getNext();
     Minverse(Rot_inv, Rtemp);
     MxM(pLeaf->m_rotate, Rot_inv, it->second.m);
     McM(Rtemp, it->second.m);
  }
  
}

void CChain::computeLinkCoords()
{
  vector<FRAME>::iterator it = m_frames.begin();
  Videntity(it->first.v);
  Midentity(it->second.m);
  cout << "it->first.v " << it->first.v[0] << " " << it->first.v[1] << " " << it->first.v[2] << endl;
  cout << "it->second.m " << it->second.m[0][0] << " " << it->second.m[0][1] << " " << it->second.m[0][2] << " " << it->second.m[1][0] << " " << it->second.m[1][1] << " " << it->second.m[1][2] << " " << it->second.m[2][0] << " " << it->second.m[2][1] << " " << it->second.m[2][2] << endl;

  REAL rot[3][3];
  REAL Rtemp[3][3];
  REAL trans[3] = {0.291, 0.005, 1.427};

  const CLeaf * pLeaf = getLink(0);

  McM(rot, pLeaf->m_rotate);
  VcV(trans, pLeaf->m_translate);

  while(pLeaf->getNext() != NULL)
    {
      it++;
      VcV(it->first.v, trans);
      McM(it->second.m, rot);

      pLeaf = pLeaf->getNext();

      MxVpV(trans, rot, pLeaf->m_translate, trans);

      McM(Rtemp, rot);
      MxM(rot, Rtemp,  pLeaf->m_rotate);
      
  Vprint(trans);cout << "it->first.v " << it->first.v[0] << " " << it->first.v[1] << " " << it->first.v[2] << endl;
      cout << "it->second.m " << it->second.m[0][0] << " " << it->second.m[0][1] << " " << it->second.m[0][2] << " " << it->second.m[1][0] << " " << it->second.m[1][1] << " " << it->second.m[1][2] << " " << it->second.m[2][0] << " " << it->second.m[2][1] << " " << it->second.m[2][2] << endl;
    }
}

void CChain::storeAngsStyle(const char * fname)
{
  ofstream fout(fname);

  if (!fout.is_open())
    {
      cout << " Could not open file " << fname << " for output!!!" << endl;
      return;
    }
  
  for (int j = 1; j < getLength(); j += 2)
  {
    fout << AA_NAMES[m_links[j]->getType()] << " ";
    REAL aa = m_links[j-1]->getAngle() * 180/M_PI;
    while (aa > 180.0)
      aa -= 360.0;
    while (aa < -180)
      aa += 360;

    fout << aa << " ";
    aa = m_links[j]->getAngle() * 180/M_PI;
    while (aa > 180.0)
      aa -= 360.0;
    while (aa < -180)
      aa += 360;

    fout << aa << " ";

    fout << m_links[j]->getRotIndex() << endl;
  }
}


void CChain::adjustBV()
{

  CRss * tmp;
  int size = 0;
  CLeaf * pLeaf = getLink(0);
  while (pLeaf->getNext() != NULL)
  {
     size = SIDECHAIN::m_aalist[pLeaf->getType()]->m_size;
     //pLeaf->m_bv = new CRss(pLeaf->getPositions(),size);
     //exit(0);
     tmp = new CRss(pLeaf->getPositions(),size);
     delete pLeaf->m_bv;
     pLeaf->m_bv = tmp;
     pLeaf = pLeaf->getNext();
//     delete pLeaf->m_bv;
//     pLeaf->m_bv->getBVvertices() = CRss->computeRss(pLeaf->getPositions(),size);
//     pLeaf->m_bv = (CRss*) new CRss(pLeaf->getPositions(),size);
//     exit(0);
  }
  /*
  SIDECHAIN * pSC;
  for (int aa = 0; aa < NUMAA; aa++)
    {
      SIDECHAIN * pSC = SIDECHAIN::m_aalist[aa];
      cout << "NUM_ROTAMERS:" << NUM_ROTAMERS[aa] << "  " << pSC->m_size << endl;
      for (int i = 0; i < NUM_ROTAMERS[aa]; i++)
      {
        ROTAMER * ss = new ROTAMER(pSC->m_size);
//#ifdef USE_RSS
        pSC->m_Urotamers[i]->m_bv = new CRss(ss->m_positions,  pSC->m_size);
//#else
        pSC->m_Urotamers[i]->m_bv = new CSphere(ss->m_positions,  pSC->m_size);
//#endif
        ROTAMER * pp = new ROTAMER(pSC->m_size);
//#ifdef USE_RSS
        pSC->m_Drotamers[i]->m_bv = new CRss(pp->m_positions,  pSC->m_size);
//#else
//        pSC->m_Drotamers[i]->m_bv = new CSphere(pp->m_positions,  pSC->m_size);
//#endif
        if (i == 0)
        {
             cout << pSC->m_Urotamers[i]->m_bv << endl;
             cout << pSC->m_Drotamers[i]->m_bv << endl;
       //    pSC->m_Drotamers[i]->m_bv->computeVolume();
           //cout << pSC->m_Drotamers[i]->m_bv->getVolume() << endl;
        }
        delete ss;
        delete pp;
        
      }
    }*/
}


// Store the chain in a PDB style file.
void CChain::store_orig_Coordinates(const char * fname, const char * fname3D)
{
  cout << "Entered storeCoordinates " << endl;
  readjustLinkCoords();
  //computeLinkCoords();
  vector<REAL> coorsx;
  vector<REAL> coorsy;
  vector<REAL> coorsz;
  load_coors(fname3D, coorsx, coorsy, coorsz);

  vector<AA> aas;
  int presize = 0;
  double det;
  AA aa;
  ATOM_ at;
  vector<FRAME>::const_iterator it = getFrames().begin();
  int i = 0, j = 0;
  int tmp1, tmp2;
  double vec[3], vec_m[3], tmp_trans[3], tmp[3];
  REAL Rot_inv[3][3];
  double tmp_mat[3][3],Rtemp[3][3];
  double t1, t2, t3;
  CLeaf * pLeaf = getLink(0);

  for (j = 0; j < getLength(); j++)
    {
      det = 0.0;
      int size = SIDECHAIN::m_aalist[m_links[j]->getType()]->m_size;
      cout << j << " Size " << size << " " << "presize " << presize << endl;

      for (int k = 0; k < size; k++)
       {

	  if (j == 0) 
	    {
	      if (k == 0)
              {
		aa.type = m_links[j+1]->getType();  //struct AA {vector<ATOM_> atoms;int type;} aa
              }
              vec[0] = coorsx[presize+k];
              vec[1] = coorsy[presize+k];
              vec[2] = coorsz[presize+k];
	    }          
          else
            {
              if (j % 2 == 0 && j != getLength() -1 && k == 2)
     	          {
		    aas.push_back(aa);   // vector<AA> aas; vector of struct AA
		    aa.atoms.clear();   // aa.atoms has name and pos
                
	            aa.type = m_links[j+1]->getType();
	          }
            }


          if (j > 0)
            {
               vec[0] = coorsx[presize+k];
               vec[1] = coorsy[presize+k];
               vec[2] = coorsz[presize+k];
               cout << "k " << k << " " <<  vec[0] << " " << vec[1] << " " << vec[2] << endl;
               

               if (k == 0)
                {
                   vec_m[0] = coorsx[presize-SIDECHAIN::m_aalist[m_links[j-1]->getType()]->m_size];
                   vec_m[1] = coorsy[presize-SIDECHAIN::m_aalist[m_links[j-1]->getType()]->m_size];
                   vec_m[2] = coorsz[presize-SIDECHAIN::m_aalist[m_links[j-1]->getType()]->m_size];

                   
                   tmp_mat[0][0] = (vec[0]-vec_m[0])/(pLeaf->m_translate[0]*3.0);
                   tmp_mat[0][1] = (vec[0]-vec_m[0])/(pLeaf->m_translate[1]*3.0);
                   tmp_mat[0][2] = (vec[0]-vec_m[0])/(pLeaf->m_translate[2]*3.0);
                   tmp_mat[1][0] = (vec[1]-vec_m[1])/(pLeaf->m_translate[0]*3.0);
                   tmp_mat[1][1] = (vec[1]-vec_m[1])/(pLeaf->m_translate[1]*3.0);
                   tmp_mat[1][2] = (vec[1]-vec_m[1])/(pLeaf->m_translate[2]*3.0);
                   tmp_mat[2][0] = (vec[2]-vec_m[2])/(pLeaf->m_translate[0]*3.0);
                   tmp_mat[2][1] = (vec[2]-vec_m[2])/(pLeaf->m_translate[1]*3.0);
                   tmp_mat[2][2] = (vec[2]-vec_m[2])/(pLeaf->m_translate[2]*3.0);
                   
                   //McM(tmp_mat, it->second.m);
                   //cout << tmp_mat[0][0]*it->first.v[0] + tmp_mat[0][1]*it->first.v[1] + tmp_mat[0][2]*it->first.v[2] + vec_m[0] << endl;
                   cout << tmp_mat[0][0]*pLeaf->m_translate[0] + tmp_mat[0][1]*pLeaf->m_translate[1] + tmp_mat[0][2]*pLeaf->m_translate[2] + vec_m[0] << " " << tmp_mat[1][0]*pLeaf->m_translate[0] + tmp_mat[1][1]*pLeaf->m_translate[1] + tmp_mat[1][2]*pLeaf->m_translate[2] + vec_m[1] << " " << tmp_mat[2][0]*pLeaf->m_translate[0] + tmp_mat[2][1]*pLeaf->m_translate[1] + tmp_mat[2][2]*pLeaf->m_translate[2] + vec_m[2] << endl;
                   cout << "pLeaf_mtranslate " << pLeaf->m_translate[0] << " " << pLeaf->m_translate[1] << " " << pLeaf->m_translate[2] << endl;
                   cout << vec_m[0] << " " << vec_m[1] << " " << vec_m[2] << endl;

                   VcV(it->first.v, vec);
                   McM(it->second.m, tmp_mat);    
                   Minverse(Rot_inv, it->second.m);
                   //VcV(it->first.v, vec);
                   //cout << SIDECHAIN::m_aalist[m_links[j]->getType()]->m_aNames[k] << " " << vec[0] << " " << vec[1] << " " << vec[2] << endl;
                   cout << "tmp_mat " << tmp_mat[0][0] << " " << tmp_mat[0][1] << " " << tmp_mat[0][2] << endl;
                   cout << "tmp_mat " << tmp_mat[1][0] << " " << tmp_mat[1][1] << " " << tmp_mat[1][2] << endl;
                   cout << "tmp_mat " << tmp_mat[2][0] << " " << tmp_mat[2][1] << " " << tmp_mat[2][2] << endl;
                   //cout << "it->second.m " << Rot_inv[0][0] << " " << Rot_inv[0][1] << " " << Rot_inv[0][2] << endl;
                   //cout << "it->second.m " << Rot_inv[1][0] << " " << Rot_inv[1][1] << " " << Rot_inv[1][2] << endl;
                   //cout << "it->second.m " << Rot_inv[2][0] << " " << Rot_inv[2][1] << " " << Rot_inv[2][2] << endl;
                   pLeaf = pLeaf->getNext();
                 }
               else
                 {
                    vec_m[0] = coorsx[presize];
                    vec_m[1] = coorsy[presize];
                    vec_m[2] = coorsz[presize];

                    tmp_trans[0] = vec[0]-vec_m[0];
                    tmp_trans[1] = vec[1]-vec_m[1];
                    tmp_trans[2] = vec[2]-vec_m[2];

                    MxV_REAL(tmp, Rot_inv, tmp_trans);
                    VcV(m_links[j]->getPositions()[k], tmp);
                    MxVpV(vec, it->second.m, m_links[j]->getPositions()[k], it->first.v);  
                    cout << SIDECHAIN::m_aalist[m_links[j]->getType()]->m_aNames[k] << " " << vec[0] << " " << vec[1] << " " << vec[2] << endl;
                    //cout << "tmp_trans " << tmp_trans[0] << " " << tmp_trans[1] << " " << tmp_trans[2] << endl;
                    //cout << "Rot_inverse " << Rot_inv[0][0] << " " << Rot_inv[0][1] << " " << Rot_inv[0][2] << endl;
                    //cout << "tmp " << tmp[0] << " " << tmp[1] << " " << tmp[2] << endl;
                    cout << "getPosition " << m_links[j]->getPositions()[k][0] << " " << m_links[j]->getPositions()[k][1] << " " << m_links[j]->getPositions()[k][2] << endl;
                    //exit(0);
                 }
            }

            if (m_links[j]->getType() == 18)
            {
              m_links[j]->getCenter()[0] = it->first.v[0];
              m_links[j]->getCenter()[1] = it->first.v[1];
              m_links[j]->getCenter()[2] = it->first.v[2];
            }
            else
            {
              m_links[j]->getCenter()[0] += m_links[j]->getPositions()[k][0];
              m_links[j]->getCenter()[1] += m_links[j]->getPositions()[k][1];
              m_links[j]->getCenter()[2] += m_links[j]->getPositions()[k][2];
            }



	    strcpy(at.name, SIDECHAIN::m_aalist[m_links[j]->getType()]->m_aNames[k]); //at.name
	    at.pos[0] = vec[0];
	    at.pos[1] = vec[1];
	    at.pos[2] = vec[2];

	    aa.atoms.push_back(at);
       }

       if (m_links[j]->getType() != 18)
       {
           m_links[j]->getCenter()[0] /= size;
           m_links[j]->getCenter()[1] /= size;
           m_links[j]->getCenter()[2] /= size;
       }
       presize += size;
      
       it++;
    }
  aas.push_back(aa);
  //writeToMOL2("protein.mol2", aas);
  writeToPDB(fname, aas);
  //refillLinkCoords();
  exit(0);
}

// Store the chain in a PDB style file.
void CChain::storeCoordinates(const char * fname)
{
  cout << "Entered storeCoordinates " << endl;
  computeLinkCoords();
  vector<AA> aas;
  AA aa;
  ATOM_ at;
  vector<FRAME>::const_iterator it = getFrames().begin();
  int i = 0, j = 0;
  const CLeaf * pLeaf = getLink(0);

  for (j = 0; j < getLength(); j++)
    {
      int size = SIDECHAIN::m_aalist[m_links[j]->getType()]->m_size;
      for (int k = 0; k < size; k++)
	{	
	  double vec[3];
	  if (j == 0) 
	    {
	      if (k == 0)
		aa.type = m_links[j+1]->getType();  
	    }
	  else
            {
              if (j % 2 == 0 && j != getLength() -1 && k == 2)
	      {
		aas.push_back(aa);
		aa.atoms.clear();
		aa.type = m_links[j+1]->getType();
	      }
             }
          if (k == 0)
            { 
               MxVpV(vec, it->second.m, pLeaf->m_translate, it->first.v);
            }
          else
               MxVpV(vec, it->second.m, m_links[j]->getPositions()[k], it->first.v);
          
          if (k == 1)
          { 
            cout << j << " m_translate2 " << it->first.v[0] << " " << it->first.v[1] << " " << it->first.v[2] << endl;
             cout << j << " m_rotate2 " << it->second.m[0][0] << " " << it->second.m[0][1] << " " << it->second.m[0][2] << endl;
          }
          
	  strcpy(at.name, SIDECHAIN::m_aalist[m_links[j]->getType()]->m_aNames[k]); //at.name 

	  at.pos[0] = vec[0];
	  at.pos[1] = vec[1];
	  at.pos[2] = vec[2];
	  aa.atoms.push_back(at);
	}
      
      it++;
      pLeaf = pLeaf->getNext();
    }
  aas.push_back(aa);
  //writeToMOL2("protein.mol2", aas);
  writeToPDB(fname, aas);
  exit(0);
}
/*
// Compute the torsion energy, which is a function of the backbone angles
REAL CChain::computeTorsionEnergy()
{
  m_torsionE = 0.0;
  for (int i = 1; i < getLength() - 2; i++)
    m_torsionE += m_links[i]->getTorsionE();

  return m_torsionE;
}

// Update the torsion energy term after the gicen change is applied
// to the backbone angles.
REAL CChain::updateTorsionEnergy(const vector<ANGLE_CHANGE> & angles)
{
  m_undoTorsionE = m_torsionE;

  for (int i = 0; i < ((int)angles.size())-1; i++)
    m_torsionE += m_links[angles[i].m_index]->getTorsionChange();
    
  return m_torsionE;
}
*/
//Compute the absolute positions of all Ca atoms in the protein chain.
void CChain::getCaPositions(POSITIONS & pos)
{
  int length = getLength() / 2;
  REAL Calphas[length][3];
  
  computeCaPositions(Calphas);

  if (pos.size() != length)
    pos.resize(length, vector<REAL>(3));

  for (int i = 0; i < length; i++)
    for (int j = 0; j < 3; j++)
      pos[i][j] = Calphas[i][j];
}

void CChain::computeCaPositions(REAL Calphas[][3])
{
  REAL rot[3][3];
  REAL Rtemp[3][3];
  REAL trans[3];
 
  const CLeaf * pLeaf = getLink(0);
  VcV(Calphas[0], pLeaf->getPositions()[4]);

  McM(rot, pLeaf->m_rotate);
  VcV(trans, pLeaf->m_translate);

  int i = 0;
  while(pLeaf->getNext()->getNext() != NULL)
    {
      pLeaf = pLeaf->getNext();

      if (pLeaf->getIndex() % 2 == 0)
	{
	  i++;
	  MxVpV(Calphas[i], rot, 
		pLeaf->getPositions()[4], trans);
	}

      MxVpV(trans, rot, pLeaf->m_translate, trans);      
      McM(Rtemp, rot);
      MxM(rot, Rtemp,  pLeaf->m_rotate);
    }
}

