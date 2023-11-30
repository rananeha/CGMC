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

#include "cgchain.H"
#include "cgconformer.H"
//#include "cgeef1.H"

//int SINGLET_IDX[TOT_SING] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
int SINGLET_IDX[TOT_SING] = {19,0,1,2,3,4,5,18,6,7,8,9,10,11,12,13,14,15,16,17};
int TRIPLET_IDX[TOT_TRIP] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119};

//char *SINGLETS[TOT_SINGLETS] = {"ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "GLY", "ALA"};

//char *TRIPLETS[TOT_TRIPLETS] = {"BBB","BAA","BBA","BAH","BAO","BAS","BAC","BAP","BAN","BHH","BBH","BHO","BHS","BHC","BHP","BHN","BOO","BBO","BOS","BOC","BOP","BON","BSS","BBS","BSC","BSP","BSN","BCC","BBC","BCP","BCN","BPP","BBP","BPN","BNN","BBN","AAA","AHH","AAH","AHO","AHS","AHC","AHP","AHN","AOO","AAO","AOS","AOC","AOP","AON","ASS","AAS","ASC","ASP","ASN","ACC","AAC","ACP","ACN","APP","AAP","APN","ANN","AAN","HHH","HOO","HHO","HOS","HOC","HOP","HON","HSS","HHS","HSC","HSP","HSN","HCC","HHC","HCP","HCN","HPP","HHP","HPN","HNN","HHN","OOO","OSS","OOS","OSC","OSP","OSN","OCC","OOC","OCP","OCN","OPP","OOP","OPN","ONN","OON","SSS","SCC","SSC","SCP","SCN","SPP","SSP","SPN","SNN","SSN","CCC","CPP","CCP","CNN","CCN","PPP","PNN","PPN","NNN"};
int TRIPLETS[TOT_TRIP] = {111,100,110,103,105,107,102,106,104,133,113,135,137,132,136,134,155,115,157,152,156,154,177,117,172,176,174,122,112,126,124,166,116164,144,114,888,833,883,835,832,836,834,855,885,857,852,856,854,877,887,872,876,874,822,882,826,824,866,886,864,844,884,333,355,335,357,352,356,354,377,337,372,376,374,322,332,326,324,366,336,364,344,334,555,577,557,572,576,574,522,552,526,524,566,556,564,544,554,777,722,772,726,724,766,776,764,744,774,222,266,226,244,224,666,644,664,444};

// Create 
//void Conformer::BoltEnCalc(const char * naccess_file)
//REAL Conformer::BoltEnCalc(CChain & m_chain, int m_length)
REAL Conformer::BoltEnCalc(const char *  nFile)
{
  // Gemotric Center of conformer
  //cout << m_length << endl;
  // Assgin residue type to residues
  int res_type[m_length];

  for (int i = 0; i < m_length; i++ )
  {
    //cout << i << "  " << m_chain.getLink(2*i)->getType() << "  " << m_chain.getLink(2*i+1)->getType() << endl; 
    //if (i == BBN || i == NTR || i == CTR)
    //if (m_chain.getLink(2*i)->getType() == BBN || m_chain.getLink(2*i)->getType() == NTR || m_chain.getLink(2*i)->getType() == CTR || m_chain.getLink(2*i)->getType() == BBP)
      //continue;
    //else
      //cout << i << " " << m_chain.getLink(2*i+1)->getType() << endl;
      if (m_chain.getLink(2*i+1)->getType() == ASP || m_chain.getLink(2*i+1)->getType() == GLU)
      {
//        cout << "Here" << endl;
        res_type[i] = a;
      }
      else if (m_chain.getLink(2*i+1)->getType() == ARG || m_chain.getLink(2*i+1)->getType() == HIS || m_chain.getLink(2*i+1)->getType() == LYS)
        res_type[i] = b;
      else if (m_chain.getLink(2*i+1)->getType() == CYS)
        res_type[i] = c;
      else if (m_chain.getLink(2*i+1)->getType() == TRP || m_chain.getLink(2*i+1)->getType() == TYR || m_chain.getLink(2*i+1)->getType() == PHE || m_chain.getLink(2*i+1)->getType() == MET || m_chain.getLink(2*i+1)->getType() == LEU || m_chain.getLink(2*i+1)->getType() == ILE || m_chain.getLink(2*i+1)->getType() == VAL)
        res_type[i] = h;
      else if (m_chain.getLink(2*i+1)->getType() == ASN || m_chain.getLink(2*i+1)->getType() == GLN)
        res_type[i] = n;
      else if (m_chain.getLink(2*i+1)->getType() == SER || m_chain.getLink(2*i+1)->getType() == THR)
        res_type[i] = o;
      else if (m_chain.getLink(2*i+1)->getType() == PRO)
        res_type[i] = p;
      else if (m_chain.getLink(2*i+1)->getType() == ALA || m_chain.getLink(2*i+1)->getType() == GLY)
      {
        //cout << "Here" << endl;
        res_type[i] = s;
      }
      else
      {
           cout << "Error in amino acid type" << endl;
           exit(0);
      }
  } 

  // Load output of naccess
  //char naccess_file[100] = "1LFH.rsa";
  
  vector<double> res_asa;
//  load_asa(const char * naccess_file, res_asa);
  load_asa(nFile, res_asa);

  tetra_center = new double*[m_length-3];
  for (int i = 0; i < m_length-3; i++)
     tetra_center[i] = new double[3];

  tetra_vert = new double**[m_length-3];
     for (i = 0; i < m_length-3; i++)
     {
       tetra_vert[i] = new double*[4];
       for (int j = 0; j < 4; j++)
         tetra_vert[i][j] = new double[3];
     }

  center_idx = new int*[m_length-3];
  for (i = 0; i < m_length-3; i++)
    center_idx[i] = new int[m_length-4];

  triplet = new int**[m_length-3];
  for (i = 0; i < m_length-3; i++)
  {
    triplet[i] = new int*[4];
    for (int j = 0; j < 4; j++)
      //triplet[i][j] = new int[m_length-4];
      triplet[i][j] = new int[m_length];
  }

  triplet_rsa = new int*[m_length-3];
  for (i = 0; i < m_length-3; i++)
     //triplet_rsa[i] = new int[m_length-4];
     triplet_rsa[i] = new int[4];

  for (i = 0; i < m_length-3; i++)
  {
    for (int j = 0; j < 4; j++)
     {
       for (int l = 0; l < m_length; l++)
          triplet[i][j][l] = 0;
     }
  }

  for (i = 0; i < m_length-3; i++)
    for (int j = 0; j < m_length-4; j++)
       center_idx[i][j] = -1;
  // Make tetrahedrons and assign points to them
  assignTetrahedronToPoints(tetra_center, tetra_vert, center_idx, triplet);
  cout << "Assigned tetrahedron" << endl;
  //cout << "2. Here" << endl;
  assignSurfaceCategory(res_asa, triplet_rsa, triplet);
  cout << "Surface assigned" << endl;
//  cout << "2. Here" << endl;
  tripB = new int[TOT_TRIP];
  quadB = new int*[TOT_TRIP];
  for (i = 0; i < TOT_TRIP; i++)
    quadB[i] = new int[TOT_SING];
  tripE = new int[TOT_TRIP];
  quadE = new int*[TOT_TRIP];
  for (i = 0; i < TOT_TRIP; i++)
    quadE[i] = new int[TOT_SING];
  tripI = new int[TOT_TRIP];
  quadI = new int*[TOT_TRIP];
  for (i = 0; i < TOT_TRIP; i++)
    quadI[i] = new int[TOT_SING];
  sing = new int[TOT_SING];

  calc_specific_quadruplets(triplet_rsa, tetra_center, res_type, triplet, 1, tripB, quadB);
  cout << "Buried quadruplet probability calculated" << endl;
  calc_specific_quadruplets(triplet_rsa, tetra_center, res_type, triplet, 2, tripE, quadE);
  calc_specific_quadruplets(triplet_rsa, tetra_center, res_type, triplet, 3, tripI, quadI);

  calc_specific_singlets(sing);

  //Scoring files read somewhere else

  int i,j;
  double score;
  for (i = 0; i < TOT_TRIP; i++)
  {
    for (j = 0; j < TOT_SING; j++)
    {
      if (m_chain.quadB_prob[(i * TOT_SING)+j] * quadB[i][j] == 0.0)
        continue;
      else
      {
        if (m_chain.tripB_prob[i]*tripB[i]*m_chain.sing_prob[j]*sing[j] == 0.0)
           score += log(m_chain.quadB_prob[(i*TOT_SING)+j]*quadB[i][j]/0.00001);
        else
           score += log(m_chain.quadB_prob[i*TOT_SING+j]*quadB[i][j]/m_chain.tripB_prob[i]*tripB[i]*m_chain.sing_prob[j]*sing[j]);
      }
    }
  }
  for (i = 0; i < TOT_TRIP; i++)
  {
    for (j = 0; j < TOT_SING; j++)
    {
      if (m_chain.quadE_prob[i*TOT_SING+j]*quadE[i][j] == 0.0)
        continue;
      else
      {
        if (m_chain.tripE_prob[i]*tripE[i]*m_chain.sing_prob[j]*sing[j] == 0.0)
           score += log(m_chain.quadE_prob[i*TOT_SING+j]*quadE[i][j]/0.00001);
        else
           score += log(m_chain.quadE_prob[i*TOT_SING+j]*quadE[i][j]/m_chain.tripE_prob[i]*tripE[i]*m_chain.sing_prob[j]*sing[j]);
      }
    }
  }
  for (i = 0; i < TOT_TRIP; i++)
  {
    for (j = 0; j < TOT_SING; j++)
    {
      if (m_chain.quadI_prob[i*TOT_SING+j]*quadI[i][j] == 0.0)
        continue;
      else
      {
        if (m_chain.tripI_prob[i]*tripI[i]*m_chain.sing_prob[j]*sing[j] == 0.0)
           score += log(m_chain.quadI_prob[i*TOT_SING+j]*quadI[i][j]/0.00001);
        else
           score += log(m_chain.quadI_prob[i*TOT_SING+j]*quadI[i][j]/m_chain.tripI_prob[i]*tripI[i]*m_chain.sing_prob[j]*sing[j]);
      }
    }
  }
  return score;

}

void Conformer::load_asa(const char * nFile, vector<double> & res_asa)
{
  //cout << nFile << endl;
  ifstream fin(nFile);
  if (!fin.is_open())
    {
      cout << "Could not open input file " << fin << endl;
      exit(0);
    }

  //vector<double> res_asa;
  //string buf;
  char buf[100];

  float abs, rel, sidech;
  int nn;
  char aa[9], ad[5], buf2[6];
  float buf3, buf4, buf5, buf6, buf7, buf8, buf9, buf10, buf11, buf12;
  fin.getline(buf, 150);
  fin.getline(buf, 150);
  fin.getline(buf, 150);
  fin.getline(buf, 150);
  while (!fin.eof())
    {
      /*
      fin >> buf2;
      cout << fin << endl;
      if (strcmp(buf2,"END") == 0)
         break;
      fin >> aa;
      fin >> ad;
      fin >> nn;
      fin >> abs;
      fin >> rel;
      fin >> sidech;
      fin >> buf3;
      fin >> buf4;
      fin >> buf5;
      fin >> buf6;
      fin >> buf7;
      fin >> buf8;
      fin >> buf9;
      fin >> buf10;
      fin >> buf11;
      fin >> buf12;
      */
      fin.getline(buf, 150);
      //exit(0);
      if (strlen(buf) == 0)
        continue;
      int res = sscanf(buf, "%s %s %s %d %f %f %f %f %f %f %f %f %f %f %f %f %f", buf2, aa, ad, &nn, &abs, &rel, &sidech, &buf3, &buf4, &buf5, &buf6, &buf7, &buf8, &buf9, &buf10, &buf11, &buf12);
      //cout << sidech << endl;
      if (strcmp(buf2,"END") == 0 || strcmp(buf2,"CHAIN") == 0 || strcmp(buf2,"TOTAL") == 0 )
      {
//        cout << " Here " << endl;
        continue;
      }
//      else
//        cout << sidech << endl;
      
      res_asa.push_back(sidech);
      
     }
  fin.close();
}

void Conformer::assignTetrahedronToPoints(double **tetra_center, double ***tetra_vert, int **center_idx, int ***triplet)
{
   int i,j,k,l;
   double x_vec, y_vec, z_vec;
   double dx, dy, dz;
   double pt_center[3];
   double d0, d1, d2, d3, d4;
   double ds1, ds2, ds3;
   float ds, ss;
   l = 0;
   d0 = 0.0;
   d1 = 0.0;
   d2 = 0.0;
   d3 = 0.0;
   d4 = 0.0;
   int sj;
   for (i = 0; i < m_length-3; i++)
   {
      //cout << m_chain.getLink(2*i+1)->getCenter()[0] << " " << m_chain.getLink(2*i+3)->getCenter()[0] << " " << m_chain.getLink(2*i+5)->getCenter()[0]+m_chain.getLink(2*i+7)->getCenter()[0] << endl;
      //exit(0);
      tetra_center[i][0] = (m_chain.getLink(2*i+1)->getCenter()[0]+m_chain.getLink(2*i+3)->getCenter()[0]+m_chain.getLink(2*i+5)->getCenter()[0]+m_chain.getLink(2*i+7)->getCenter()[0])/4.0;
      tetra_center[i][1] = (m_chain.getLink(2*i+1)->getCenter()[1]+m_chain.getLink(2*i+3)->getCenter()[1]+m_chain.getLink(2*i+5)->getCenter()[1]+m_chain.getLink(2*i+7)->getCenter()[1])/4.0;
      tetra_center[i][2] = (m_chain.getLink(2*i+1)->getCenter()[2]+m_chain.getLink(2*i+3)->getCenter()[2]+m_chain.getLink(2*i+5)->getCenter()[2]+m_chain.getLink(2*i+7)->getCenter()[2])/4.0;
      
   }
   for (i = 0; i < m_length-3; i++)
   {
      for (j = 0; j < 4; j++)
      {
         x_vec = m_chain.getLink((2*i+1)+2*j)->getCenter()[0] - tetra_center[i][0];
         y_vec = m_chain.getLink((2*i+1)+2*j)->getCenter()[1] - tetra_center[i][1];
         z_vec = m_chain.getLink((2*i+1)+2*j)->getCenter()[2] - tetra_center[i][2];
         //printf("%.5f   %.5f   %.5f\\n",x_vec,y_vec,z_vec);
         //printf("%.5f   \\n",sqrt(x_vec*x_vec+y_vec*y_vec+z_vec*z_vec));
         x_vec = x_vec/sqrt(x_vec*x_vec+y_vec*y_vec+z_vec*z_vec);
         y_vec = y_vec/sqrt(x_vec*x_vec+y_vec*y_vec+z_vec*z_vec);
         z_vec = z_vec/sqrt(x_vec*x_vec+y_vec*y_vec+z_vec*z_vec);
         //printf("%.5f   %.5f   %.5f\\n",x_vec,y_vec,z_vec);   

         tetra_vert[i][j][0] =  tetra_center[i][0] + 8.0*x_vec;
         tetra_vert[i][j][1] =  tetra_center[i][1] + 8.0*y_vec;
         tetra_vert[i][j][2] =  tetra_center[i][2] + 8.0*z_vec;
         //printf("%.3f   %.3f   %.3f\\n",tetra_center(i,0),tetra_center(i,1),tetra_center(i,2));
      }
      //exit(0);
      for (k = 0; k < m_length; k++)
      {
         if (k >= i && k < i+4)
            continue;
         dx = m_chain.getLink(2*k+1)->getCenter()[0] - tetra_center[i][0];
         dy = m_chain.getLink(2*k+1)->getCenter()[1] - tetra_center[i][1];
         dz = m_chain.getLink(2*k+1)->getCenter()[2] - tetra_center[i][2];
         if (sqrt(dx*dx+dy*dy+dz*dz) <= 8.0)
         {
            //cout << "--------------------------" << i << " " << l << " " << k << " " << sqrt(dx*dx+dy*dy+dz*dz) << endl;
            //exit(0);
            center_idx[i][l] = k;
            //cout << center_idx[i][l] << endl;
            l += 1;
         }
         else
         {
            //cout << i << " " << l << " " << k << " " << sqrt(dx*dx+dy*dy+dz*dz) << endl;
            //cout << center_idx[i][l] << endl;
            continue;
         }
      }
      l = 0;
   }
   /*
   for (i = 0; i < m_length-3 ; i++)
      for (j = 0; j < 4; j++)
         for (k = 0; k < m_length-4; k++)
           cout << triplet[i][j][k] << endl;
   exit(0);
   */
   
   //cout << "********************" << endl; 
   /*
   for (i = 0; i < m_length-3; i++)
   {
      for (l = 0; l < m_length-4; l++)
         //cout << i << " " << l << " " << center_idx[i][l] << endl;
         cout <<  l << endl;

      exit(0);   
   }
   */
   //cout << "1. Here" << endl;
   //exit(0);
   for (i = 0; i < m_length-3; i++)
   {
     for (l = 0; l < m_length-4; l++)
     {

      ss = 999.0;
      sj = 0;
      //cout << center_idx[i][l] << endl;
      if (center_idx[i][l] == -1)
         continue;// All the l's are filled first and rest is left with -1. No point of 'for loop'
      //cout << "2. Here " << i << " " << l << " " << center_idx[i][l] << endl;
      k = center_idx[i][l];
      //cout << "2. Here " << i << " " << l << " " << k << endl;
      pt_center[0] = m_chain.getLink(2*k+1)->getCenter()[0];
      pt_center[1] = m_chain.getLink(2*k+1)->getCenter()[1];
      pt_center[2] = m_chain.getLink(2*k+1)->getCenter()[2];
      //printf("%d ////// %d\\n",i,k);
      //printf("%.3f   %.3f   %.3f    \\n",pt_center[0],pt_center[1],pt_center[2]);
      // Ref. http://mathworld.wolfram.com/Line-PlaneIntersection.html
      // Ref. http://steve.hollasch.net/cgindex/geometry/ptintet.html
      //Tetrahedron 1 :- center, i, i+1, i+2.
      //Tetrahedron 2 :- center, i+1, i+2, i+3.
      //Tetrahedron 3 :- center, i+2, i+3, i.
      //Tetrahedron 3 :- center, i+3, i, i+1.
      for (j = 0; j < 4; j++)
      {
      dx = m_chain.getLink((2*i+1)+2*(j%4))->getCenter()[0] - pt_center[0];
      dy = m_chain.getLink((2*i+1)+2*(j%4))->getCenter()[1] - pt_center[1];
      dz = m_chain.getLink((2*i+1)+2*(j%4))->getCenter()[2] - pt_center[2];
      ds1 = sqrt(dx*dx+dy*dy+dz*dz);
      dx = m_chain.getLink((2*i+1)+2*((j+1)%4))->getCenter()[0] - pt_center[0];
      dy = m_chain.getLink((2*i+1)+2*((j+1)%4))->getCenter()[1] - pt_center[1];
      dz = m_chain.getLink((2*i+1)+2*((j+1)%4))->getCenter()[2] - pt_center[2];
      ds2 = sqrt(dx*dx+dy*dy+dz*dz);
      dx = m_chain.getLink((2*i+1)+2*((j+2)%4))->getCenter()[0] - pt_center[0];
      dy = m_chain.getLink((2*i+1)+2*((j+2)%4))->getCenter()[1] - pt_center[1];
      dz = m_chain.getLink((2*i+1)+2*((j+2)%4))->getCenter()[2] - pt_center[2];
      ds3 = sqrt(dx*dx+dy*dy+dz*dz);
      ds = ds1+ds2+ds3;
      //printf("%.3f  \\n",ds);
      //cout << "3. Here " << ds << endl;
      d0 = tetra_vert[i][j%4][0]*(tetra_vert[i][(j+1)%4][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][(j+1)%4][2]*tetra_vert[i][(j+2)%4][1])-tetra_vert[i][(j+1)%4][0]*(tetra_vert[i][j%4][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][j%4][2]*tetra_vert[i][(j+2)%4][1])+tetra_vert[i][(j+2)%4][0]*(tetra_vert[i][j%4][1]*tetra_vert[i][(j+1)%4][2]-tetra_vert[i][j%4][2]*tetra_vert[i][(j+1)%4][1]);
      d0 = d0 - (tetra_center[i][0]*(tetra_vert[i][(j+1)%4][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][(j+1)%4][2]*tetra_vert[i][(j+2)%4][1])-tetra_vert[i][(j+1)%4][0]*(tetra_center[i][1]*tetra_vert[i][(j+2)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+2)%4][1])+tetra_vert[i][(j+2)%4][0]*(tetra_center[i][1]*tetra_vert[i][(j+1)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+1)%4][1]));
      d0 = d0 + tetra_center[i][0]*(tetra_vert[i][j%4][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][j%4][2]*tetra_vert[i][(j+2)%4][1])-tetra_vert[i][j%4][0]*(tetra_center[i][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][(j+2)%4][1]*tetra_center[i][2])+tetra_vert[i][(j+2)%4][0]*(tetra_center[i][1]*tetra_vert[i][j%4][2]-tetra_center[i][2]*tetra_vert[i][j%4][1]);
      d0 = d0 - (tetra_center[i][0]*(tetra_vert[i][j%4][1]*tetra_vert[i][(j+1)%4][2]-tetra_vert[i][j%4][2]*tetra_vert[i][(j+1)%4][1])-tetra_vert[i][j%4][0]*(tetra_center[i][1]*tetra_vert[i][(j+1)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+1)%4][1])+tetra_vert[i][(j+1)%4][0]*(tetra_center[i][1]*tetra_vert[i][j%4][2]-tetra_center[i][2]*tetra_vert[i][j%4][1]));
      //printf("%.3f   \\n",d0);
//      cout << "d0 " << d0 << endl;
   ///////////////////////////////////////////////////////////
      d1 = tetra_vert[i][j%4][0]*(tetra_vert[i][(j+1)%4][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][(j+1)%4][2]*tetra_vert[i][(j+2)%4][1])-tetra_vert[i][(j+1)%4][0]*(tetra_vert[i][j%4][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][j%4][2]*tetra_vert[i][(j+2)%4][1])+tetra_vert[i][(j+2)%4][0]*(tetra_vert[i][j%4][1]*tetra_vert[i][(j+1)%4][2]-tetra_vert[i][j%4][2]*tetra_vert[i][(j+1)%4][1]);
      d1 = d1 - (pt_center[0]*(tetra_vert[i][(j+1)%4][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][(j+1)%4][2]*tetra_vert[i][(j+2)%4][1])-tetra_vert[i][(j+1)%4][0]*(pt_center[1]*tetra_vert[i][(j+2)%4][2]-pt_center[2]*tetra_vert[i][(j+2)%4][1])+tetra_vert[i][(j+2)%4][0]*(pt_center[1]*tetra_vert[i][(j+1)%4][2]-pt_center[2]*tetra_vert[i][(j+1)%4][1]));
      d1 = d1 + pt_center[0]*(tetra_vert[i][(j+0)%4][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][(j+0)%4][2]*tetra_vert[i][(j+2)%4][1])-tetra_vert[i][(j+0)%4][0]*(pt_center[1]*tetra_vert[i][(j+2)%4][2]-pt_center[2]*tetra_vert[i][(j+2)%4][1])+tetra_vert[i][(j+2)%4][0]*(pt_center[1]*tetra_vert[i][(j+0)%4][2]-pt_center[2]*tetra_vert[i][(j+0)%4][1]);
      d1 = d1 - (pt_center[0]*(tetra_vert[i][(j+0)%4][1]*tetra_vert[i][(j+1)%4][2]-tetra_vert[i][(j+0)%4][2]*tetra_vert[i][(j+1)%4][1])-tetra_vert[i][(j+0)%4][0]*(pt_center[1]*tetra_vert[i][(j+1)%4][2]-pt_center[2]*tetra_vert[i][(j+1)%4][1])+tetra_vert[i][(j+1)%4][0]*(pt_center[1]*tetra_vert[i][(j+0)%4][2]-pt_center[2]*tetra_vert[i][(j+0)%4][1]));
      //printf("%.3f   \\n",d1);
//      cout << "d1 " << d1 << endl;
   ///////////////////////////////////////////////////////////
      d2 = pt_center[0]*(tetra_vert[i][(j+1)%4][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][(j+1)%4][2]*tetra_vert[i][(j+2)%4][1])-tetra_vert[i][(j+1)%4][0]*(pt_center[1]*tetra_vert[i][(j+2)%4][2]-pt_center[2]*tetra_vert[i][(j+2)%4][1])+tetra_vert[i][(j+2)%4][0]*(pt_center[1]*tetra_vert[i][(j+1)%4][2]-pt_center[2]*tetra_vert[i][(j+1)%4][1]);
      d2 = d2 - (tetra_center[i][0]*(tetra_vert[i][(j+1)%4][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][(j+1)%4][2]*tetra_vert[i][(j+2)%4][1])-tetra_vert[i][(j+1)%4][0]*(tetra_center[i][1]*tetra_vert[i][(j+2)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+2)%4][1])+tetra_vert[i][(j+2)%4][0]*(tetra_center[i][1]*tetra_vert[i][(j+1)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+1)%4][1]));
      d2 = d2 + tetra_center[i][0]*(pt_center[1]*tetra_vert[i][(j+2)%4][2]-pt_center[2]*tetra_vert[i][(j+2)%4][1])-pt_center[0]*(tetra_center[i][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][(j+2)%4][1]*tetra_center[i][2])+tetra_vert[i][(j+2)%4][0]*(tetra_center[i][1]*pt_center[2]-tetra_center[i][2]*pt_center[1]);
      d2 = d2 - (tetra_center[i][0]*(pt_center[1]*tetra_vert[i][(j+1)%4][2]- pt_center[2]*tetra_vert[i][(j+1)%4][1])- pt_center[0]*(tetra_center[i][1]*tetra_vert[i][(j+1)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+1)%4][1])+tetra_vert[i][(j+1)%4][0]*(tetra_center[i][1]* pt_center[2]-tetra_center[i][2]* pt_center[1]));
      //printf("%.3f   \\n",d2);
//      cout << "d2 " << d2 << endl;
   ///////////////////////////////////////////////////////////
      d3 = tetra_vert[i][(j+0)%4][0]*(pt_center[1]*tetra_vert[i][(j+2)%4][2]-pt_center[2]*tetra_vert[i][(j+2)%4][1])-pt_center[0]*(tetra_vert[i][(j+0)%4][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][(j+2)%4][1]*tetra_vert[i][(j+0)%4][2])+tetra_vert[i][(j+2)%4][0]*(tetra_vert[i][(j+0)%4][1]*pt_center[2]-tetra_vert[i][(j+0)%4][2]*pt_center[1]);
      d3 = d3 - (tetra_center[i][0]*(pt_center[1]*tetra_vert[i][(j+2)%4][2]-pt_center[2]*tetra_vert[i][(j+2)%4][1])-pt_center[0]*(tetra_center[i][1]*tetra_vert[i][(j+2)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+2)%4][1])+tetra_vert[i][(j+2)%4][0]*(tetra_center[i][1]*pt_center[2]-tetra_center[i][2]*pt_center[1]));
      d3 = d3 + tetra_center[i][0]*(tetra_vert[i][(j+0)%4][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][(j+0)%4][2]*tetra_vert[i][(j+2)%4][1])-tetra_vert[i][(j+0)%4][0]*(tetra_center[i][1]*tetra_vert[i][(j+2)%4][2]-tetra_vert[i][(j+2)%4][1]*tetra_center[i][2])+tetra_vert[i][(j+2)%4][0]*(tetra_center[i][1]*tetra_vert[i][(j+0)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+0)%4][1]);
      d3 = d3 - (tetra_center[i][0]*(tetra_vert[i][(j+0)%4][1]*pt_center[2]-tetra_vert[i][(j+0)%4][2]*pt_center[1])-tetra_vert[i][(j+0)%4][0]*(tetra_center[i][1]*pt_center[2]-tetra_center[i][2]*pt_center[1])+pt_center[0]*(tetra_center[i][1]*tetra_vert[i][(j+0)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+0)%4][1]));
      //printf("%.3f   \\n",d3);
//      cout << "d3 " << d3 << endl;
   ///////////////////////////////////////////////////////////
      d4 = tetra_vert[i][(j+0)%4][0]*(tetra_vert[i][(j+1)%4][1]*pt_center[2]-tetra_vert[i][(j+1)%4][2]*pt_center[1])-tetra_vert[i][(j+1)%4][0]*(tetra_vert[i][(j+0)%4][1]*pt_center[2]-tetra_vert[i][(j+0)%4][2]*pt_center[1])+pt_center[0]*(tetra_vert[i][(j+0)%4][1]*tetra_vert[i][(j+1)%4][2]-tetra_vert[i][(j+0)%4][2]*tetra_vert[i][(j+1)%4][1]);
      d4 = d4 - (tetra_center[i][0]*(tetra_vert[i][(j+1)%4][1]*pt_center[2]-tetra_vert[i][(j+1)%4][2]*pt_center[1])-tetra_vert[i][(j+1)%4][0]*(tetra_center[i][1]*pt_center[2]-tetra_center[i][2]*pt_center[1])+pt_center[0]*(tetra_center[i][1]*tetra_vert[i][(j+1)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+1)%4][1]));
      d4 = d4 + tetra_center[i][0]*(tetra_vert[i][(j+0)%4][1]*pt_center[2]-tetra_vert[i][(j+0)%4][2]*pt_center[1])-tetra_vert[i][(j+0)%4][0]*(tetra_center[i][1]*pt_center[2]-pt_center[1]*tetra_center[i][2])+pt_center[0]*(tetra_center[i][1]*tetra_vert[i][(j+0)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+0)%4][1]);
      d4 = d4 - (tetra_center[i][0]*(tetra_vert[i][(j+0)%4][1]*tetra_vert[i][(j+1)%4][2]-tetra_vert[i][(j+0)%4][2]*tetra_vert[i][(j+1)%4][1])-tetra_vert[i][(j+0)%4][0]*(tetra_center[i][1]*tetra_vert[i][(j+1)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+1)%4][1])+tetra_vert[i][(j+1)%4][0]*(tetra_center[i][1]*tetra_vert[i][(j+0)%4][2]-tetra_center[i][2]*tetra_vert[i][(j+0)%4][1]));
      //printf("%.3f   \\n",d4);
      //cout << "d0 " << d0 << " d1 " << d1 << " d2 " << d2 << " d3 " << d3 << " d4 " << d4 << endl;      
   ///////////////////////////////////////////////////////////
      //cout << "4. Here" << ds << " " << ss << endl;
      if ((d0 <= 0.0 && d1 <= 0.0 && d2 <= 0.0 && d3 <= 0.0 && d4 <= 0.0) || (d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0 && d3 >= 0.0 && d4 >= 0.0))
      {
         //cout << ds << " " << ss << endl;
         //cout << i << " " << sj << " " << j << " " << k << endl;
         //cout << "triplet : " << triplet[i][j][k] << endl;
         //dx = tetra_vert[i][(j+1)%4][0] - pt_center[0];
         //dy = tetra_vert[i][(j+1)%4][1] - pt_center[1];
         //dz = tetra_vert[i][(j+1)%4][2] - pt_center[2];
         //ds1 = sqrt(dx*dx+dy*dy+dz*dz);
         
         //point is inside center, i, i+1, i+2
         //printf("%.3f   %.3f\\n",ds,ss);
         if (ds < ss)
         {
            //cout << "ss " << ss << endl;
            //cout << "ds " << ds << endl;
            //cout << "sj " << sj << endl;
            //cout << "i  " << i << endl;
            //cout << "j  " << j << endl;
            //cout << "k  " << k << endl;
            //cout << "A:  " << triplet[i][sj][k] << endl;
            triplet[i][sj][k] = 0;
            //cout << "B:  " << triplet[i][sj][k] << endl;
            ss = ds;
            sj = j;
            triplet[i][j][k] = 1;
            //cout << "C:  " << triplet[i][j][k] << endl;
         }
         //printf("%d   %d  %d  %d  \\n",i,j,k,triplet(i,j,k));
         d0 = 0.0;
         d1 = 0.0;
         d2 = 0.0;
         d3 = 0.0;
         d4 = 0.0;
      }
      else
      {
         //point is outside
         d0 = 0.0;
         d1 = 0.0;
         d2 = 0.0;
         d3 = 0.0;
         d4 = 0.0;
         //continue;
      }
//      cout << "4. Here " << endl;
//      if (i == 5 || i == 6)
//        cout << "J completed: " << j << endl;
      }
      ss = 999.99;
      sj = 0;
     }
   }
   //////////////////////////////////////////////////////////
}

void Conformer::assignSurfaceCategory(vector<double> & res_asa, int  **triplet_rsa, int ***triplet)
{

  //cout << "3. Here " << endl;
  for (int i = 0; i < m_length-3; i++)
  {
     for (int j = 0; j < 4; j++)
     {
        for (int k = 0; k < m_length-4; k++)
        //for (int k = 0; k < m_length; k++)
        {
           if (triplet[i][j][k] == 1)
           {
              if (res_asa[i+j] < 20.0 and res_asa[i+(j+1)%4] < 20.0 and res_asa[i+(j+2)%4] < 20.0)
                  triplet_rsa[i][j] = 1;
              else if (res_asa[i+j] >= 20.0 and res_asa[i+(j+1)%4] >= 20.0 and res_asa[i+(j+2)%4] >= 20.0)
                  triplet_rsa[i][j] = 2;
              else
                  triplet_rsa[i][j] = 3;
           }
        }
      }
      //cout << i << endl;
   }
}

void Conformer::calc_specific_quadruplets(int **triplet_rsa, double **tetra_center, int *res_type, int ***triplet, int triplet_prop, int *trip, int **quad)
{

  //char res1[2], res2[2], res3[2], res4[2], res5[2], res6[2];
  int res1, res2, res3, res4, res5, res6;
  int triplet_type, triplet_type1, triplet_type2, triplet_type3, triplet_type4, triplet_type5, triplet_type6;
  int res_class;
  int si;

  for ( int i = 0; i < m_length-3; i++)
  {
      cout << i << endl;
      //for (int j = 0; j < m_length; j++)
      for (int j = 0; j < 4; j++)
      {
         //cout << "j " << j << endl;
         if (triplet_rsa[i][j] == triplet_prop)
         {
            //cout << i << " " << j << " " << triplet_rsa[i][j] << endl;
            res1 = res_type[i+j];
            res2 = res_type[i+(j+1)%4];
            res3 = res_type[i+(j+2)%4];
            //cout << res1 << " " << res2 << " " << res3 << endl;

            if (res1 == res2 and res2 == res3 and res1 == res3)
            {
              triplet_type = res1*100+res2*10+res3;
              for (int k = 0; k < TOT_TRIP; k++)
                //if (triplet_type == *TRIPLETS[k])
                if (triplet_type == TRIPLETS[k])
                {
                 res_class = TRIPLET_IDX[k];
                 trip[res_class] += 1;
                }
              continue;
            }
            else if (res1 == res2 or res2 == res3 or res1 == res3)
            {
              if (res1 != res2 and res2 == res3)
                 {
                   triplet_type1 = res1*100+res2*10+res3;
                   triplet_type2 = res2*100+res3*10+res1;
                   triplet_type3 = res2*100+res1*10+res3;
                 }
              else if (res1 != res2 and res1 == res3)
                 {
                   triplet_type1 = res1*100+res2*10+res3;
                   triplet_type2 = res1*100+res3*10+res2;
                   triplet_type3 = res2*100+res1*10+res3;
                   //cout << triplet_type1 << " " << triplet_type2 << " " << triplet_type3 << endl;
                 }
              else
                 {
                   triplet_type1 = res1*100+res2*10+res3;
                   triplet_type2 = res1*100+res3*10+res2;
                   triplet_type3 = res3*100+res1*10+res2;
                 }
                   for (int k = 0; k < TOT_TRIP; k++)
                     if (triplet_type1 == TRIPLETS[k])
                     {
                       res_class = TRIPLET_IDX[k];
                       trip[res_class] += 1;
                     }
                   for (int k = 0; k < TOT_TRIP; k++)
                     if (triplet_type2 == TRIPLETS[k])
                     {
                       res_class = TRIPLET_IDX[k];
                       trip[res_class] += 1;
                     }
                   for (int k = 0; k < TOT_TRIP; k++)
                     if (triplet_type3 == TRIPLETS[k])
                     {
                       res_class = TRIPLET_IDX[k];
                       trip[res_class] += 1;
                     }
                  
 
              }
            else
            {
                  triplet_type1 = res1*100+res2*10+res3;
                  triplet_type2 = res1*100+res3*10+res2;
                  triplet_type3 = res3*100+res1*10+res2;
                  triplet_type4 = res2*100+res1*10+res3;
                  triplet_type5 = res2*100+res3*10+res1;
                  triplet_type6 = res3*100+res2*10+res1;
                  //cout << triplet_type1 << " " << triplet_type2 << " " << triplet_type3 << " " << triplet_type4 << " " << triplet_type5 << " " << triplet_type6 << endl;
                  for (int k = 0; k < TOT_TRIP; k++)
                     if (triplet_type1 == TRIPLETS[k])
                     {
                       res_class = TRIPLET_IDX[k];
                       trip[res_class] += 1;
                     }
                  for (k = 0; k < TOT_TRIP; k++)
                     if (triplet_type2 == TRIPLETS[k])
                     {
                       res_class = TRIPLET_IDX[k];
                       trip[res_class] += 1;
                     }
                   for (k = 0; k < TOT_TRIP; k++)
                     if (triplet_type3 == TRIPLETS[k])
                     {
                       res_class = TRIPLET_IDX[k];
                       trip[res_class] += 1;
                     }
                  for (k = 0; k < TOT_TRIP; k++)
                     if (triplet_type4 == TRIPLETS[k])
                     {
                       res_class = TRIPLET_IDX[k];
                       trip[res_class] += 1;
                     }
                  for (k = 0; k < TOT_TRIP; k++)
                     if (triplet_type5 == TRIPLETS[k])
                     {
                       res_class = TRIPLET_IDX[k];
                       trip[res_class] += 1;
                     }
                  for (k = 0; k < TOT_TRIP; k++)
                     if (triplet_type6 == TRIPLETS[k])
                     {
                       res_class = TRIPLET_IDX[k];
                       trip[res_class] += 1;
                     }
                  //cout << "res_class " << res_class << endl;

            }

              for (int l = 0; l < m_length-4; l++)
              //for (int l = 0; l < m_length; l++)
              {
                //cout << "l " << l << " " << triplet[i][j][l] << endl;
                //cout << triplet[i][j][l] << endl;
                if (triplet[i][j][l] == 1)
                {
                  //cout << i << " " << j << " " << l << " " << triplet[i][j][l] << endl;
                  //for (int k = 0; k < TOT_SING; k++)
                  //{
                    //if (m_links[]->getType() == *SINGLETS[k])
                    int mc = m_chain.getLink(2*l+1)->getType();
                    //cout << "mc " << mc << " " << res_class << endl;
                    int si = SINGLET_IDX[mc];
                  //}
                  quad[res_class][si] += 1;
                }
              }
         }
      }
  }
}

void Conformer::calc_specific_singlets(int *sing)
{

   for (int i = 0; i < m_length; i++)
   {
     int k = m_chain.getLink(2*i+1)->getType();
     sing[SINGLET_IDX[k]] += 1;
      //specific_sing[singlet_sets.index(residues[i])] += 1;
   }
}
