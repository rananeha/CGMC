#ifndef _CGLIGCONFORMER_H
#define _CGLIGCONFORMER_H

#include <vector>
#include <cmath>

//#include "cgchain.H"
//#include "cgleaf.H"
#include "cgeef1.H" 

#define TOT_TRIP 120
#define TOT_SING 4

#define a 0
#define b 1
#define c 2
#define h 3
#define n 4
#define o 5
#define p 6
#define s 7

//extern int SING_IDX[TOT_SING];
//extern int TRIP_IDX[TOT_TRIP];
//extern char SINGLETS[TOT_SINGLETS];
//extern char * SINGLETS;
//extern char TRIPLETS[TOT_TRIPLETS];
//class CChain;

class Ligconformer
{
public:
   Ligconformer(CChain & chain, int clength) : m_chain(chain), m_length(clength), m_energy(0.0) {}
     REAL BoltEnCalc(REAL (*cxs)[3], const char * nFile);
//   REAL BoltEnCalc(CChain & m_chain, int m_length);
//   ~Conformer();

   void assignTetrahedronToPoints(double **tetra_center, double ***tetra_vert, 
				int **center_idx, int ***triplet, int, REAL (*)[3]);
   void assignSurfaceCategory(vector<double> & res_asa, int  **triplet_rsa, 
				int ***triplet, int);
   void calc_specific_quadruplets(int **triplet_rsa, double **tetra_center, 
				int *res_type, int ***triplet, int triplet_prop,
				int *trip, int **quad, int);

   void calc_specific_singlets(int *sing, int);
   
   double** getTC() const
   { return tetra_center; }
   double*** getTV() const
   { return tetra_vert; }
   int** getCI()
   { return center_idx; }
   int*** getTriplet() const
   { return triplet; }
   int** getTriplet_rsa()
   { return triplet_rsa; }
   int* getTripB()
   { return tripB; }
   int** getQuadB()
   { return quadB; }
   int* getTripE()
   { return tripE; }
   int** getQuadE()
   { return quadE; }
   int* getTripI()
   { return tripI; }
   int** getQuadI()
   { return quadI; }
   int* getSing()
   { return sing; }
   int getCEnergy()
   { return m_energy; }
   
private:

   //void BoltEnCalc();
   void load_asa(const char * nFile, vector<REAL> & res_asa);

   CChain & m_chain;
   int m_length;
   REAL m_energy;

   double **tetra_center; //array of pointer to double
   //tetra_center = new double*[AAtypes.size()-3];
   //for (int i = 0; i < AAtypes.size()-3; i++)
     //tetra_center[i] = new double[3];
   
   double ***tetra_vert;
   //tetra_vert = new double**[AAtypes.size()-3];
   //for (i = 0; i < AAtypes.size()-3; i++)
   //{
     //tetra_vert[i] = new double*[4];
     //for (int j = 0; j < 4; j++)
       //tetra_vert[i][j] = new double[3];
   //}

   int **center_idx;
   //center_idx = new int*[AAtypes.size()-3];
   //for (i = 0; i < AAtypes.size()-3; i++)
     //center_idx[i] = new int[AAtypes.size()-4];

   int ***triplet;
   //triplet = new int**[AAtypes.size()-3];
   //for (i = 0; i < AAtypes.size()-3; i++)
   //{
      //triplet[i] = new int*[4];
      //for (j = 0; j < 4; j++)
         //triplet[i][j] = new int[AAtypes.size()-4];
   //}

   int **triplet_rsa;
   //triplet_rsa = new int*[AAtypes.size()-3];
   //for (i = 0; i < AAtypes.size()-3; i++)
     //triplet_rsa = new int[AAtypes.size()-4];

   int *tripB;
   //tripB = new int[TOT_TRIP];

   int **quadB;
   //quadB = new int*[TOT_TRIP];
   //for (i = 0; i < TOT_TRIP; i++)
     //quadB[i] = new int[TOT_SING];

   int *tripE;
   //tripE = new int[TOT_TRIP];

   int **quadE;
   //quadE = new int*[TOT_TRIP];
   //for (i = 0; i < TOT_TRIP; i++)
     //quadE[i] = new int[TOT_SING];

   int *tripI;
   //tripI = new int[TOT_TRIP];

   int **quadI;
   //quadI = new int*[TOT_TRIP];
   //for (i = 0; i < TOT_TRIP; i++)
     //quadI[i] = new int[TOT_SING];

   int *sing;
   //sing = new int[TOT_SING];
};


#endif
