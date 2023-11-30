#ifndef _CGEEF1_H
#define _CGEEF1_H

#include <cmath>
#include <cassert>

// for new g++ compilers use:
#include <ext/hash_map>
// for older versions use
//#include <hash_map>

#include "MatVec.h"

#define T_0 298.15
#define MAX_ROTAMER_SIZE 13 // size of ARG
#define MAX_ROTAMER_PAIRS (MAX_ROTAMER_SIZE*(MAX_ROTAMER_SIZE-1))
extern double CUTOFF_DISTANCE;
extern double CUTOFF_DISTANCE_2; 
extern double CLASH_CUTOFF_DISTANCE;

// The Boltzman constant
#define Kb 0.00198576

//Atom types
#define NUM_ATYPES 19
#define NUM_HEAVY_TYPES 17

#define  C    0
#define  CR   1
#define  CH1E 2
#define  CH2E 3
#define  CH3E 4
#define  CR1E 5
#define  NH1  6
#define  NR   7
#define  NH2  8
#define  NH3  9
#define  NC2  10
#define  N    11
#define  OH1  12
#define  O    13
#define  OC   14
#define  S    15
#define  SH1E 16
#define  H    17
#define  HC   18

// Amino acid types
#define NUMAA 26
#define NUMLIG 1
#define NUMCONF 5

#define IND -1
#define ARG 0
#define ASN 1
#define ASP 2
#define CYS 3 
#define GLN 4
#define GLU 5
#define HIS 6
#define ILE 7 
#define LEU 8 
#define LYS 9
#define MET 10
#define PHE 11 
#define PRO 12
#define SER 13
#define THR 14
#define TRP 15
#define TYR 16
#define VAL 17
#define GLY 18
#define ALA 19
#define NTR 20 // N-terminal cap
#define CTR 21 // C-terminal cap
#define BBN 22 // C - N backbone piece
#define BBP 23 // C - N backbone piece for Proline
#define BBG 24 // C - N backbone piece for Glycine
#define BCN 25 // C - N backbone

extern char AA_NAMES[][5];
extern int NUM_ROTAMERS[];
extern int ROTAMER_START[];
//extern float ROTAMER_VALUE[][4];

extern double epsilon[NUM_ATYPES];
extern double SIGMA[NUM_ATYPES];

// The factor used to determine when we have a steric clash
#define SIGMA_FACTOR 0.5

#define TOT_TRIP 120
#define TOT_SING 20
extern const double DIELECTRIC;

int getAA(const char * aa);

inline
double compute_ES(double charge1, double charge2, double dist)
{
  // We assume the "dist" is the squared distance.
  return (DIELECTRIC * (charge1 * charge2)/(dist));
}

inline
bool isStericClash(int type1, int type2, double dist, bool is14)
{
  double sig1 = SIGMA[type1];
  double sig2 = SIGMA[type2];

  if (is14)
    {
      if (type1 < NH1)
	  sig1 = 1.9;

      if (type2 < NH1)
	  sig2 = 1.9;
    }

  double sig = SIGMA_FACTOR * (sig1 + sig2);

  return (dist < sig*sig);
}

inline
double compute_vdW(int type1, int type2, double dist, bool is14)
{
  // We assume "dist" is the distance squared.
  if (dist > CUTOFF_DISTANCE_2)
    return 0.0;

  double sig1 = SIGMA[type1];
  double eps1 = epsilon[type1];
  double sig2 = SIGMA[type2];
  double eps2 = epsilon[type2];

  if (is14)
    {
      if (type1 < NH1)
	{
	  sig1 = 1.9;
	  eps1 = -0.1;
	}

      if (type2 < NH1)
	{
	  sig2 = 1.9;
	  eps2 = -0.1;
	}
    }

  double eps = sqrt(eps1*eps2);
  double sig = (sig1 + sig2);
 
  double rat2 = (sig*sig)/dist;
  double rat6 = rat2*rat2*rat2;
		   
  double res = (eps*(rat6*rat6 - 2*rat6));
  return res;
}


// Atom volumes
extern double volume[NUM_ATYPES];
extern double deltaG_ref[NUM_ATYPES];
extern double deltaG_free[NUM_HEAVY_TYPES];
extern double deltaH_ref[NUM_HEAVY_TYPES];
extern double deltaCp_ref[NUM_HEAVY_TYPES];

inline
double compute_deltaG_ref(double T, int aType)
{
  assert(aType < NUM_HEAVY_TYPES && aType >= 0);

  double deltaS_ref = (deltaH_ref[aType] - deltaG_ref[aType]) / T_0;
  
  double res = deltaG_ref[aType] - deltaS_ref*(T - T_0) - 
    deltaCp_ref[aType]*T*log(T/T_0) + deltaCp_ref[aType]*(T-T_0);

  return res;
}

inline
double getLambda(int AAtype, int index, int aType)
{
  if (aType == NH3 || aType == NC2 || aType == OC)
    return 6.0;

  if ((AAtype == ARG && index >= 2) ||
      (AAtype == LYS && index >= 3) ||
      (AAtype == ASP) ||
      (AAtype == GLU && index >= 1) ||
      (AAtype == NTR) ||
      (AAtype = CTR))
    return 6.00;
  else 
    return 3.50;
}

extern const double SOLVATION_K;
inline
double computeSolventEffect(int aType1, int aType2, double dist, double lambda1,
			  double lambda2, bool is14)
{
  if (dist > CUTOFF_DISTANCE_2)
    return 0.0;

  double t1 = 0.0, t2 = 0.0;

  if (aType1 < NUM_HEAVY_TYPES && aType2 < NUM_HEAVY_TYPES)
  {
    double sig1 = SIGMA[aType1], sig2 = SIGMA[aType2];

  if (is14)
    {
      if (aType1 < NH1)
	sig1 = 1.9;

      if (aType2 < NH1)
	sig2 = 1.9;
    }

    // We assume "dist" is the squared distance.
    double d = sqrt(dist);

    double X12 = (d - sig1)/lambda1;
    t1 = SOLVATION_K * deltaG_free[aType1] * exp(-X12*X12) * volume[aType2] 
      / (lambda1 * dist);

    double X21 = (d - sig2)/lambda2;
    t2 = SOLVATION_K * deltaG_free[aType2] * exp(-X21*X21) * volume[aType1] 
      / (lambda2 * dist);
  }

  return -(t2 + t1);
}

// dihedral constants
#define NUM_DIHEDRALS 6

extern const double E_0[NUM_DIHEDRALS];
extern const double theta_0[NUM_DIHEDRALS];
extern const double dihed_n[NUM_DIHEDRALS];

inline
double compute_dihedral(double angle, int type)
{
  assert(type >= 0 && type < NUM_DIHEDRALS);

  return E_0[type] * (1 + cos(dihed_n[type]*angle - theta_0[type]));
}

//class CBV;
class CRss;
// A structure to hold rotamer information that can be switched quickly
// into a leaf.
struct ROTAMER
{

    //double **m_positions;
  ROTAMER(int size) : m_energy(0.0), m_bv(NULL)
  {
    m_positions = new double[size][3];
  }

  enum ROT_TYPE {UP = 1, DOWN = -1};
  double (*m_positions)[3];
  double m_sc[3];
  CRss* m_bv;
  double m_energy;
  double m_phi;
  double m_psi;
  double m_prob;
  double m_chi1mean;
  double m_chi2mean;
  double m_chi3mean;
  double m_chi4mean;
  double m_chi1std;
  double m_chi2std;
  double m_chi3std;
  double m_chi4std;
};

double compute_rotamer_energy(int type, const ROTAMER *rot, int index, float ROTAMER_VALUE[][4]);

struct LCONFORMER
{

  LCONFORMER(int size) : m_bv(NULL)
  {
    m_positions = new double[size][3];
    m_centers = new double[size/8][3];
  }

  double (*m_positions)[3];
  double (*m_centers)[3];
  CRss* m_bv;
};

struct LIGAND
{
  LIGAND(int size, int nConformers, int nGroups, int nChis) :
    m_size(size), m_nConformers(nConformers), m_nGroups(nGroups), m_nChis(nChis) {}

  static void createLigconformers(const char * dir);
  static LIGAND * m_liglist[NUMLIG];

  LCONFORMER * ll;
  std::vector<LCONFORMER *> m_lconformer;

  int m_size;
  int m_nConformers;
  int m_nGroups;
  int * m_groups;
  int m_nChis;
  int * m_aTypes;
  char ** m_aNames;
};

// A structure to hold all sidechain information
struct SIDECHAIN
{
  SIDECHAIN(int size, int nRotamers, int nChis) : 
    m_size(size), m_nRotamers(nRotamers), 
    m_nChis(nChis)
  {
    //ROTAMER *m_Urotamers;
    //m_Urotamers = new ROTAMER[nRotamers](size);
    //m_Drotamers = new ROTAMER[nRotamers](size);
    //printf("ccc %d %d\n", nRotamers, size);
    //std::vector<ROTAMER> m_Urotamers(nRotamers,size);
    //printf("ddd %p\n", m_Urotamers[0]);
    //std::vector<ROTAMER> m_Drotamers(nRotamers,size);
    /*
     m_groups = new int[nGroups];
    if (nChis > 0)
      m_chiTypes = new int[nChis];
    else 
      m_chiTypes = NULL;
    */
  }
  
  const ROTAMER & getRotamer(int index, ROTAMER::ROT_TYPE rType)
  {
    if (rType == ROTAMER::UP)
      return *m_Urotamers[index];
    else if  (rType == ROTAMER::DOWN)
      return *m_Drotamers[index];
    else
      assert(false);
  } 
  
    //printf("ddd2%p\n", m_Urotamers[0]);
  static void createSidechains(const char * dir);
  static SIDECHAIN * m_aalist[NUMAA];

    //printf("ddd3%p\n", m_Urotamers[0]);
  ROTAMER * ss;
  ROTAMER * pp;
  std::vector<ROTAMER *> m_Urotamers;
  std::vector<ROTAMER *> m_Drotamers;
  //ROTAMER * m_Drotamers;
    //printf("ddd4%p\n", m_Urotamers[0]);
  int m_nRotamers;
  int m_size;
  int m_nGroups;
  int * m_groups;
  double * m_charges;
  int * m_aTypes;
  int m_nChis;
  int * m_chiTypes;
  char ** m_aNames;

};

// Axis of rotatin for PHI and PSI angles
//extern const double NTR_PHI_AXIS[3];
//extern const double U_PHI_AXIS[3];
//extern const double D_PHI_AXIS[3];

//extern const double U_PSI_AXIS[3];
//extern const double D_PSI_AXIS[3];

// Translation between origins of the different kinds of links.
//extern const double N_C_NTR_TRANS[3];

//extern const double U_CA_C_TRANS[3];
//extern const double D_CA_C_TRANS[3];

//extern const double UD_C_CA_TRANS[3];
//extern const double DU_C_CA_TRANS[3];

// Exclusion list
inline
int to_index(int diff, int t1, int i, int t2, int j)
{
  return (diff*1024*1024 + t1*1024*32 + i*1024 + t2*32 + j);
}

// for new g++ versions use:
typedef __gnu_cxx::hash_map<int, int> EXCLUSIONS;
// for older versions use
//typedef hash_map<int, int> EXCLUSIONS;

extern EXCLUSIONS exclusion_list;
#define EXCLUDED 1
#define PAIR1_4  2
#define NOT_EXCLUDED 0

inline
int isExcluded(int diff, int type1, int id1, int type2, int id2)
{
  if (diff > 3) 
    return NOT_EXCLUDED;
  else if (diff == 0 && id1 >= id2)
    return EXCLUDED;

  EXCLUSIONS::const_iterator ex = exclusion_list.find(to_index(diff,type1,id1,type2,id2));

  if (ex == exclusion_list.end())
    return NOT_EXCLUDED;
  else if (ex->second < 3)
    return EXCLUDED;
  else if (ex->second == 3)
    return PAIR1_4;
  else
    assert(false);
}

void Initialize(const char * dir);
#endif
