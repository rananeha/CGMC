
#ifndef _CGLIGAND_H
#define _CGLIGAND_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
//#include "utilities.h"
//#include "parameters.h"
//#include "ligand.h"
#include "cgleaf.H"
#include "cgchain.H"
#include "MatVec.h"

#define MC_MAX_RMSD 3.0
#define TMP_LOCAL_CM 0.1
#define TMP_LOCAL_ROT 1.0
#define SIGMAA 1.0

struct LPOSE {
    REAL (*save_pxs)[3];
    REAL (*save_cxs)[3];
};
/*
typedef struct conformation
{
    float           x[3];
    //float           charge;
    LPList          *lplist;

} conformation;
*/
/*
typedef struct atom
{
    int             atom_no; //real atom number start from 1
    int             Num; // atom number start from 0
    char            name[8];
    char            sybyl_type[10];

    struct atom     *bond[6];
	int		        bondType[6];
	int             numBonds;
	int		        bt[6];

	int             element;
	float           mass;
	int             hyb;
	int             mark;
	int             var;
	float			charge;				//charge from amber top file (timed by 18.2223, not partial charge)

    int             DonAcc;
    int             DonAccType;
    int             SolType;
    int				numHB;
    int             numH;
    int             numLP;
    int				LP_index;
    int				atmH[3];
    int				Don;
    int				Acc;

	conformation    *conf;
    conformation    *sol;
    conformation	*optimized;
	conformation    *MCsol;
	conformation	*MCtrj;

	char            FFType[8];

    int             Connect;

    int				symmatom; // counterpart of the symm atom
	int				symmgroup;

	int				ATscoring;
	int             hydroflag;
	float			VDWradius;
	
	int				scoretype;
	
	float			hydro_multi;

	int				vinatype;	//1. HBD; 2, HBA; 3, HBDA; 4, Polar; 5, Hydrophobic
	float			radius;		//atomic radius according to Xscore
} atom;
*/
/*
typedef struct sybylbond
{
    int            at1;
    int            at2;
    char           bondtype[10];
} sybylbond;
*/
/*
typedef struct ligand
{
    int             numconf;
    int		    *numsols;
    int             numatoms;
    atom            *atm;
	int				CA;		// Centroid atom, count from 0


} ligand;
*/
class MonteCarloLigand
{
public:
   MonteCarloLigand(const char * fname, REAL pointx, REAL pointy, REAL pointz) : m_pointx(pointx), m_pointy(pointy), m_pointz(pointz) 
   { load_lig_coors(fname, xcoors, ycoors, zcoors); }
   //void MC_Ligand(poses, res_list);
   //void MC_Ligand(CChain & m_chain, int, vector<int >, REAL (*)[3], REAL (*)[3]);
   void MC_Ligand(CChain & m_chain, int);
   void undoLastMove(vector<LPOSE> & accept_coords);
   //void findClash();
   //void findNeighborLeaf(CChain & m_chain, int, vector<int >);
   void findNeighborLeaf(CChain & m_chain, int);
   void WriteLigCoordinates(vector<LPOSE > & accept_coords, int num_poses, char * ofile);
   void writeLine(ofstream & fout, int index, const char * aname, const char * resname,
		const char * chainid, int resnum, REAL x, REAL y, REAL z);
   
   //CLeaf * getLigLeaf() const 
   //{  return m_pLeaf; }
   CLeaf* getLink(int i) const
   { return m_lLeaf[i]; }

   REAL get_Xcoors (int z)
   { return xcoors[z]; }
   REAL get_Ycoors (int z)
   { return ycoors[z]; }
   REAL get_Zcoors (int z)
   { return zcoors[z]; }
   vector<int > get_res_index() const
   { return res_index; }
   COORDS get_pxs() const
   { return pxs; }
   COORDS get_cxs() const
   { return cxs; }

private:

   void load_lig_coors(const char * fname, vector<REAL> & xcoors, vector<REAL > & ycoors, vector<REAL> & zcoors);
   CLeaf * m_pLeaf;
   vector<CLeaf *> m_lLeaf;
   REAL m_translate[3];
   REAL m_rotate[3][3];
   REAL m_pointx;
   REAL m_pointy;
   REAL m_pointz;
   vector<REAL> xcoors, ycoors, zcoors;
   vector<int > res_index;
   REAL (*pxs)[3];
   REAL (*cxs)[3];
   //vector<lpose *> accept_coords;

};

/*
void AssignLigandProperties(ligand *);
void AssignElements(ligand *);
void AssignHybrid(ligand *);
void AssignBondTypesAllH(ligand *);
void CopyUnmarkBonds(ligand *, int);
void ScanNextAtom(ligand *, int, int, int);
void MarkAllScannedBondsCyclic(ligand *);
void AssignRings(ligand *);
void MarkNextCyclic(ligand *, int, int);
ring* AssignFusedRings(ligand *, ring *);
void ScanNextRing(ligand *, int, int, int, float *);
void MarkScannedRingBonds(ligand *);
void UnmarkRingBonds(ligand *, int);
void RingNormals(ligand *);
void DetermineCMAtom(ligand *);
void AssignDonorAcceptor(ligand *);
void AssignForceFieldTypes(ligand *);
void AssignTorsions(ligand *);
void AssignTorsionParameters(ligand *);
void CreateTorsionConnectTrees(ligand *);
void ConnectedAtoms(atomList *, atom *);
void AssignRingtoTorsionTrees(ligand *);
void DetermineGlobalProbability4MC(ligand *);
void InitializeLonePairVectors(ligand *);
void LonePairVectors(ligand *);
void AssignScoringType(ligand *);
void AssignPairList(ligand *);
void AssignVinaAtomicRadius(ligand *);
void AssignVinaTypeLigand(ligand *);
*/
#endif
