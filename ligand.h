
#ifndef _LIGAND_H
#define _LIGAND_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utilities.h"
#include "parameters.h"
#include "global.h"

typedef struct LPList
{
	float				x[3];
	int					atnA;	//atom number of Donor
	int					num;
	struct LPList		*next;
} LPList;

typedef struct conformation
{
    float           x[3];
    //float           charge;
    LPList          *lplist;

} conformation;

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

typedef struct sybylbond
{
    int            at1;
    int            at2;
    char           bondtype[10];
} sybylbond;

typedef struct LPE
{
    float           x[3];		//coordinate	
	float			v[3];
} LPE;

typedef struct LEdge			
//edge formed by two ligand pharmacophores, please noticing the difference between this with the ligand-prot pharmacophore pairs (in the input ligand mol2 file, this edge are called pair)
{
    int             type;
    float			mindist, maxdist;
} LEdge;

typedef struct atomList
{
   atom		 			*ai;
   struct atomList		*next;
} atomList;

typedef struct ring	// structures for rings and fused rings
{
	atomList		 	*list;
	conformation		*center;
	conformation		**normal;
	conformation        **cen_sol;
	conformation        ***norm_sol;
	conformation		*cen_MC;
	conformation        **norm_MC;
	int					numatom;
	struct ring			*next;
} ring;

typedef struct Water
{
	int				id;
	float			Ocoor[3];
	float			Hcoor[2][3];
	float			LPcoor[2][3];
	float			occupancy;
	int				on;
} Water;

typedef struct Pose
{
	int				id;
	float			score;
	float			WaterScore;
	float			**coor;
	Water			*wtr;
	int				numMediatingWater;

} Pose;


typedef struct MC_sol
{
    float			stability;
    char			nearest_group[8];
    int				nearest_num;
    int				good;
    int				same;
    float			pharm_score;	//score before MC search
    float			LigIntEn0, LigHDA0, LigHYDR0, LigVDW0, LigRing0, LigElec0, LigTor0, LigNonBond0;// starting score
	float			MC_score;	// minimum score after 1st MC search
	float			MC_score2; 	// minimum score after 2nd MC search
	int				Numcf;		// number of conformations contributing to calculation of entropy
	float			HDAscore, HYDRscore, VDWscore, Ringscore, Elecscore, Metascore;
	float			LigIntEn, LigHDA, LigHYDR, LigVDW, LigRing, LigElec, LigTor, LigNonBond;
    float			Gibbs;   	//score with entropy included
    int				prot_conf;

    float           **lpe_coor;
    int             *lpe_type;
    int             numlpe;
    int             *lpe_list;

    float			rmsd_MC;   	//rmsd after MC refinement
    float			dist_exp, dist_near;
    float			Prob; 		//metabolic probability
} MC_sol;

typedef struct torsion
{
	atom				*at1;
	atom				*at2;
	atom				*at3;
	atom				*at4;
	atom				*atm_move_end;		//move end of the torsion 
	atom				*atm_fix_end;		//fixed end of the torsion
    int					NumMovingAtom;	//number of atoms in the moving part

	float				nonbonded_gradient;		// gradient for the nonbonded interaction
	float				bonded_gradient;		// gradient for the bonded term (torsion energy)
	atomList			*tree0;			//moving side of the torsion
//	atomList			*tree1;
//	atomList			*treelib;
	float				Phase[4];
	float				Multi[4];
	float				Depth[4];
	float				k[4];			// force constant
	int					isDone;
	int					numTorsionForBond[4];
	float				actval, Xrayval;
    int					num;
    int					FFnum;
    int					MaxMultiNum;	//number of the Maximum multiplicity
    int					amberindex;		//when reading prot.top file for dihedral...

	int					numrings;		//number of rings in the moving part of the torsion (tree0)
	int					*ringindex;		//ring index for the rings in the moving part
	float				local_tors_size;
	float				global_tors_prob;
	struct torsion		*next;
} torsion;

typedef struct InternalPairList
{
	int				atn1;
	int				atn2;
	float			Prefact;
	
	struct InternalPairList		*next;
} InternalPairList;

typedef struct HBList
{
	int				atnD;	// donor
	int				atnA;	// acceptor
	int				atnH;	// donor-H
	int				LPnum;	// lone pair number
	float			Prefact;
	
	struct HBList		*next;
} HBList;

typedef struct ligand
{
    int             numconf;
    int		    *numsols;
    int             numatoms;
    int				numhatms;	//number of heavy atoms
    int             numlps;  // number of total lone pairs of the ligand
    int             numrings;
    atom            *atm;
	int				CA;		// Centroid atom, count from 0

	int				numsymmaps; //number of symm atom pairs
	int				numsymmgrs; //number of symm groups
	int				numsymmtri;
	int				numsymmquo;
	int				**symmtrilist;
	int				**symmquolist;

	float			TotalCharge;
	float			MW;			//molecular weight
    int             res_no;
    char            res_name[20];
    char			meta_group[10][8];
    int				num_metas; // number of meta groups
    int				meta_num[10];  // atom number of the meta_group (count from 0)

	int				lenPL;			//length of the pairlist 
	int				*PairList;		//nonbonded pairlist
	int				*PairList_0;	//starting index of the pairlist for atom i
	int				*PairList_1;	//ending index of the pairlist for atom i
	float			*C6_List;		//C6 parameter list for Lennard-Jones calc
	float			*C12_List;		//C12 parameter list for Lennard-Jones calc
	float			*qq_List;		//charge-charge list for elec calc

    int             numbonds;		//number of bonds
	int				*bond_i;		//atom involved in bond "i"
	int				*bond_j;		//atom involved in bond "i"
	int				*bond_amber;	//index into parameter arrays RK and REQ
	int				*bond_withH;	//whether the bond involves hydrogen (1=yes, 0=no)
	float			*bond_Fconst;
	float			*bond_l0;

	int				numAngles;
	int				*angle_i;
	int				*angle_j;
	int				*angle_k;
	int				*angle_amber;	//index into parameter arrays TK and TEQ for angle
	float			*angle_Fconst;	
	float			*angle_t0;

	int				numTor;
	int				numimTor;

	int				numDihedrals;
	int				*dih_i, *dih_j, *dih_k, *dih_l;
	int				*dih_amber;
	float			*dih_Depth, *dih_Multi, *dih_Phase;
	
	int				numImpropers;
	int				*imp_i, *imp_j, *imp_k, *imp_l;
	int				*imp_amber;
	float			*imp_Depth, *imp_Multi, *imp_Phase;

	float			unbound_Int_E[8];	// minimum internal energy among unbound form ligands
	float			*omega2Energy;
    int             num_Ligand_Pharmacophores;
    LPE             **Ligand_Pharmacophore;	
	int				*Ligand_Pharmacophore_type;
	int				*Ligand_Pharmacophore_heavyatom;
	float			*Ligand_Pharmacophore_scale;
	int				**Ligand_Pharmacophore_H;
	int				num_hda_pharms;
	int				pharms4rmsd[4];

	LEdge			**Ligand_Edges;   
	int				num_pharm_type[num_lig_pharm_type];	//number of pharmacophores for each type
	
    ring            *Ring0;
	torsion			*Torsion0;
	torsion			*ImproperTorsion0;
	InternalPairList	*PairList0;
	HBList			*PairListDA0;

	int				*degree_counter;
	float			**LSposes;//uses ligandscout scoring function
	float			*LSscore;
	float			**LSscore_component;

	Pose			*pose;
	int				NumMCSols;		// number of solutions get from combination of different conformations sols fot MC
    MC_sol			*MCsol;
    float			*optimized_score;
    float			*G_score;
	int				NumSols;
    float			xrayscore[13];
    
    float			Rv[3];			//rotation vector in internal coordinate minimizer
	float			RotGradient;	//rotation gradient in internal coordinate minimizer
	float			TransForce[3];

} ligand;


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
#endif
