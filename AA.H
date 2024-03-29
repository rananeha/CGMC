//********************************************************
// Amino Acids parameters (also serve as leaf node types).
//********************************************************

// Arginine
const int ARG_size = 13;
const int ARG_nGroups = 4;
const int ARG_nChis = 4;

int ARG_aTypes[ARG_size] = {CH1E, CH2E, CH2E, CH2E, NH1, H, CR, NC2, HC, HC, 
			    NC2, HC, HC };
char * ARG_aNames[ARG_size] = {"CA", "CB", "CG", "CD", "NE", "HE", "CZ", "NH1",
			       "1HH1", "2HH1", "NH2", "1HH2", "2HH2"};
REAL ARG_charges[ARG_size] = {0.10, 0.00, 0.00, 0.10, -0.40, 0.15, 0.25, -0.85,
			      0.40, 0.40, -0.85, 0.40, 0.40};  
int ARG_groups[ARG_nGroups] = {1, 3, 7,10};
int ARG_chiTypes[ARG_nChis] = {2, 2, 2, 3};

// Asparagine
const int ASN_size = 7;
const int ASN_nGroups = 3;
const int ASN_nChis = 2;

int ASN_aTypes[ASN_size] = {CH1E, CH2E, C, O, NH2, H, H};
char * ASN_aNames[ASN_size] = {"CA", "CB", "CG", "OD1", "ND2", "1HD2", "2HD2"};
REAL ASN_charges[ASN_size] = {0.10, 0.00, 0.55, -0.55, -0.60, 0.30, 0.30};  
int ASN_groups[ASN_nGroups] = {1, 2, 4};
int ASN_chiTypes[ASN_nChis] = {2, 4};

// Asparatic Acid
const int ASP_size = 5;
const int ASP_nGroups = 1;
const int ASP_nChis = 2;

int ASP_aTypes[ASP_size] = {CH1E, CH2E, C, OC, OC};
char * ASP_aNames[ASP_size] = {"CA", "CB", "CG", "OD1", "OD2"};
REAL ASP_charges[ASP_size] = {0.10, -0.15, 1.35, -0.60, -0.60};  
int ASP_groups[ASP_nGroups] = {1};
int ASP_chiTypes[ASP_nChis] = {2, 4};

// Cysteine
const int CYS_size = 3;
const int CYS_nGroups = 1;
const int CYS_nChis = 1;

int CYS_aTypes[CYS_size] = {CH1E, CH2E, SH1E};
char * CYS_aNames[CYS_size] = {"CA", "CB", "SG"};
REAL CYS_charges[CYS_size] = {0.10, 0.19, -0.19};  
int CYS_groups[CYS_nGroups] = {1};
int CYS_chiTypes[CYS_nChis] = {2};

// Glutamine
const int GLN_size = 8;
const int GLN_nGroups = 3;
const int GLN_nChis = 3;

int GLN_aTypes[GLN_size] = {CH1E, CH2E, CH2E, C, O, NH2, H, H};
char * GLN_aNames[GLN_size] = {"CA", "CB", "CG", "CD", "OE1", "NE2", "1HE2",
			       "2HE2"};
REAL GLN_charges[GLN_size] = {0.10, 0.00, 0.00, 0.55, -0.55, -0.60, 0.30, 
			      0.30};  
int GLN_groups[GLN_nGroups] = {1, 4, 5};
int GLN_chiTypes[GLN_nChis] = {2, 2, 4};

// Glutamic Acid
const int GLU_size = 6;
const int GLU_nGroups = 2;
const int GLU_nChis = 3;

int GLU_aTypes[GLU_size] = {CH1E, CH2E, CH2E, C, OC, OC};
char * GLU_aNames[GLU_size] = {"CA", "CB", "CG", "CD", "OE1", "OE2"};
REAL GLU_charges[GLU_size] = {0.10, 0.00, -0.15, 1.35, -0.60, -0.60};  
int GLU_groups[GLU_nGroups] = {1, 2};
int GLU_chiTypes[GLU_nChis] = {2, 2, 4};

// Histidine
const int HIS_size = 8;
const int HIS_nGroups = 3;
const int HIS_nChis = 2;

int HIS_aTypes[HIS_size] = {CH1E, CH2E, CR, NH1, H, CR1E, NR, CR1E};
char * HIS_aNames[HIS_size] = {"CA", "CB", "CG", "ND1", "CD2", "NE2",
			       "CE1", "HE2"};
REAL HIS_charges[HIS_size] = {0.10, 0.00, 0.10, -0.40, 0.30, 0.10, -0.40, 
			      0.30};  
int HIS_groups[HIS_nGroups] = {1, 2, 5};
int HIS_chiTypes[HIS_nChis] = {2, 4};

// Isoleucine
const int ILE_size = 5;
const int ILE_nGroups = 2;
const int ILE_nChis = 2;

int ILE_aTypes[ILE_size] = {CH1E, CH2E, CH3E, CH2E, CH3E};
char * ILE_aNames[ILE_size] = {"CA", "CB", "CG2", "CG1", "CD1"};
REAL ILE_charges[ILE_size] = {0.10, 0.00, 0.00, 0.00, 0.00};  
int ILE_groups[ILE_nGroups] = {1,3};
int ILE_chiTypes[ILE_nChis] = {2, 2};

// Leucine
const int LEU_size = 5;
const int LEU_nGroups = 1;
const int LEU_nChis = 2; 

int LEU_aTypes[LEU_size] = {CH1E, CH2E, CH1E, CH3E, CH3E};
char * LEU_aNames[LEU_size] = {"CA", "CB", "CG", "CD1", "CD2"};
REAL LEU_charges[LEU_size] = {0.10, 0.00, 0.00, 0.00, 0.00};  
int LEU_groups[LEU_nGroups] = {1};
int LEU_chiTypes[LEU_nChis] = {2, 2};

// Lysine
const int LYS_size = 9;
const int LYS_nGroups = 2;
const int LYS_nChis = 4;

int LYS_aTypes[LYS_size] = {CH1E, CH2E, CH2E, CH2E, CH2E, NH3, HC, HC, HC};
char * LYS_aNames[LYS_size] = {"CA", "CB", "CG", "CD", "CE", "NZ", "1HZ",
			       "2HZ", "3HZ"};
REAL LYS_charges[LYS_size] = {0.10, 0.00, 0.00, 0.00, 0.00, -1.35, 0.45, 
			      0.45, 0.45};  
int LYS_groups[LYS_nGroups] = {1, 4};
int LYS_chiTypes[LYS_nChis] = {2, 2, 2, 2};

// Methionine
const int MET_size = 5;
const int MET_nGroups = 2;
const int MET_nChis = 3;

int MET_aTypes[MET_size] = {CH1E, CH2E, CH2E, S, CH3E};
char * MET_aNames[MET_size] = {"CA", "CB", "CG", "SD", "CE"};
REAL MET_charges[MET_size] = {0.10, 0.00, 0.06, -0.12, 0.06};  
int MET_groups[MET_nGroups] = {1, 2};
int MET_chiTypes[MET_nChis] = {2, 2, 5};

// Phenylalanine
const int PHE_size = 8;
const int PHE_nGroups = 2;
const int PHE_nChis = 2;

int PHE_aTypes[PHE_size] = {CH1E, CH2E, CR, CR1E, CR1E, CR1E, CR1E, CR1E};
char * PHE_aNames[PHE_size] = {"CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2",
			       "CZ"};
REAL PHE_charges[PHE_size] = {0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}; 
int PHE_groups[PHE_nGroups] = {1, 4};
int PHE_chiTypes[PHE_nChis] = {2, 4};

// Proline
const int PRO_size = 4;
const int PRO_nGroups = 1;
const int PRO_nChis = 2;

int PRO_aTypes[PRO_size] = {CH1E, CH2E, CH2E, CH2E};
char * PRO_aNames[PRO_size] = {"CA", "CB", "CG", "CD"};
REAL PRO_charges[PRO_size] = {0.10, 0.00, 0.00, 0.1};  
int PRO_groups[PRO_nGroups] = {1};
int PRO_chiTypes[PRO_nChis] = {2, 2};

// Serine
const int SER_size = 4;
const int SER_nGroups = 1;
const int SER_nChis = 1;

int SER_aTypes[SER_size] = {CH1E, CH2E, OH1, H};
char * SER_aNames[SER_size] = {"CA", "CB", "OG", "HG"};
REAL SER_charges[SER_size] = {0.10, 0.25, -0.65, 0.40};  
int SER_groups[SER_nGroups] = {1};
int SER_chiTypes[SER_nChis] = {2};

// Threonine
const int THR_size = 5;
const int THR_nGroups = 2;
const int THR_nChis = 1;

int THR_aTypes[THR_size] = {CH1E, CH2E, OH1, H, CH3E};
char * THR_aNames[THR_size] = {"CA", "CB", "OG1", "1HG", "CG2"};
REAL THR_charges[THR_size] = {0.10, 0.25, -0.65, 0.40, 0.00};  
int THR_groups[THR_nGroups] = {0, 3};
int THR_chiTypes[THR_nChis] = {2};

// Tryptophan
const int TRP_size = 12;
const int TRP_nGroups = 4;
const int TRP_nChis = 2;

int TRP_aTypes[TRP_size] = {CH1E, CH2E, CR, CR, CR, CR1E, CR1E, NH1, H, CR1E,
			    CR1E, CR1E};
char * TRP_aNames[TRP_size] = {"CA", "CB", "CG", "CD2", "CE2", "CE3", "CD1",
			       "NE1", "HE1", "CZ2", "CZ3", "CH2"};
REAL TRP_charges[TRP_size] = {0.10, 0.00, -0.03, 0.10, -0.04, -0.03, 0.06, 
			      -0.36, 0.30, 0.00, 0.00, 0.00};  
int TRP_groups[TRP_nGroups] = {1, 2, 6, 9};
int TRP_chiTypes[TRP_nChis] = {2, 4};

// Tyrosine
const int TYR_size = 10;
const int TYR_nGroups = 4;
const int TYR_nChis = 2;

int TYR_aTypes[TYR_size] = {CH1E, CH2E, CR, CR1E, CR1E, CR1E, CR1E, C, OH1, H};
char * TYR_aNames[TYR_size] = {"CA", "CB", "CG", "CD1", "CE1", "CD2", "CE2",
			       "CZ", "OH", "HH"};
REAL TYR_charges[TYR_size] = {0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -0.25,
			      -0.65, 0.40};  
int TYR_groups[TYR_nGroups] = {1, 3, 5, 7};
int TYR_chiTypes[TYR_nChis] = {2, 4};

// Valine
const int VAL_size = 4;
const int VAL_nGroups = 1;
const int VAL_nChis = 1;

int VAL_aTypes[VAL_size] = {CH1E, CH2E, CH3E, CH3E};
char * VAL_aNames[VAL_size] = {"CA", "CB", "CG1", "CG2"};
REAL VAL_charges[VAL_size] = {0.10, 0.00, 0.00, 0.00};  
int VAL_groups[VAL_nGroups] = {1};
int VAL_chiTypes[VAL_nChis] = {2};

// Glycine (the only atom is a dummy atom. no sidechain!!!)
const int GLY_size = 1;
const int GLY_nGroups = 1;
const int GLY_nChis = 0;

int GLY_aTypes[GLY_size] = {CH1E};
char * GLY_aNames[GLY_size] = {"CA"};
REAL GLY_charges[GLY_size] = {0.10};
int GLY_groups[GLY_nGroups] = {0};
int * GLY_chiTypes = NULL;

// Alanine
const int ALA_size = 2;
const int ALA_nGroups = 1;
const int ALA_nChis = 0;

int ALA_aTypes[ALA_size] = {CH1E, CH3E};
char * ALA_aNames[ALA_size] = {"CA", "CB"};
REAL ALA_charges[ALA_size] = {0.10, 0.00};
int ALA_groups[ALA_nGroups] = {0};
int * ALA_chiTypes = NULL;

// N-terminal
const int NTR_size = 4;
const int NTR_nGroups = 1;

int NTR_aTypes[NTR_size] = {NH3, HC, HC, HC};
char * NTR_aNames[NTR_size] = {"N", "H1", "2H", "3H"};
REAL NTR_charges[NTR_size] = {-1.35, 0.45, 0.45, 0.45}; 
int NTR_groups[NTR_nGroups] = {0};
//const int NTR_size = 2;
//const int NTR_nGroups = 1;

//int NTR_aTypes[NTR_size] = {NH3, CH1E};
//char * NTR_aNames[NTR_size] = {"N", "CA"};
//REAL NTR_charges[NTR_size] = {-1.35, 0.25}; 
//int NTR_groups[NTR_nGroups] = {0};

// C-terminal
const int CTR_size = 3;
const int CTR_nGroups = 1;

int CTR_aTypes[CTR_size] = {C, OC, OC};
char * CTR_aNames[CTR_size] = {"C", "O", "OXT"};
REAL CTR_charges[CTR_size] = {1.20, -0.60, -0.60}; 
int CTR_groups[CTR_nGroups] = {0};
//const int CTR_size = 1;
//const int CTR_nGroups = 1;
//int CTR_aTypes[CTR_size] = {C};
//char * CTR_aNames[CTR_size] = {"C"};
//REAL CTR_charges[CTR_size] = {1.20}; 
//int CTR_groups[CTR_nGroups] = {0};

// Backbone -C, N, Ca
const int BBN_size = 4;
const int BBN_nGroups = 1;

int BBN_aTypes[BBN_size] = {C, O, NH1, H};
char * BBN_aNames[BBN_size] = {"C", "O", "N", "H"};
REAL BBN_charges[BBN_size+1] = {0.55, -0.55, -0.35, 0.25};
int BBN_groups[BBN_nGroups] = {0};
//const int BBN_size = 2;
//const int BBN_nGroups = 2;


// Backbone -C, N for Proline only (no H, but with CD)
const int BBP_size = 3;
const int BBP_nGroups = 1;

int BBP_aTypes[BBP_size] = {C, O, N};
char * BBP_aNames[BBP_size] = {"C", "O", "N"};
REAL BBP_charges[BBP_size] = {0.55, -0.55, -0.2};  
int BBP_groups[BBP_nGroups] = {0};
//const int BBP_size = 2;
//const int BBP_nGroups = 2;


// Backbone -C, N for Glycine only (the CA is a diferent atom type)
const int BBG_size = 2;
const int BBG_nGroups = 1;

int BBG_aTypes[BBG_size] = {NH1, H};
char * BBG_aNames[BBG_size] = {"N", "H"};
REAL BBG_charges[BBG_size+1] = {-0.35, 0.25};
int BBG_groups[BBG_nGroups] = {0};
//const int BBG_size = 2;
//const int BBG_nGroups = 2;

// Backbone -C, N for Glycine only (the CA is a diferent atom type)
const int BCN_size = 2;
const int BCN_nGroups = 1;

int BCN_aTypes[BCN_size] = {C, O};
char * BCN_aNames[BCN_size] = {"C", "O"};
REAL BCN_charges[BCN_size+1] = {0.55, -0.55};
int BCN_groups[BCN_nGroups] = {0};
