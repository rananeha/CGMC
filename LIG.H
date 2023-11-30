//********************************************************
// Ligand parameters (also serve as leaf node types).
//********************************************************
#ifndef _LIG_H_
#define _LIG_H_

#define LIG 0
#define NUMLIG 1

// LIGAND

const int LIG_size = 13;
const int LIG_nGroups = 2;
const int LIG_nChis = 0;

char * LIG_aNames[LIG_size] = {"CAA", "CAB", "CAC", "CAD", "CAE", "CAF", "CAG", "CAH", "CAI", "CAJ", "CAL", "CAM", "NAK"};
REAL LIG_charges[LIG_size] = {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
//int LIG_groups[LIG_nGroups] = {0,1};  //{0,0};
int LIG_groups[LIG_nGroups] = {4,11};
int LIG_aTypes[LIG_nGroups] = {3,3};  // Ligand functional group types. Note: one cluster point can have more than types of function
int * LIG_chiTypes = NULL;

#endif
