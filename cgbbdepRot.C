#include <iostream>
#include <fstream>

#include "rss.H"
#include "cgbbdepRot.H"

using namespace std;

#define BB_DEP_ROT_LIB "rotamer_lib"

char AA_NAMES[NUMAA][5]={"ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "GLY", "ALA", "NTR", "CTR", "BBN", "BBP","BBG"};
int NUM_ROTAMERS[NUMAA] = {81, 9, 9, 3, 27, 27, 6, 9, 9, 81, 27, 6, 2, 3, 3, 9, 6, 3, 1, 1, 1, 1, 1, 1, 1};

void SIDECHAIN::createSidechains(const char * dir)
{
  ifstream fin("%s/%s", dir ,BB_DEP_ROT_LIB);
  if (fin.fail())
  {
    cerr << "Unable to open rotamer_coords\n";
    exit(8);
  }

  
