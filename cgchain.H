#ifndef _CGCHAIN_H
#define _CGCHAIN_H

#include <vector>
#include "cgpdb.H"
#include "cgleaf.H"

#define ANGS_FILE 1
#define STARTING_PT = {"0.292","0.005","1.433"}
struct REAL3 {
  REAL v[3];
};

struct REAL33 {
  REAL m[3][3];
};

typedef pair<REAL3, REAL33> FRAME;

//class CBV;
class CRss;
class CSkipList;

class CChain
{
public:
  CChain(const char * fname, const char * fnameAX, const char * fname3D, const char * fn_quadB, const char * fn_quadE, const char * fn_quadI, const char * fn_tripB, const char * fn_tripE, const char * fn_tripI, const char * fn_sing, const char * fl_quadB, const char * fl_quadE, const char * fl_quadI, const char * fl_tripB, const char * fl_tripE, const char * fl_tripI, const char * fl_sing, int type);

  void computeLinkCoords();
  void readjustLinkCoords();
  void refillLinkCoords();
  void store_orig_Coordinates(const char * fname, const char * fname3D);
  void storeCoordinates(const char * fname);
  void load_quad_buried_prob(const char * fn_quadB, vector<REAL> & quadB_prob);
  void load_quad_exposed_prob(const char * fn_quadE, vector<REAL> & quadE_prob);
  void load_quad_intermediate_prob(const char * fn_quadI, vector<REAL> & quadI_prob);
  void load_trip_buried_prob(const char * fn_tripB, vector<REAL> & tripB_prob);
  void load_trip_exposed_prob(const char * fn_tripE, vector<REAL> & tripE_prob);
  void load_trip_intermediate_prob(const char * fn_tripI, vector<REAL> & tripI_prob);
  void load_sing_prob(const char * fn_sing, vector<REAL> & sing_prob);
  void load_quadlig_buried_prob(const char * fl_quadB, vector<REAL> & quadligB_prob);
  void load_quadlig_exposed_prob(const char * fl_quadE, vector<REAL> & quadligE_prob);
  void load_quadlig_intermediate_prob(const char * fl_quadI, vector<REAL> & quadligI_prob);
  void load_triplig_buried_prob(const char * fl_tripB, vector<REAL> & tripligB_prob);
  void load_triplig_exposed_prob(const char * fl_tripE, vector<REAL> & tripligE_prob);
  void load_triplig_intermediate_prob(const char * fl_tripI, vector<REAL> & tripligI_prob);
  void load_singlig_prob(const char * fl_sing, vector<REAL> & singlig_prob);

  void load_coors(const char * fname3D, vector<REAL> & coorsx, vector<REAL> & coorsy, vector<REAL> & coorsz);
  void storeAngsStyle(const char * fname);
  void adjustBV();
  //REAL computeTorsionEnergy();
  //REAL updateTorsionEnergy(const vector<ANGLE_CHANGE> & angles);

  // Undo the change to the torsion energy caused by the latest move.
  void undoTorsionChange()
  { m_torsionE = m_undoTorsionE; }

  void getCaPositions(POSITIONS & pos);
  
  int getLength() const
  { return m_links.size(); }
  const vector<FRAME> & getFrames() const
  { return m_frames; }
  const vector<CLeaf*> & getLinks() const
  { return m_links; }
  CLeaf* getLink(int i) const
  { return m_links[i]; }
  //REAL getTorsionEnergy() const
  //{ return m_torsionE; }
  //REAL getSolventE() const
  //{ return m_solventE; }
  int getChainLength() const
  { return m_links.size()/2.0-1; }

  vector<REAL> quadB_prob;
  vector<REAL> quadE_prob;
  vector<REAL> quadI_prob;
  vector<REAL> tripB_prob;
  vector<REAL> tripE_prob;
  vector<REAL> tripI_prob;
  vector<REAL> sing_prob;
  vector<REAL> quadligB_prob;
  vector<REAL> quadligE_prob;
  vector<REAL> quadligI_prob;
  vector<REAL> tripligB_prob;
  vector<REAL> tripligE_prob;
  vector<REAL> tripligI_prob;
  vector<REAL> singlig_prob;

private:
  void create(const vector<int> AAtypes, const vector<REAL> & phis, 
  	      const vector<REAL> & psis, const vector<int> & rots,
              const vector<REAL> & axes_x, const vector<REAL> & axes_y, const vector<REAL> & axes_z);
              //const vector<REAL> & coorsx, const vector<REAL> & coorsy, const vector<REAL> & coorsz);
  void load_angs(const char * fname, vector<int> & AAtypes, vector<REAL> & phis, vector<REAL> & psis, vector<int> & rots);
  void load_axes(const char * fnameAX, vector<REAL> & axes_x, vector<REAL> & axes_y, vector<REAL> & axes_z);
  //void load_coors(const char * fname3D, vector<REAL> & coorsx, vector<REAL> & coorsy, vector<REAL> & coorsz);
  //REAL computeSolventRefEnergy();
  void computeCaPositions(REAL Calphas[][3]);

  vector<FRAME> m_frames;
  vector<CLeaf*> m_links;
  REAL m_torsionE, m_undoTorsionE;
  //REAL m_solventE;
};

#endif
