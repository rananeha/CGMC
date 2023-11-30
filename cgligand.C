/***************************************************************************
 *   Copyright (C) 2010 by xu57   *
 *   xu57@gorilla.mcmp.purdue.edu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdlib>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <string>
#include <sstream>
#include <set>
#include <cmath>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include "dataType.h"
//#include "memory_usage.h"

#ifndef _cplusplus
#define _cplusplus
//#include "limocscore.h"
#endif

#define MAXINDEX 1000
#define SCALAR_FACTOR 10
//#define _LIMOC_DEBUG

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstdlib>
#include <vector>
#include "cgeef1.H"
#include "cgchain.H"
#include "cgligand.H"
#include "cgpdb.H"

using namespace std;

//inline float *FindClusterCord(int ** const clusters, int ** const cluster_pos, const ClusterT * const pCl, const int &NumRes, const ResCA * const resca, int &atom_index, float ** const prot_cord);

//inline float *FindClusterMemberCord(int ** const clusters, int ** const cluster_pos, const ClusterT * const pCl, const int &NumRes, const ResCA * const resca, int &atom_index, float ** const prot_cord, const int n);
//=================================================================================================
// Monte Carlo Ligand class
//=================================================================================================
/*
MonteCarloLigand::MonteCarloLigand(const char * fname, REAL pointx, REAL pointy, REAL pointz)
{
   //vector<REAL> xcoors, ycoors, zcoors;
   load_lig_coors(fname, xcoors, ycoors, zcoors);
}
*/
//=================================================================================================
// Load Ligand conformer
//=================================================================================================
void MonteCarloLigand::load_lig_coors(const char * fname, vector<REAL> & xcoors, vector<REAL> & ycoors, vector<REAL> & zcoors)
{

  ifstream fin(fname);
  if (!fin.is_open())
    {
      cout << "Could not open coordinate file " << fname << endl;
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


        xcoors.push_back(cx);
        ycoors.push_back(cy);
        zcoors.push_back(cz);
      }
   fin.close();
   cout << "loaded  lig" << endl;

}
/*
//=================================================================================================
// output ligand Amber crd
//=================================================================================================
void WriteCoordinates(int num_flex_atoms, float **prot_cords, ligand *lig, float **pxs, FILE *fout)
{
    int		         i, j, n;
    //int				na, nb, nc, nr, bnum;
    //bondList		*pb;

    n = 0;
    for(i = 0; i < num_flex_atoms; i++)
    {
        for(j = 0; j < 3; j++)
        {
            if(n >= 10)
            {
                fprintf(fout, "\n");
                n = 0;
            }
            if(isnan(prot_cords[i][j]))
            {
                printf("Coordinates are nan.\n");
                exit(0);
            }
            fprintf(fout, "% 8.3f", prot_cords[i][j]);
            n++;
        }
    }
    for(i = 0; i < lig->numatoms; i++)
    {
        for(j = 0; j < 3; j++)
        {
            if(n >= 10)
            {
                fprintf(fout, "\n");
                n = 0;
            }
            if(isnan(pxs[i][j]))
            {
                printf("Coordinates are nan.\n");
                exit(0);
            }
            fprintf(fout, "% 8.3f", pxs[i][j]);
            n++;
        }
    }
    if(n >= 1)
    {
        fprintf(fout, "\n");
    }
    //if(par->BoxOn == 1)
    //{
        //fprintf(fout, "% 8.3f%8.3f%8.3f\n", BoxX[0], BoxX[1], BoxX[2]);
    //}
    //else
    //{
        //		fprintf(fout, "% 8.3f%8.3f%8.3f\n", 0.0, 0.0, 0.0);
    //}
    fflush(fout);

}
*/

//=================================================================================================
// Write lig coordinates
//=================================================================================================
//void WriteLigCoordinates(ligand *lig, float **pxs, FILE *fout)
//void WriteLigCoordinates(float *pxs, FILE *fout)
void MonteCarloLigand::WriteLigCoordinates(vector<LPOSE> & accept_coords, int num_poses, char * ofile)
{
    int		         i, j, n, l;
    //int				na, nb, nc, nr, bnum;
    //bondList		*pb;
    ofstream fout;
    char buf[100];
    //sprintf(buf, "%s/%s.out", dir, ofile);
    sprintf(buf, "%s.out", ofile);
    fout.open(buf);

    if (!fout.is_open())
    {
      cout << "Could not open output file " << buf << endl;
      exit(1);
    }

  sprintf(buf,"HEADER    %-70s", ofile);
  fout << buf << endl;
  sprintf(buf, "%-80s", "COMPND");
  fout << buf << endl;
  sprintf(buf, "%-80s", "SOURCE");
  fout << buf << endl;

  int m;
  //vector<AA>::const_iterator aat = aas.begin();
  AA aa;
  ATOM_ at;

    //cout << "Reached here" << endl;
    int size = LIGAND::m_liglist[0]->m_size;
    double (*pxs)[3];
    pxs = new double[size][3];
    n = 0;

    for (l = 0; l < num_poses; l++)
    {
       //pxs = accept_coords[l].save_pxs;

       for(i = 0; i < size; i++)
       {
          strcpy(at.name, LIGAND::m_liglist[0]->m_aNames[i]);
          //cout << i << " " << at.name << endl;
          for(j = 0; j < 3; j++)
           {
            pxs[i][j] = accept_coords[l].save_pxs[i][j];
            //if(n >= 10)
            //{
                //fprintf(fout, "\n");
                //n = 0;
            //}
            at.pos[j] = pxs[i][j];
           }
            if(isnan(pxs[i][j]))
            {
                printf("Coordinates are nan.\n");
                exit(0);
            }

            //cout << "Reached here" << endl;
            //writeLine(fout, i++, aat->atoms[k].name, AA_NAMES[aat->type], " ", j, pxs[i][0], pxs[i][1], pxs[i][2]);
            writeLine(fout, n++, at.name, "LIG", " ", 1, at.pos[0], at.pos[1], at.pos[2]);

            //if (aat == aas.end()-1)
                //writeLine(fout, i++, aat->atoms[size-1].name, AA_NAMES[aat->type], " ", j, aat->atoms[size-1].pos[0], aat->atoms[size-1].pos[1], aat->atoms[size-1].pos[2]);

            //sprintf(buf,"TER   %5d%6s%3s  %4d%54s",i,"",AA_NAMES[(aat-1)->type],j-1," ");
            //sprintf(buf, "MASTER    %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%10s", 0,0,0,0,0,0,0,0,i-1,1,0,0," ");
            //fout << buf << endl;

            //sprintf(buf, "% 8.3f", pxs[i][j]);
            //fout << buf << endl;
            //fprintf(fout, "% 8.3f", pxs[i][j]);
            //n++;
        }
        sprintf(buf,"%-80s","TER");
        fout << buf << endl;
        sprintf(buf,"%-80s","END");
        fout << buf << endl;
        
     }
    cout << "Coordinates written" << endl;
    //if(n >= 1)
    //{
        //fprintf(fout, "\n");
    //}
    //if(par->BoxOn == 1)
    //{
        //fprintf(fout, "% 8.3f%8.3f%8.3f\n", BoxX[0], BoxX[1], BoxX[2]);
    //}
    //else
    //{
        //		fprintf(fout, "% 8.3f%8.3f%8.3f\n", 0.0, 0.0, 0.0);
    //}
    //fflush(fout);

}
//=================================================================================================
// output ligand mol2
//=================================================================================================
void MonteCarloLigand::writeLine(ofstream & fout, int index, const char * aname, const char * resname,
               const char * chainid, int resnum,
               REAL x, REAL y, REAL z)
{
  fout << "ATOM  ";

  fout.width(5);
  //fout << right;
  fout.setf(ios::right, ios::adjustfield);
  fout << index;

  fout << ' ';

  fout.setf(ios::left, ios::adjustfield);
  if (aname[0] > '0' && aname[0] <= '9')
    {
      fout.width(4);
      fout << aname;
    }
  else
    {
      fout << ' ';
      fout.width(3);
      fout << aname;
    }

  fout << ' ';

  fout << resname;

  fout << ' ';

  fout << chainid;

  fout.setf(ios::right, ios::adjustfield);
  fout.width(4);
  fout << resnum;

  fout << "    ";

  char buf[50];
  sprintf(buf, "%8.3f%8.3f%8.3f%26s", x,y,z," ");
  fout << buf << endl;;
}
/*
//=================================================================================================
// output ligand mol2
//=================================================================================================
void OutputLigand(ligand *lig, int lig_nr, float **pxs, FILE *fp){

	char line[512];
	FILE *fi;
	int i, atomnum, bondnum;

	fi = fopen(par.ListOfLigands[lig_nr],"r");
	fprintf(fp, "@<TRIPOS>MOLECULE\n****\n");
	atomnum = lig->numatoms;
	bondnum = lig->numbonds;
	fprintf(fp, "% 5d % 5d     0     0     0\n", atomnum, bondnum);
	fprintf(fp, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");
			
	for(i = 0; i < lig->numatoms; i++){
		fprintf(fp, "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 SUB       % 8.4f\n",
		lig->atm[i].atom_no, lig->atm[i].name,
		pxs[i][0], pxs[i][1], pxs[i][2],
		lig->atm[i].sybyl_type, lig->atm[i].charge);
	}
	fprintf(fp, "@<TRIPOS>BOND\n");
	while(!feof(fi)){
		fgets(line, lString, fi);
		if(strstr(line, "@<TRIPOS>BOND"))
			break;
	}
	for(i = 0; i < lig->numbonds; i++)
	{
		fgets(line, lString, fi);
		fprintf(fp, line);
	}	fclose(fi);   
}*/

//=================================================================================================
// DetermineInertiaTensor
//=================================================================================================
//void DetermineInertiaTensor(int size, float * cen, float **coor, double II[3][3])
void DetermineInertiaTensor(int size, double * cen, float (*coor)[3], double II[3][3])
{
	double	xx, yy, zz, xy, xz, yz;
	int		i, cm;
	xx = 0.0;
	yy = 0.0;
	zz = 0.0;
	xy = 0.0;
	xz = 0.0;
	yz = 0.0;
//	cm = lig->CA;

	for(i = 0; i < size; i++)
	{
		//if(lig->atm[i].element != eH)
		//{
			xx += pow((coor[i][0]-cen[0]),2);
			yy += pow((coor[i][1]-cen[1]),2);
			zz += pow((coor[i][2]-cen[2]),2);
			xy += (coor[i][0]-cen[0])*(coor[i][1]-cen[1]);
			xz += (coor[i][0]-cen[0])*(coor[i][2]-cen[2]);
			yz += (coor[i][1]-cen[1])*(coor[i][2]-cen[2]);
		//}
	}
	II[0][0] = yy + zz;
	II[1][1] = xx + zz;
	II[2][2] = yy + xx;
	II[0][1] = -xy;
	II[1][0] = -xy;
	II[0][2] = -xz;
	II[2][0] = -xz;
	II[1][2] = -yz;
	II[2][1] = -yz;
}
/*
//=================================================================================================
// monte_carlo
//=================================================================================================
 void MonteCarloLimoc(ligand * lig)
{
        int i;
	//par.MC_steps = par.MC_steps_1st;
	//par.MC_steps = 0;
        cout<<"Total number of poses need to be optimized: "<<lig->NumSols<<endl;
	//lig->NumSols = 100;

	//par.MC_steps = par.MC_steps_1st;

        //srand( time(NULL) );
	srand(4);	
	//srand(2);
	//srand(3);
	//sran();

	cout<<"--------------------------------------------"<<endl;
	//MC_Limoc(lig, s, prot, Gridmin, Griddelta, numGP, prot_cord, clusters, cluster_pos, resca, gridscore, NumGridScore, Comb, NumCombs, NumRes, sol_cl, 1, SolutionComb+s, outfile, lig_nr);

	
}
*/
//=================================================================================================
// monte_carlo_Limoc
//=================================================================================================

/***********************************************
* ----------  MC algorithm ------------------- *
* 1. perform move                              *
* 2. assign configuration                      *
* 3. calc score                                *
* 4. save new coordinates                      *
************************************************/
//void MonteCarloLigand::MC_Ligand(int pose, vector<CLeaf *> res_list)
//void MonteCarloLigand::MC_Ligand(CChain & m_chain, int pose, vector<int > res_index, REAL (*pxs)[3], REAL (*cxs)[3])
void MonteCarloLigand::MC_Ligand(CChain & m_chain, int pose)
{
	//float		**xs, **save_xs, **orig_xs, **min_xs;
        //float           (*xs)[3], (*orig_xs)[3];
        float           xs[LIGAND::m_liglist[0]->m_size][3];
        float           orig_xs[LIGAND::m_liglist[0]->m_size][3];
        float           * save_xs[3], (*min_xs)[3];
        float           *xa;
        double          cen[3], orig_cen[3];
	REAL		CM[3], Pre_CM[3], Rot[4], Pre_Rot[4];
	double		I[3][3], d[3], v[3][3];
	int			n_rot, ss=0;
	int			i, j, k, l, t, ii, st, m, numLP, a, lp, atm2, atm3;


//	float		av[3],bv[3],cv[3],dv[3],center[3];
    	float		c, s, uxx, uyy, uzz, uxy1, uxy2, uxz1, uxz2, uyz1, uyz2;
        REAL		R[3][3];
//    	float		min_local_CM, min_local_Rot,  min_local_Tors;
//	float		tmp_local_CM, tmp_local_Rot, tmp_local_Tors, sum;

	int         hatm, accept_counter_1 = 0, non_accepted_counter = 0, non_accepted_counterB = 0, rmsdmax_counter = 0, save_counter = 0; //YY add
	float		diff[3], y[3];
	float		fact, random_nr, rmsd_max, pe; 
        double          random, tev;
	//int			MMn;

	int			tt, rmsx, Numcf, cf, *Sflag, MINST = 0;
        const gsl_rng_type * T;
        gsl_rng * r;
        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        
        int ligsize = LIGAND::m_liglist[0]->m_nGroups;
        cout << ligsize << endl;
        pxs = new double[LIGAND::m_liglist[0]->m_size][3];
        cxs = new double[LIGAND::m_liglist[0]->m_nGroups][3];
        /*for (int i = 0; i < ligsize; i++)
           for (int j = 0; j < 3; j++)
           {
              cout << i << " " << j << endl;
              cxs[i][j] = 0.0;
           }*/
/*
	FILE		*fi, *fo, *fp, *fp2, *fp3;
	int     	atomnum, bondnum;
	char		line[kString];

	vector<pair<long double, int> >   accept_solutions;

	xa = (float*) calloc(3*lig->numatoms, sizeof(float));
*/
//        vector<float** > accept_cords;


//	FILE	*fj,*lig_trj,*lig_score;//YY add

	//lig_trj = fopen("ligand_trj.mol2","w");
	//lig_score = fopen("ligand_score","w");
	//fj = fopen(par.ListOfLigands[0],"r");


    /************************************************
    * start monte carlo
    ************************************************/


    //vector<REAL> xcoors, ycoors, zcoors;
    //load_lig_coors(fname, xcoors, ycoors, zcoors);
    //assign coordinates
    //xs = (float**)calloc(lig->numatoms,sizeof(float*));
    //save_xs = (float**)calloc(lig->numatoms,sizeof(float*));
    //orig_xs = (float**)calloc(lig->numatoms,sizeof(float*));
    //min_xs = (float**)calloc(lig->numatoms,sizeof(float*));


    // store starting coordinates
    //for (i=0;i<lig->numatoms;i++)
    //{
	/*
        xs[i]=(float*)calloc(3,sizeof(float));
        save_xs[i]=(float*)calloc(3,sizeof(float));
        orig_xs[i]=(float*)calloc(3,sizeof(float));
        min_xs[i]=(float*)calloc(3,sizeof(float));*/
        /*
        for (j=0;j<3;j++)
        {
            xs[i][j] = lig->atm[i].sol[solnum].x[j];
            save_xs[i][j] = lig->atm[i].sol[solnum].x[j];
            orig_xs[i][j] = lig->atm[i].sol[solnum].x[j];
            min_xs[i][j] = lig->atm[i].sol[solnum].x[j];
        }*/
    //}
    /*
    for(i = 0; i < lig->numatoms; i++)
    {
	xa[3*i] = xs[i][0];
	xa[3*i+1] = xs[i][1];
	xa[3*i+2] = xs[i][2];
    }
    */
    /*
    if(level == 1)
    {

	//YY add
	X0 = par.C1[0]-par.C0[0];
	Y0 = par.C1[1]-par.C0[1];
	Z0 = par.C1[2]-par.C0[2];
	A0 = sqrt( X0*X0 + Y0*Y0 + Z0*Z0 );

        float *cen1 = getCentroid(xa, lig);       
        X1 = cen1[0]-par.C0[0];
        Y1 = cen1[1]-par.C0[1];
        Z1 = cen1[2]-par.C0[2];
        A1 = (X1*X0+Y1*Y0+Z1*Z0) / sqrt( X0*X0 + Y0*Y0 + Z0*Z0 );

		//YY add
		fprintf(lig_trj, "\n\n@<TRIPOS>MOLECULE\n****\n");
		atomnum = lig->numatoms;
		bondnum = lig->numbonds;
		fprintf(lig_trj, "% 5d % 5d     0     0     0\n", atomnum, bondnum);
		fprintf(lig_trj, "SMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n"); //?

	
		for(i = 0; i < lig->numatoms; i++)
		{
			fprintf(lig_trj, "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 SUB       % 8.4f\n",
					lig->atm[i].atom_no, lig->atm[i].name,
					xa[3*i+0], xa[3*i+1], xa[3*i+2],
					lig->atm[i].sybyl_type, lig->atm[i].charge);
		}

		fprintf(lig_trj, "@<TRIPOS>BOND\n");
		while(!feof(fj))
		{
			fgets(line, lString, fj);
			if(strstr(line, "@<TRIPOS>BOND"))
				break;
		}
		for(i = 0; i < lig->numbonds; i++)
		{
			fgets(line, lString, fj);
			fprintf(lig_trj, line);
		}
		rewind(fj);
		//YY end



    }*/


    //DetermineInertiaTensor(lig, save_xs, I);
	
    //jacobi3(I, d, v, &n_rot);
    /*	
    for (i = 0; i < 3; i++)
    {
        Rot[i] = v[i][0];
        Pre_Rot[i] = v[i][0];
    }
    */
    // set starting configuration
    for (i = 0; i < 3; i++)
    {
        CM[i] = 0.0;
        Pre_CM[i] = 0.0;
        cen[i] = 0.0;
    }

    Rot[3] = 0.0;
    Pre_Rot[3] = 0.0;

    //set local start values
    //tmp_local_CM   = par.local_CM;
    //tmp_local_Rot  = par.local_Rot;
    //tmp_local_Tors = par.local_Tors;
    //min_local_CM = par.min_local_CM;
    //min_local_Rot = par.min_local_Rot;
    //min_local_Tors = par.min_local_Tors;

    tt = 0;
    rmsx = 0;
    Numcf = 1;

    //YY add -allocate a two dimentional array
    //float *cen2;
    //float *save_cen;
    //save_cen = (float *)calloc(par.MC_steps,sizeof(float));
    float max_move;
	
    // start monte carlo search
    int size = LIGAND::m_liglist[0]->m_size;
    for (st = 1; st < 2; st++) //par.MC_steps
    {
        // 1. Translation: Centroid of ligand // Gaussian(Normal) Distribution
        //int size = LIGAND::m_liglist[0]->m_size;
	//random = ((double) rand() / (RAND_MAX)) ;//YY random
        //random = gsl_rng_uniform (r);
        //gsl_rng_free(r);
        ///cout << "ligand size " << size << endl;
        //cout << "pose " << pose << endl;
        for (i = 0; i < size; i++)
          {
            for (j=0;j<3;j++)
            {
                xs[i][j] = LIGAND::m_liglist[0]->m_lconformer[pose]->m_positions[i][j]; //position of conformer in 0-0-0 space
                //cout << xs[i][j] << endl;
                orig_xs[i][j] = xs[i][j];
                //cout << orig_xs[i][j] << endl;
                //xs[i][j] =  0.0;
                //cout << LIGAND::m_liglist[0]->m_lconformer[pose]->m_positions[i][j] << endl;
                //cout << xs[i][j] << endl; 
            }
            cen[0] += xs[i][0];
            cen[1] += xs[i][1];
            cen[2] += xs[i][2];
          }
        cen[0] /= size;
        cen[1] /= size;
        cen[2] /= size;
        orig_cen[0] = cen[0];
        orig_cen[1] = cen[1];
        orig_cen[2] = cen[2];
        //cout << cen[0] << " " << cen[1] << " " << cen[2] << endl;
        //cen = getCentroid(xs, lig);
        //exit(0);
        
        //cout << "Pose" << pose << endl;
        for (i = 0; i < 3; i++)
        {
            //cout << (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) << endl;
            //int TMP_LOCAL_CM = gsl_ran_gaussian(r, SIGMAA);
            //cout << tmp << endl;
            //CM[i] = Pre_CM[i] + 2.0*tmp_local_CM*( ((double)rand() / (RAND_MAX)) - 0.5);
            //CM[i] = Pre_CM[i] + gsl_ran_gaussian(r, SIGMAA); // random translation in Gaussian distribution
            //cout << TMP_LOCAL_CM*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX))-0.5 << endl;
            //exit(0);
            CM[i] = Pre_CM[i] + 2.0*TMP_LOCAL_CM*( (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5);
            //cout << CM[i] << endl;
            //cout << gsl_ran_gaussian(r, SIGMAA) << endl;
        }
        //gsl_rng_free (r);
        // 2. Rotation:
        // a, Determine Inertia Tensor
        DetermineInertiaTensor(size, cen, xs, I); 
        //cout << I[0][0] << " " << I[0][1] << " " << I[0][2] << endl;
	// b, find eigenvector for Inertia Tensor
	//jacobi3(I, d, v, &n_rot);
        Meigen(v, d, I);
	// c, get unit vector
		
		
	if(random < 0.333) j = 0;
	else if(random > 0.667) j = 2;
	else j = 1;

        for (i = 0; i < 3; i++)
        {
            Rot[i] = v[i][j];
        }
        //cout << d[0] << " " << d[1] << " " << d[2] << endl;
        // d, Initialize angle
        tev = d[0]+d[1]+d[2];
        pe = d[j]/tev;
        Rot[3] = Pre_Rot[3] + 2.0*TMP_LOCAL_ROT*pe*exp(10.0/sqrt(tev))*( (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5);
        //cout << 2.0*TMP_LOCAL_ROT*pe*exp(10.0/sqrt(tev))*( (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5) << endl;
		
        // Translation
        for (i = 0; i < size; i++)
        {
            /*
            for (j = 0; j < 3; j++)
            {
                xs[i][j] = xs[i][j] + CM[j];
            }*/
            //cout << get_Xcoors(i) << endl;
            xs[i][0] = get_Xcoors(i) + CM[0]; //position of ligand inside tunnel.
            xs[i][1] = get_Ycoors(i) + CM[1];
            xs[i][2] = get_Zcoors(i) + CM[2];
            cen[0] += xs[i][0];
            cen[1] += xs[i][1];
            cen[2] += xs[i][2];
        }

        for (j = 0; j < 3; j++)
        {
            cen[j] /= size;
            
        }
		
        // Rotation
        c = cos(Rot[3]);
        s = sin(Rot[3]);

        uxx  = Rot[0]*Rot[0];
        uyy  = Rot[1]*Rot[1];
        uzz  = Rot[2]*Rot[2];
        uxy1 = Rot[0]*Rot[1]*(1-c);
        uxy2 = Rot[2]*s;
        uxz1 = Rot[0]*Rot[2]*(1-c);
        uxz2 = Rot[1]*s;
        uyz1 = Rot[1]*Rot[2]*(1-c);
        uyz2 = Rot[0]*s;

        R[0][0] = uxx + c*(1-uxx);
        R[0][1] = uxy1 - uxy2;
        R[0][2] = uxz1 + uxz2;
        R[1][0] = uxy1 + uxy2;
        R[1][1] = uyy + c*(1-uyy);
        R[1][2] = uyz1 - uyz2;
        R[2][0] = uxz1 - uxz2;
        R[2][1] = uyz1 + uyz2;
        R[2][2] = uzz + c*(1-uzz);
		

        for (i = 0; i < size; i++)
        {
            for (j = 0; j < 3; j++)
            {
                diff[j] = xs[i][j] - cen[j];
            }

            //y = MatrixMult(R, diff);
            MatrixMult(R, diff, y);

            for (j = 0; j < 3; j++)
            {
                xs[i][j] = cen[j] + y[j];
            }
        }
        //cout << "Till Here, everything's fine" << endl;
        //for (i = 0; i < LIGAND::m_liglist[0]->m_nGroups; i ++)
           //for (j = 0; j < 3; j++)
                //centroid_xs[i][j] = xs[LIGAND::m_liglist[0]->m_groups[i]][j];
        /*
	for(i = 0; i < size; i++)
	{
		xa[3*i] = xs[i][0];
		xa[3*i+1] = xs[i][1];
		xa[3*i+2] = xs[i][2];
	}
        */




	  //max_move = sqrt(   ( cent[0]-(par.C1[0]+par.C0[0])/2 )*( cent[0]-(par.C1[0]+par.C0[0])/2 ) + ( cent[1]-(par.C1[1]+par.C0[1])/2 )*( cent[1]-(par.C1[1]+par.C0[1])/2 ) + ( cent[2]-(par.C1[2]+par.C0[2])/2 )*( cent[2]-(par.C1[2]+par.C0[2])/2 )  );
           //max_move = sqrt(   (cen[0]- m_pointx)*(cen[0]- m_pointx) + (cen[1]- m_pointx)*(cen[1]- m_pointx) + (cen[2]- m_pointx)*(cen[2]- m_pointx) );
           max_move = sqrt(   (cen[0]- orig_cen[0])*(cen[0]- orig_cen[0]) + (cen[1]- orig_cen[1])*(cen[1]- orig_cen[1]) + (cen[2]- orig_cen[2])*(cen[2]- orig_cen[2]) );
           //cout << max_move << endl;
		
           // Assign new best configuration
           // dont accept conformations if rmsd_max > MC_max_rmsd
        
            if (max_move > MC_MAX_RMSD)
            {
                 cout << "max_move greater than maximum allowed movement" << endl;
                 rmsx++;
                 //reset translation, rotation and torsion
                 rmsdmax_counter++;
                 non_accepted_counter = 0;   //newly added
                 non_accepted_counterB = 0;  //newly added
                 for (i = 0; i < 3; i++)
                 {
                    Pre_CM[i] = 0.;
                 }
                 for (i = 0; i < 4; i++)
                 {
                    Pre_Rot[i] = 0.;
                 }
                 Pre_Rot[2]=1.0;
		 /*	
                 //reassign start coordinates and LP coordinates         
                 //for (i=0;i<lig->numatoms;i++)
                 for (i=0; i<size; i++)
                 {
                   for (j=0; j<3; j++)
                   {
                    save_xs[i][j]=orig_xs[i][j];
                   }
                 }
                 */
                 st--;
                 //cout << st << endl;
            }
            else
            {
                //MINST = st;

                //increase accepted solution counter
		//accept_solutions.push_back(make_pair(tmp_score,accept_counter_1));

		accept_counter_1++ ;
/*
		fprintf(lig_score, "MC Step %d: %f\n", st, tmp_score); //YY remove
//YY add
		fprintf(lig_trj, "@<TRIPOS>MOLECULE\n****\n");
		atomnum = lig->numatoms;
		bondnum = lig->numbonds;
		fprintf(lig_trj, "% 5d % 5d     0     0     0\n", atomnum, bondnum);
		fprintf(lig_trj, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n"); 

	
		for(i = 0; i < lig->numatoms; i++)
		{
			fprintf(lig_trj, "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 SUB       % 8.4f\n",
					lig->atm[i].atom_no, lig->atm[i].name,
					xa[3*i+0], xa[3*i+1], xa[3*i+2],
					lig->atm[i].sybyl_type, lig->atm[i].charge);
		}

		fprintf(lig_trj, "@<TRIPOS>BOND\n");
		while(!feof(fj))
		{
			fgets(line, lString, fj);
			if(strstr(line, "@<TRIPOS>BOND"))
				break;
		}
		for(i = 0; i < lig->numbonds; i++)
		{
			fgets(line, lString, fj);
			fprintf(lig_trj, line);
		}
		rewind(fj);
//YY end
*/
                non_accepted_counter = 0;   //this is tricky
                non_accepted_counterB = 0;

                //cout << "Till here, everything's fine" << endl;
                //store ligand cords
                //pxs = (float**)calloc(lig->numatoms,sizeof(float*));
                //pxs = new double[size][3];
                //cxs = new double[ligsize][3];
                for (i=0;i<size;i++)
                {
		        //pxs[i] = (float*)calloc(3,sizeof(float));
                        //pxs[i] = (float*)calloc(3,sizeof(float));
                    for (j=0;j<3;j++)
                    {
                        pxs[i][j] = xs[i][j];
                    //pxs.push_back(xs[i]);
                    }
                    
                }
                //cout << pxs[0][0] << " " << pxs[0][1] << " " << pxs[0][2] << endl; 

                for (i=0;i<ligsize;i++)
                {
                     int mark = LIGAND::m_liglist[0]->m_groups[i];
                     //cout << i << " " << mark << endl;
                     for (j=0;j<3;j++)
                        cxs[i][j] = xs[mark-1][j];
                        //cout << xs[mark-1][j] << endl;
                     
                }
                //cxs[0][0] = xs[3][0];
                //accept_coords.push_back(pxs);
                CLeaf * m_pLeaf = new CLeaf(0, pxs, 99);
                //cout << m_pLeaf->getLeafSize() << endl;			
                //exit(0);	
                VcV(m_pLeaf->m_translate, CM);
                      McM(m_pLeaf->m_rotate, R);
                            m_pLeaf->m_bv = new CRss(pxs,size);
                m_lLeaf.push_back(m_pLeaf);

                //cout << "m_orient" << endl;
                //cout << m_pLeaf->m_bv->getMOrient()[0][0] << " " << m_pLeaf->m_bv->getMOrient()[0][1] << " " << m_pLeaf->m_bv->getMOrient()[0][2] << endl;
                //cout << m_pLeaf->m_bv->getMOrient()[1][0] << " " << m_pLeaf->m_bv->getMOrient()[1][1] << " " << m_pLeaf->m_bv->getMOrient()[1][2] << endl;
                //cout << m_pLeaf->m_bv->getMOrient()[2][0] << " " << m_pLeaf->m_bv->getMOrient()[2][1] << " " << m_pLeaf->m_bv->getMOrient()[2][2] << endl;
                //cout << m_pLeaf->m_bv->getMPose(0) << " " << m_pLeaf->m_bv->getMPose(1) << " " << m_pLeaf->m_bv->getMPose(2) << endl;

                //cout << CM[0] << " " << CM[1] << " " << CM[2] << endl;
		//cout << m_pLeaf->m_translate[0] << " " << m_pLeaf->m_translate[1] << " " << m_pLeaf->m_translate[2] << endl;
                //cout << R[0][0] << " " << R[0][1] << " " << R[0][2] << endl;	
                //cout << m_pLeaf->m_rotate[0][0] << " " << m_pLeaf->m_rotate[0][1] << " " << m_pLeaf->m_rotate[0][2] << endl;
		//findNeighborLeaf(m_chain, pose, res_index);
                findNeighborLeaf(m_chain, pose);
                //cout << "res_index " << res_index[0] << endl;
                //reassign Pre values
                for (i = 0; i < 3; i++)
                {
                    Pre_CM[i] = 0.;
                }
                for (i = 0; i < 4; i++)
                {
                    Pre_Rot[i] = 0.;
                }
                Pre_Rot[2]=1.0;
            }


        /*
        //check how many solutions accepted and update local tor, rot and CM values
        //reassign local values every 100 steps depending on number accepted solutions
        if ((st%100)==0 && st!=0)
        {
            acc_factor = (float)accept_counter_1/(float)st;
            //betwenn 0.4 - 0.6 do nothing, else update
            if (acc_factor < 0.4 || acc_factor > 0.6)
            {
                //new factor = 2*acc if too small 0.5
                new_fact=2*acc_factor;
                if (new_fact<0.5)
                {
                    new_fact=0.5;
                }
                //update and check if values get too small, then set back too min values
                tmp_local_CM *= new_fact;
                tmp_local_Rot *= new_fact;

                if (tmp_local_CM< min_local_CM)
                {
                    tmp_local_CM = min_local_CM;
                }
                if (tmp_local_Rot< min_local_Rot)
                {
                    tmp_local_Rot = min_local_Rot;
                }
            }

        }

        if (((st+1)%1000)==0 && st!=0){
	    cout<<st+1<<" steps CM, accept rate:"<<(float)accept_counter_1/(float)(st+1)<<" Rot and Tors:"<<tmp_local_CM<<" "<<tmp_local_Rot<<" "<<tmp_local_Tors<<endl;
	}
        */
	//if(1 == level)
		//cout <<"CM:"<<tmp_local_CM<<endl;
       //cout << "Till Here, everything's fine" << endl;	
    }//end of monte carlo step loop
    //cout << "out of the loop" << endl;



	//calc final free energy
    //if 1: only save best solution
    // store min coordinates in MCsol
    /*
	for (i=0;i<lig->numatoms;i++)
	{
		for (j=0;j<3;j++)
		{
			lig->atm[i].optimized[solnum].x[j] = min_xs[i][j];
			lig->atm[i].sol[solnum].x[j] = min_xs[i][j];
		}
	}
	//printf("xmin: %f %f %f\n", min_xs[0][0], min_xs[0][1], min_xs[0][2]);
	//printf("minscore final %f\n", min_score);
	lig->optimized_score[solnum] = min_score;


 
    long double H_score = 0.00;
    long double G_score = 0.00;
    long double total_weight  = 0.00;
    long double weight;
    long double beta = 1.67841; 
    long double Weight = 1.00;
    long double piRMSD = 0.00, piwRMSD = 0.00;
    int Ndij, Ntotal = accept_counter_1;
    long double e = 2.718281828;
    long double omega_squ = 4.00;
    long double RMSD_THRESHOLD = 0.5;
    long double cluster_min = 99999.00;
    long double cluster_min_w = 99999.00;
    long double tmp_score1;
    bool flag_cluster;
  

    if(accept_solutions.size() != accept_cords.size() || accept_counter_1 != accept_solutions.size()){
	cout<<"Size doesn't match!"<<endl;
	exit(0);
    }
 

    float **RMSD = (float **)calloc(accept_counter_1, sizeof(float *));
    for(i=0; i<accept_counter_1; i++){
	RMSD[i] = (float *)calloc(accept_counter_1, sizeof(float));
    }

    
    float **pxs1;
    for(i=0; i<accept_counter_1; i++){
	pxs = accept_cords[i];
	for(j=i+1; j<accept_counter_1; j++){
	    pxs1 = accept_cords[j];
            rmsd_max = calc_rmsd(lig, pxs, pxs1);
	    RMSD[i][j] = rmsd_max;
	    RMSD[j][i] = rmsd_max;
        }
    }
    //cout <<"after RMSD"<<endl;
    

    int after_cluster = 0;
    vector<pair<long double, int> > af_cl_sol;
    sort(accept_solutions.begin(), accept_solutions.end());
    //cout <<"after sort, num of accepted sols: "<<accept_counter_1<<endl;

    //char tmpname2[256];
    
    //char *fname = strrchr(par.ListOfLigands[lig_nr], '/');

    //sprintf(tmpname2, "%s_Limoc_level%d.mol2", par.ligname[lig_nr], level);
    //fp = fopen(tmpname2, "w");

    if(accept_counter_1 != 0){
        af_cl_sol.push_back(accept_solutions[0]);
        //WriteCoordinates(num_flex_atoms, accept_prot_cords[0], lig, accept_cords[0], fp3);
        //WriteLigCoordinates(lig, accept_cords[0], fp4);
        //OutputLigand(lig, lig_nr, accept_cords[0], fp);
        //OutputLigand(lig, lig_nr, accept_cords[accept_solutions[0].second], fp);
        after_cluster++;
    }
    
    
    for(i=1; i<accept_counter_1; i++){

	
        flag_cluster = true;
        for(j=0; j<after_cluster; j++){
            if(RMSD[accept_solutions[i].second][af_cl_sol[j].second] < RMSD_THRESHOLD){
		flag_cluster = false;
                break;
            }
        }
        if(flag_cluster == true){
            after_cluster++;
	    af_cl_sol.push_back(accept_solutions[i]);
        }
	
        //WriteProtCoordinates(num_flex_atoms, accept_prot_cords[i], fp3);
        //WriteCoordinates(num_flex_atoms, accept_prot_cords[i], lig, accept_cords[i], fp3);
        //WriteLigCoordinates(lig, accept_cords[i], fp4);
        //OutputLigand(lig, lig_nr, accept_cords[i], fp);
        //OutputLigand(lig, lig_nr, accept_cords[accept_solutions[i].second], fp);
    }
    //cout <<"after clustering"<<endl;



    for(i=0; i<af_cl_sol.size(); i++){
	tmp_score = af_cl_sol[i].first;
	//cout<<"score: " <<tmp_score<<endl;
	weight = pow(2.71828, -beta*tmp_score);
	total_weight += weight;
	H_score += tmp_score * weight;
        //OutputLigand(lig, lig_nr, accept_cords[af_cl_sol[i].second], fp);
        //WriteCoordinates(num_flex_atoms, accept_prot_cords[af_cl_sol[i].second], lig, accept_cords[af_cl_sol[i].second], fp2);
        //WriteLigCoordinates(lig, accept_cords[af_cl_sol[i].second], fp5);
    }
    //cout <<"Z: "<<total_weight<<endl;


    //fclose(fp);


    if(total_weight != 0){
    	H_score /= total_weight;
	G_score = -0.5958 * log(total_weight);
    }
    else{
	H_score = min_score;
        G_score = min_score;
    }


    lig->G_score[solnum] = G_score;
    cout <<"number accepted before cluster:"<<accept_counter_1<<"\tafter:"<<after_cluster<<endl;
    cout<<"min_score: "<<min_score<<"\tH: "<<H_score<<"\tG: "<<G_score<<" TdeltaS: "<<H_score-G_score<<endl;

    fprintf(outfile, "%d\t%d\t%f\t%f\t%d\t%d\t%f\n",solnum, level, min_score, G_score, accept_counter_1, after_cluster, H_score-G_score);
    fflush( outfile );
    
    //cout <<"number accepted steps:"<<accept_counter_1<<"\tmin_score: "<<min_score<<endl;

    //YY add
    FILE *savefile;
    savefile = fopen("centroid", "w");
    fprintf(savefile, "#Piece of Reaction Coordinates\n# %f\n",A0);
    for (i=0;i<save_counter;i++)
    {
	fprintf(savefile, "%f\t%f\n", i+0.00001, save_cen[i]);
    }*/
	//fprintf(savefile, "#END\n");
	//fclose(savefile);
	//fclose(lig_trj); //YY add
	//fclose(lig_score); //YY remove

    //free memory
    /*
    for(i=0; i<accept_counter_1; i++){
        free(RMSD[i]);
    }
    free(RMSD);

    free(save_cen); //YY add


    for(i=0; i < accept_cords.size(); i++){
        pxs = accept_cords[i];
	for(j=0 ; j< lig->numatoms; j++){
	    free(pxs[j]);
        }
 	free(pxs);
    }
    
    for(i=0;i < lig->numatoms;i++)
    {
        free(xs[i]);
        free(save_xs[i]);
        free(orig_xs[i]);
        free(min_xs[i]);
    }
    free(xs);
    free(save_xs);
    free(orig_xs);
    free(min_xs);
    xs=NULL;
    save_xs=NULL;
    orig_xs=NULL;
    min_xs=NULL;
	

    if (lig->numTor > 0)
    {
        free(Pre_Tors);
        free(Save_Tors);
        free(Tors);
	free(Orig_Tors);
	free(Min_Tors);
    }
    free(xa);
    */
    

}

/*
//=================================================================================================
// Assign PLP type to protein atoms: 1, steric; 2, donor; 3, acceptor; 4, both.
//=================================================================================================
void AssignPLPtype(protein *prot)
{
	typedef char string4[5];
	int  NUMPLPDONORS = 38;
	static string4 plp_res_donor[] = {"ALA", "ARG", "ARG", "ARG", "ARG", "ASN", "ASN", "ASP", "ASH", "CYS", "CYX", "GLN", "GLN", "GLU", "GLH", "GLY",
                              "HIS", "HIS", "HIE", "HIE", "HID", "HID", "HIP", "HIP", "HIP", "ILE", "LEU", "LYS", "LYS", "MET",
                              "PHE", "PRO", "SER", "THR", "TRP", "TRP", "TYR", "VAL"};
	static string4 plp_atm_donor[] = {"N"  , "N"  , "NE" , "NH1", "NH2", "N",   "ND2", "N"  , "N"  , "N"  , "N"  , "N"  , "NE2", "N"  , "N"  , "N"  ,
                              "N"  , "ND1", "N"  , "NE2", "N"  , "ND1", "N"  , "ND1", "NE2", "N"  , "N"  , "N"  , "NZ" , "N"  ,
                              "N"  , "N"  , "N"  , "N"  , "N"  , "NE1", "N"  , "N"  };
	int NUMPLPACCEPTORS  = 37;
	static string4 plp_res_accep[] = {"ALA", "ARG", "ASN", "ASN", "ASP", "ASP", "ASP", "ASH", "ASH", "CYS", "CYX", "GLN", "GLN", "GLU", "GLU", "GLU", "GLH", "GLH", "GLY",
                              "HIS", "HIS", "HIE", "HIE", "HID", "HID", "HIP", "ILE", "LEU", "LYS", "MET",
                              "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
	static string4 plp_atm_accep[] = {"O"  , "O"  , "O"  , "OD1", "O"  , "OD1", "OD2", "O"  , "OD1", "O"  , "O"  , "O"  , "OE1", "O"  , "OE1", "OE2", "O"  , "OE1", "O"  ,
                              "O"  , "NE2", "O"  , "ND1", "O"  , "NE2", "O"  , "O"  , "O"  , "O"  , "O"  ,
                              "O"  , "O"  , "O"  , "O"  , "O"  , "O"  , "O"  };
	int NUMPLPDONACC  = 5;
	static string4 plp_res_donac[] = {"ASH", "GLH", "SER", "THR", "TYR"};
	static string4 plp_atm_donac[] = {"OD2", "OE2", "OG" , "OG1", "OH"};
	int		a, i;
	for(a = 0; a < prot->num_atoms; a++)
	{
		prot->atomlist[a].plptype = -1;
		if(prot->atomlist[a].name[0] != 'H')
		{
			for(i = 0; i < NUMPLPDONORS; i++)
			{				
				if(!strcmp(prot->atomlist[a].resname, plp_res_donor[i]) && !strcmp(prot->atomlist[a].name, plp_atm_donor[i]))
				{
					prot->atomlist[a].plptype = 2;
					goto L10;
				}
			}
		
			// check for acceptors
			for(i = 0; i < NUMPLPACCEPTORS; i++)
			{
				if(!strcmp(prot->atomlist[a].resname, plp_res_accep[i]) && !strcmp(prot->atomlist[a].name, plp_atm_accep[i]))
				{
					prot->atomlist[a].plptype = 3;
					goto L10;
				}
			}
		
			// check for donor-acceptors
			for(i = 0; i < NUMPLPDONACC; i++)
			{
				if(!strcmp(prot->atomlist[a].resname, plp_res_donac[i]) && !strcmp(prot->atomlist[a].name, plp_atm_donac[i]))
				{
					prot->atomlist[a].plptype = 4;
					goto L10;
				}
			}
		
			// assign steric groups
			if(!strncmp(prot->atomlist[a].name, "C", 1))
			{
				prot->atomlist[a].plptype = 1;
			}
			else if(!strncmp(prot->atomlist[a].name, "S", 1))
			{
				prot->atomlist[a].plptype = 1;
			}
			
			L10:;
		}
	}
}
*/

/*
template<typename T> 
bool Compare_Triple(const Triple<T> &i, const Triple<T> &j){
      return(i.data[2] < j.data[2]);
} 

template<typename V> 
bool Compare_Pair(const pair<vector<Triple<V> >, V> &i, const pair<vector<Triple<V> >, V> &j){
      return(i.second < j.second);
} 

template<typename T, typename V> 
bool Compare_Pair1(const pair<T , V> &i, const pair<T , V> &j){
      return(i.second < j.second);
}

template<typename V> 
inline unsigned Compare_Vector(const vector<Triple<V> > &i, const vector<Triple<V> > &j){
	unsigned num = 0;
	for(size_t iter=0; iter<i.size(); iter++){
		if(i[iter] != j[iter]){
			++num;
		}
	}
	return num;
}

template<typename T, typename V>
bool BinarySearch(PairVecF<Triple<T>, Triple<T>, V> &M, vector<Triple<T> > &Res_pos, unsigned start, unsigned end, const T &r1, const T &f1, const T &c1, const T &r2, const T &f2, const T &c2, float &strength, bool flag){
	T hash1 = r1 * 10000 + f1 * 100 + c1;
	T hash2 = r2 * 10000 + f2 * 100 + c2;
	T tempHashA, tempHashB;

        if(flag == true){
            T index;
            int i = 0,k;
            int j = Res_pos.size();
            if(r1 < r2){
                index = r1;
            }
            else{
                index = r2;
            }

            while(1){
                k = (i + j) / 2;
                if(Res_pos[k].data[0] == index){
                    start = Res_pos[k].data[1];
                    end = Res_pos[k].data[2];
                    break;
                }
                else if(Res_pos[k].data[0] > index){
                    j = k;
                }
                else{
                    i = k;
                }
            }
        }
        
	if(start == end){
		tempHashA = M[start].first.data[1] * 10000 + M[start].first.data[0] * 100 + M[start].first.data[2];
        tempHashB = M[start].second.data[1] * 10000 + M[start].second.data[0] * 100 + M[start].second.data[2];
		if(hash1 == tempHashA && hash2 == tempHashB || hash1 == tempHashB && hash2 == tempHashA){
                        strength += M.GetValue(start);
			return true;
                }
		else{
			return false;
                }
	}
	else{
	        unsigned mid = (start + end) / 2;
		T temp, tmpr1, tmpr2;
		if(hash1 > hash2){
			temp = hash1;
			hash1 = hash2;
			hash2 = temp;
 		}
		if(r1 < r2){
			tmpr1 = r1;
			tmpr2 = r2;
 		}
		else{
			tmpr1 = r2;
			tmpr2 = r1;		    
		}

		if(tmpr2 > M[mid].second.data[1]){
	            return BinarySearch(M, Res_pos, mid+1, end, r1, f1, c1, r2, f2, c2, strength,false);
		}
		else if(tmpr2 <  M[mid].second.data[1]){
		    return BinarySearch(M, Res_pos, start, mid, r1, f1, c1, r2, f2, c2, strength,false);
		}
		else{
		        //int lowerBound = mid, upperBound = mid;
			tempHashA = M[mid].first.data[1] * 10000 + M[mid].first.data[0] * 100 + M[mid].first.data[2];
                        if(tempHashA > hash1){
				return BinarySearch(M, Res_pos, start, mid, r1, f1, c1, r2, f2, c2, strength,false);
			}
                        else if(tempHashA < hash1){
				return BinarySearch(M, Res_pos, mid+1, end, r1, f1, c1, r2, f2, c2, strength,false);
                        }
                        else{
				tempHashB = M[mid].second.data[1] * 10000 + M[mid].second.data[0] * 100 + M[mid].second.data[2];
				if(tempHashB > hash2){
					return BinarySearch(M, Res_pos, start, mid, r1, f1, c1, r2, f2, c2, strength,false);
				}
                                else if(tempHashB < hash2){
					return BinarySearch(M, Res_pos, mid+1, end, r1, f1, c1, r2, f2, c2, strength,false);
				}
				else{
					strength += M.GetValue(mid);
					return true;
				}
			}
		}
	}
}

template<typename T, typename V>
vector<pair<vector<Triple<V> >, V> > CombineSolutions(PairVec<T,T> &P, PairVecF<Triple<T>, Triple<T>, V> &M, vector<Triple<T> > &Res_pos, vector<pair<vector<Triple<V> >, V> > &A, vector<pair<vector<Triple<V> >, V> > &B){
    vector<pair<vector<Triple<V> >, V> > tempV;
    size_t i, j, k, m, n, p, q, size_A, size_B, size_Ai, size_Bj, size_P, size_M;
    T res_id_A, traj_id_A, cluster_id_A;
    T res_id_B, traj_id_B, cluster_id_B;
    T tempA1, tempA2, tempA3, tempB1, tempB2, tempB3;
    bool flag_isPair, flag_hasMatch, b1, b2, b3, b4;
    float strength;
    vector<pair<vector<Triple<V> >, V> > Results;
    vector<Triple<V> > subV;
    size_A = A.size();
    size_B = B.size();
    size_M = M.size();
    if(size_A == 0 || size_B == 0){
	cout<<"No partical solution\n "<<endl;
	exit(0);
    }

    vector<pair<unsigned,unsigned> > index_pair;
    //cout<<"check solution: \n";
    for(i=0; i<A[0].first.size(); i++){
        for(j=0; j<B[0].first.size(); j++){
	    //cout<<(A[0].first)[i].data[0]<<"\t"<<(B[0].first)[j].data[0]<<endl;
            for(k=0; k<P.size();k++){
                if(P[k].first == (A[0].first)[i].data[0] && P[k].second == (B[0].first)[j].data[0] || P[k].first == (B[0].first)[j].data[0] && P[k].second == (A[0].first)[i].data[0]){
                        index_pair.push_back(make_pair(i,j));
                        break;
                }
            }
        }
    }
    for(i=0; i<size_A; i++){
        size_Ai = A[i].first.size();
	#ifdef _LIMOC_DEBUG
	if(i % 10 == 0) cout<<i<<"\n";
	#endif
        for(j=0; j<size_B; j++){
			//cout<<j<<"\n";
            size_Bj = B[j].first.size();
            //flag whether two partial solutions are competible with each other
            flag_hasMatch = true;
            strength = 0.0;
            for(k=0; k<index_pair.size(); k++){
		    m = index_pair[k].first;
		    n = index_pair[k].second;
		    //cout<<"m = "<<m<<"\t"<<"n = "<<n<<endl;
                    res_id_A = (T)(A[i].first)[m].data[0];
                    res_id_B = (T)(B[j].first)[n].data[0];


                    //size_P = P.size();
                    //flag_isPair = false;
                    //whether two residues are neighbor?
                    
                    //for(p=0; p<size_P; p++){
                        //if(P[p].first == res_id_A && P[p].second == res_id_B || P[p].second == res_id_A && P[p].first == res_id_B){
                            //flag_isPair = true;
                            //break;
                        //}
                    //}//for(p=0)
                    

                    //two residues are neighbors
                    //flag_hasMatch = false;

                    traj_id_A = (T)(A[i].first)[m].data[1];
                    cluster_id_A = (T)(A[i].first)[m].data[2];
                    traj_id_B = (T)(B[j].first)[n].data[1];
                    cluster_id_B = (T)(B[j].first)[n].data[2];
                    //are the traj id and cluster id of two solutions matching to each other 

	            flag_hasMatch = BinarySearch(M,Res_pos,0,size_M-1,res_id_A,traj_id_A,cluster_id_A,res_id_B,traj_id_B,cluster_id_B, strength,true);
                    if(flag_hasMatch == false){
                        //if(res_id_A >= 108 && res_id_A <= 115 && res_id_B >= 116 && res_id_B <= 122){
                            //cout<<res_id_A<<"\t"<<traj_id_A<<"\t"<<cluster_id_A<<"\t"<<res_id_B<<"\t"<<traj_id_B<<"\t"<<cluster_id_B<<endl;
                        //}
                        break;
                    }
            }//for(k=0)
            //find a solution
            if(flag_hasMatch == true){
		//cout<<"strength = "<<strength<<endl;
                subV.clear();
                for(m=0; m<size_Ai; m++){
                    subV.push_back((A[i].first)[m]);
                }
                for(n=0; n<size_Bj; n++){
                    subV.push_back((B[j].first)[n]);
                }
                //cout<<'\t'<<A[i].second<<'\t'<<B[j].second<<endl;
                Results.push_back(make_pair(subV,strength+A[i].second+B[j].second));
            }
        }//for(j=0)
    }//for(i=0)

    #ifdef _LIMOC_DEBUG
    cout<<"check fininshed"<<endl;
    #endif
    
    //A.clear();
    //for(i=0; i<Results.size(); i++){
        //A.push_back(Results[i]);
    //}
    
    #ifdef _LIMOC_DEBUG
    cout<<"combined solutions stored"<<endl;
    #endif

    //Results.clear();
    index_pair.clear();
    return Results;

}

template<typename T, typename V>
vector<pair<vector<Triple<V> >, V> > CombinationSearch(PairVec<T, vector<Triple<V> > > & PV, PairVec<T,T>& Pairs, PairVecF<Triple<T>, Triple<T>, V> &Matches, vector<Triple<T> > &Res_Match_pos, const vector<T>& MaxNum, unsigned augment, const size_t _begin, const size_t _end, const unsigned total, unsigned TruncateLevel, float mutation){
    size_t i, j;
    if(_begin > _end || _begin < 0 || _end < 0){ cout<<"Wrong range!\n"; exit(-1);}
    else if(_begin == _end){
        vector<pair<vector<Triple<V> >, V> > tempVV;
        vector<Triple<V> >subV;
        Triple<V> tempT;
        j = PV[_begin].second.size();
        if(j > (MaxNum[_begin] + augment))
            j = MaxNum[_begin] + augment;
        for(i=0; i<j; i++){
            subV.clear();
            tempT.data[0] = PV[_begin].first;
            tempT.data[1] = (PV[_begin].second)[i].data[0];
            tempT.data[2] = (PV[_begin].second)[i].data[1];
            subV.push_back(tempT);
            tempVV.push_back(make_pair(subV,0.0));
        }
        //cout<<PV[_begin].first<<"\t\n";
        //for(i=0; i<tempVV.size(); i++)
            //for(j=0; j<tempVV[i].size(); j++)
            //cout<<"\t\t"<<tempVV[i][j].data[0]<<"\t"<<tempVV[i][j].data[1]<<"\t"<<tempVV[i][j].data[2]<<endl;
        return tempVV;
    }
    else{
        vector<pair<vector<Triple<V> >, V> > A1;
        vector<pair<vector<Triple<V> >, V> > A2;
        unsigned mid = (_begin + _end) / 2;
        unsigned m = augment, n = augment;
        unsigned local = TruncateLevel;
        do{
            if(m > 20 || n > 20)
                break;
            if(A1.size() > 0 || A2.size() > 0){
                A1.clear();
                A2.clear();
            }
            while(A1.size() == 0){
                A1 = CombinationSearch(PV, Pairs, Matches, Res_Match_pos, MaxNum, m, _begin, mid, total, local, mutation);
                ++m;
                if(m > 20)
                    break;
            }
            while(A2.size() == 0){
                A2 = CombinationSearch(PV, Pairs, Matches, Res_Match_pos, MaxNum, n, mid+1, _end, total, local, mutation);
                ++n;
                if(n > 20)
                    break;
            }
            cout<<"-----------------------------------------"<<endl;
            cout<<"combining ["<<_begin<<","<<mid<<"]"<<"\t"<<A1.size()<<"\t"<<"\t"<<"["<<mid+1<<","<<_end<<"]"<<"\t"<<A2.size()<<"\t"<<"Truncate Level: "<<TruncateLevel<<"\tm = "<<m-1<<"\tn = "<<n-1<<endl;

            A1 = CombineSolutions<T, V>(Pairs, Matches, Res_Match_pos, A1, A2);
            if(m > 4){
                ++local;
            }
        }
        while(A1.size() == 0);

	unsigned size = A1.size();
        #ifdef _LIMOC_DEBUG
        cout<<"size of partial solution is: "<<size<<endl;
        #endif
        sort(A1.begin(), A1.end(), Compare_Pair<float>);
        #ifdef _LIMOC_DEBUG
	cout<<"sort complete"<<endl;
        #endif

	unsigned num_retain = TruncateLevel * 20 * (_end - _begin + 1);
        if(_end - _begin < total && size > num_retain){
        //if(size > num_retain){
	    vector<pair<vector<Triple<V> >, V> > Diverse;
	    unsigned size_D = 1;
	    bool flag;
	    //mutation rate
	    float threshold;
	    if(mutation > 0.0){
		threshold = A1[0].first.size() * mutation;
	    }
	    else{
		threshold = A1[0].first.size() * 0.2;
	    }
	    Diverse.push_back(A1[0]);
	    for(i=1; i<size; i++){
		flag = true;
		for(j=0; j<size_D; j++){
			if(Compare_Vector<V>(A1[i].first,Diverse[j].first) < threshold){
				flag = false;
				break;
			}
		}
		if(Diverse.size() < num_retain && flag == true){
			Diverse.push_back(A1[i]);
			++size_D;
		}
		else if(flag == true){
			break;
		}
	    }
            //A1.erase(A1.begin() + TruncateLevel * 20 * (_end - _begin + 1), A1.end());
	    A1.clear();
            #ifdef _LIMOC_DEBUG
	    cout<<"Diverse check complete"<<endl;
            #endif
	    return Diverse;
        }

        return A1;
    }
}
*/
/*
int LimocSearch(int argc, char *argv[], protein *prot, float ***cord, AtomR ** AR, ResCA ** resca, int ***clusters, int ***cluster_pos, int *NumTotal, GridScore **gridscore, int *NumGrid, Combination **comb, int *NumCombs, const int NC, const float mu, int ***sol_cl)
{
    unsigned i, j, k, m, n, counter;
    ifstream in_file;
    const size_t MAX = 256;
    char input_file[MAX];
    static size_t NumRes;

    char receiver[10];

    unsigned input_sign = 0;
    unsigned var1, var2;
    //intial number of clusters for each residue
    unsigned NumClusterRetain = 10;
    //Only consider residues within the Distance Threshold (angstrom)
    float DistanceThres = 8.0;
    string temp_str1, temp_str2, temp_str3,temp_str4;

    char buffer[MAX];
    bool flag;
    float fff, strength;
    vector<string> tempS;
    pair<unsigned, string> pus;
    Triple<float> _triple;
    Triple<unsigned> _tripleU;
    vector<Triple<float> > tempT;
    vector<Triple<unsigned> > tempU, tempU2;
    vector<unsigned> temp; 
    vector< pair< vector< Triple<float> >, float> > Combinations;
    cout<<"--------------------------------------------"<<endl;
    cout<<"Limoc Search"<<endl;

    for(j = 1; j < argc; j+=2)
    {
        
        //if(!strcmp(argv[j], "-n"))
        //{
            //strcpy(receiver, argv[j+1]);
            //NumRes = atoi(receiver);
            //cout<<"number of residues = "<<NumRes<<endl;
        //}
        
        if(!strcmp(argv[j], "-i"))
        {
            strcpy(input_file, argv[j+1]);
            input_sign = 1;

        }
        else if(!strcmp(argv[j], "-h"))
        {
            cout<<"Usage: Limoc_Scoring \n\t-n [Number of flexible residues]"<<endl;
            exit(0);
        }
        else if(!strcmp(argv[j], "-cl"))
        {
            strcpy(receiver, argv[j+1]);
            NumClusterRetain = atoi(receiver);
            //cout<<"number of starting clusters for each residue = "<<NumClusterRetain<<endl;
        }

    }
	if(NC > 0){
        	NumClusterRetain = NC;
                cout<<"number of starting clusters for each residue = "<<NumClusterRetain<<endl;
	}

    

        pair<unsigned, unsigned> varP;
        vector<pair<unsigned, unsigned> > tempP;
        PairVec<unsigned, unsigned> Patches_Pair;
        //PairVec<Triple<unsigned>, Triple<unsigned> > Matches_States;
        PairVecF<Triple<unsigned>, Triple<unsigned>, float> Matches_States;
        vector<Triple<unsigned> > Res_Match_pos;
        PairVec<unsigned, vector<unsigned>  > Patch_Partners;
        //changed because add a new attribute 
        PairVec<unsigned, vector<Triple<float> > > Patch_Candidate_States;

        vector<Triple<unsigned> >Threshold;

        //Read Patches and matching states bwtween Patches
        in_file.open("matching_Patches.txt");
        if(!in_file){
            cout<<"File open failed: cannot find matching_Patches.txt"<<endl;
            return EXIT_FAILURE;
        }
    
        vector<unsigned> BPRes;
        //BPRes[10] = 1;
        //BPRes[12] = 1;
        //BPRes[16] = 1;
        i = 0;
        flag = false;
        while(getline(in_file,temp_str1)){
            if(temp_str1[0] == '$'){
                if(temp_str1[1] == 'R'){
                    flag = true;
                    continue;
                }
                else{
                    break;
                }
            }
            if(true == flag){
                istringstream iss(temp_str1);
                iss >> temp_str2;
                iss >> temp_str3;
                iss >> fff;
		BPRes.push_back(0);
                if(fff < DistanceThres){
                    BPRes[i] = 1;
                }
                ++i;
            }
        }
        NumRes = i;
        cout<<"NumRes: "<<NumRes<<endl;
        *NumTotal = NumRes;
    
        unsigned NumBPRes = 0;
        for(i=0; i<NumRes; i++){
            if(BPRes[i] == 1){
                NumBPRes++;
            }
        }

        vector<unsigned> MaxNum(NumBPRes,NumClusterRetain);

        in_file.clear();
        in_file.seekg(0, ios::beg);
        while(in_file.getline(buffer,MAX)){
            if(buffer[0] == '$' && buffer[1] == 'P'){  //Patches
                copy(istream_iterator<unsigned>(in_file), istream_iterator<unsigned>(), back_inserter(temp));
                assert(temp.size()%2 == 0);
                vector<unsigned>::const_iterator it(temp.begin());
                while(it < temp.end() ){ 
                    if(BPRes[*it] == 1 && BPRes[*(it+1)] == 1){
                        Patches_Pair.PushBack(*it,*(it+1)); 
                    }
                    it += 2;
                }
                temp.clear();
            }
        }
        in_file.clear();
        in_file.seekg(0, ios::beg);
        while(in_file.getline(buffer,MAX)){
            if(buffer[0] == '$' && buffer[1] == 'M'){ //Matches
                temp.clear();
                Triple<unsigned> TTT;
                Triple<unsigned> TTT1;
	        while(in_file){
                    in_file >> TTT;
                    in_file >> TTT1;
	            in_file >> fff;
                    if(BPRes[TTT.data[1]] == 1 && BPRes[TTT1.data[1]] == 1){
		        Matches_States.PushBack(TTT,TTT1,fff);
                    }
	        }
            }
        }
        //in_file.close();

        
        //for(i=0; i<Patches_Pair.size(); i++){
            //cout<<"Pair "<<Patches_Pair[i].first<<"\t"<<Patches_Pair[i].second<<endl;
        //}
        
        //for(i=0; i<Matches_States.size(); i++){
            //if(Matches_States[i].second.data[1] == 95)
                //cout<<Matches_States[i].first<<"\t"<<Matches_States[i].second<<"\t"<<Matches_States.GetValue(i)<<endl;
        //}
        
        
        //cout<<"Num Match states: "<<Matches_States.size()<<endl;
        //cout<<"Num Pairs: "<<Patches_Pair.size()<<endl;
        //select all candidate states for each Patch: common group of states that appear in 
        //Matching states groups between a certain Patch and all other Patches 

        for(j=0; j<NumRes; j++){
            if(BPRes[j] == 1){
                temp.clear();
                Patch_Partners.PushBack(j, temp); 
                for(k=0; k<Patches_Pair.size(); k++){
                    if(Patches_Pair[k].first == j || Patches_Pair[k].second == j){
                        if(Patches_Pair[k].first == j)
                            i = Patches_Pair[k].second;
                        else
                            i = Patches_Pair[k].first;
                        //cout<<j<<" "<<i<<endl;
                        Merge(Patch_Partners, j, i);
                    }
                }
            }
        }

	#ifdef _LIMOC_DEBUG
        for(j=0; j<Patch_Partners.size(); j++){
            cout<<Patch_Partners[j].first<<endl;
            for(k=0; k<Patch_Partners[j].second.size(); k++){
                cout<<"\t\t"<<(Patch_Partners[j].second)[k]<<endl;
            } 
        }
	#endif

        unsigned index;
        for(j=0; j<NumBPRes; j++){
            //if(BPRes[j] == 1){
                index = Patch_Partners[j].first;
                //cout<<"index = "<<index<<endl;
                for(k=0; k<Patch_Partners[j].second.size(); k++){
                    tempT.clear();
                    temp.clear();
                    //partner residue number of residue j
                    i = (Patch_Partners[j].second)[k];
                    for(m=0; m<Matches_States.size(); m++){
                        var1 = Matches_States[m].first.data[1];
                        var2 = Matches_States[m].second.data[1]; 
                        if(var1 == index && var2 == i || var1 == i && var2 == index){
                            if(var1 == index)
                                varP = make_pair(Matches_States[m].first.data[0],Matches_States[m].first.data[2]);
                            else
                                varP = make_pair(Matches_States[m].second.data[0],Matches_States[m].second.data[2]);
                            //strength of the pair relationship between two clusters
                            strength = Matches_States.GetValue(m);
                            //if(index == 95){ cout<<"strength = "<<strength<<endl; }
                            bool flag = false;
                            for(n=0; n<tempT.size(); n++){
                                if(tempT[n].data[0] == varP.first && tempT[n].data[1] == varP.second){
                                    flag = true;
                                    //use the maximum strength between one cluster of residue j and all connect clusters of partner residue i 
                                    //if(strength > tempT[n].data[2])
                                    tempT[n].data[2] += strength;
                                    temp[n] += 1;
                                    break;
                                }
                            }
                            if(flag == false){
                                _triple.data[0] = (float)varP.first;
                                _triple.data[1] = (float)varP.second;
                                _triple.data[2] = strength;
                                tempT.push_back(_triple);
                                temp.push_back(1);
                            }
                        }
                    }//for(m=0
                    for(n=0; n<tempT.size(); n++){
                        
                        //if(index == 95){
                            //cout<<"temp: "<<temp[n]<<"\t"<<"tempT:"<<tempT[n]<<endl;
                        //}
                        
                        tempT[n].data[2] /= (float)temp[n];
                        
                        //if(index == 95){
                            //cout<<"tempT: "<<tempT[n].data[0]<<"\t"<<tempT[n].data[1]<<"\t"<<tempT[n].data[2]<<endl;
                        //}
                        
                    }
                    if(k == 0)
                        Patch_Candidate_States.PushBack(index, tempT);
                    else{
                        Intersect_T(Patch_Candidate_States, index, tempT);
                        
                        //if(index == 95){
                            //size_t tmp_size = Patch_Candidate_States.size();
                            //for(n=0; n<Patch_Candidate_States[tmp_size-1].second.size();n++)
                                //cout<<"instance: "<<i<<"\t"<<(Patch_Candidate_States[tmp_size-1].second)[n].data[0]<<"\t"<<(Patch_Candidate_States[tmp_size-1].second)[n].data[1]<<"\t "<<(Patch_Candidate_States[tmp_size-1].second)[n].data[2]<<endl;
                        //}
                        
                        
	            }	
                }//for(k=0
                if(Patch_Partners[j].second.size() == 0){
                    tempT.clear();
                    temp.clear();

        	    in_file.clear();
        	    in_file.seekg(0, ios::beg);
        	    while(in_file.getline(buffer,MAX)){
                        if(buffer[0] == '$' && buffer[1] == 'C'){ //Matches
                            temp.clear();
			    flag = false;
	                    while(getline(in_file,temp_str1)){
				if(temp_str1[0] == '$' && temp_str1[1] == 'P'){
					break;
				}
				else if(temp_str1[0] == 'r'){
                			istringstream iss1(temp_str1);
					unsigned cur_resid;
                			iss1 >> temp_str2;
                			iss1 >> cur_resid;
                			iss1 >> temp_str2;
                			iss1 >> var1;
					if(cur_resid == index){
						flag = true;
						//cout<<"var1"<<var1<<endl;
					}
					else if(flag == true){
						flag = false;
						break;
					}
				}
				else if(temp_str1[0] == 'c' && temp_str1[1] == 'l'){
					if(flag == true){
                				istringstream iss1(temp_str1);
                				iss1 >> temp_str2;
                				iss1 >> var2;
						//cout<<"var2"<<endl;
                                		_triple.data[0] = (float)var1;
                                		_triple.data[1] = (float)var2;
                                		_triple.data[2] = 0.0;
                                		tempT.push_back(_triple);		
					}
				}
	                    }
                        } 
                    }
*/
		    /*
                    //cout<<"index = "<<index<<endl;
                    for(m=0; m<Matches_States.size(); m++){
                        var1 = Matches_States[m].first.data[1];
                        var2 = Matches_States[m].second.data[1]; 
                        if(var1 == index || var2 == index){
                            //cout<<var1<<" "<<var2<<endl;
                            if(var1 == index)
                                varP = make_pair(Matches_States[m].first.data[0],Matches_States[m].first.data[2]);
                            else
                                varP = make_pair(Matches_States[m].second.data[0],Matches_States[m].second.data[2]);
                            //strength of the pair relationship between two clusters
                            strength = Matches_States.GetValue(m);
                            bool flag = false;
                            for(n=0; n<tempT.size(); n++){
                                if(tempT[n].data[0] == varP.first && tempT[n].data[1] == varP.second){
                                    flag = true;
                                    //use the maximum strength between one cluster of residue j and all connect clusters of partner residue i 
                                    //if(strength > tempT[n].data[2])
                                    tempT[n].data[2] += strength;
                                    temp[n] += 1;
                                    break;
                                }
                            }
                            if(flag == false){
                                _triple.data[0] = (float)varP.first;
                                _triple.data[1] = (float)varP.second;
                                _triple.data[2] = strength;
                                tempT.push_back(_triple);
                                temp.push_back(1);
                            }
                        }
                    }//for(m=0
                    for(n=0; n<tempT.size(); n++){
                        tempT[n].data[2] /= (float)temp[n];
                    }
		    */
/*
                    Patch_Candidate_States.PushBack(index, tempT);
                }//if(Patch_Partners[j].second.size() == 0)
            //}//if(BPRBPRes[j])
        }//for(j=0

        in_file.close();
*/        //exit(0);
        /*
        for(j=0; j<NumBPRes; j++){
            cout<<Patch_Candidate_States[j].first<<endl;
            for(k=0; k<Patch_Candidate_States[j].second.size(); k++){
                (Patch_Candidate_States[j].second)[k].data[2] /= (float)Patch_Partners[j].second.size();
                cout<<"\t\t"<<(Patch_Candidate_States[j].second)[k].data[0]<<"\t"<<(Patch_Candidate_States[j].second)[k].data[1]<<"\t"<<(Patch_Candidate_States[j].second)[k].data[2]<<endl;
            }
        }
        */
        //rank Patch_candidate states of each residue according to their average strength
        //A comparing function that compares the three value of two triples
/*        for(j=0; j<NumBPRes; j++){
            #ifdef _LIMOC_DEBUG
            cout<<Patch_Candidate_States[j].first<<endl;
	    #endif
            sort((Patch_Candidate_States[j].second).begin(), (Patch_Candidate_States[j].second).end(), Compare_Triple<float>);
	    if(Patch_Candidate_States[j].second.size() == 0){
		cout<<"couldn't identify any candidate cluster after matching check for residue: "<<j<<endl;
		exit(0);
	    }
            #ifdef _LIMOC_DEBUG
            for(k=0; k<Patch_Candidate_States[j].second.size(); k++){
		if(k<NumClusterRetain)
                cout<<"\t\t"<<(Patch_Candidate_States[j].second)[k].data[0]<<"\t"<<(Patch_Candidate_States[j].second)[k].data[1]<<"\t"<<(Patch_Candidate_States[j].second)[k].data[2]<<endl;

            }
	    #endif
        }
        //exit(0);

        //determine the number of clusters to retain for each residue
        unsigned cur_A =  Matches_States[0].first.data[1], cur_E = cur_A, cur_B = Matches_States[0].second.data[1], cur_C, cur_D;
        _tripleU.data[0] = cur_A;     _tripleU.data[1] = cur_B;     _tripleU.data[2] = 1;
        Threshold.push_back(_tripleU);
        _tripleU.data[0] = cur_E;     _tripleU.data[1] = 0;     _tripleU.data[2] = 0;
        Res_Match_pos.push_back(_tripleU);
        counter = 0;
        for(m=1; m<Matches_States.size(); m++){
            cur_C = Matches_States[m].first.data[1];
            cur_D = Matches_States[m].second.data[1];
            if(cur_C != cur_E){
                Res_Match_pos[Res_Match_pos.size()-1].data[2] = m - 1; 
                _tripleU.data[0] = cur_C;     _tripleU.data[1] = m;     _tripleU.data[2] = 0;
                Res_Match_pos.push_back(_tripleU); 
                cur_E = cur_C;
            }
            if(cur_C != cur_A || cur_D != cur_B){
                 _tripleU.data[0] = cur_C;
                 _tripleU.data[1] = cur_D;
                 cur_A = cur_C;
                 cur_B = cur_D;
                 _tripleU.data[2] = 1;
                 Threshold.push_back(_tripleU); 
                 counter += 1;
            }
            else{
                Threshold[counter].data[2] += 1;
            }
        }
        Res_Match_pos[Res_Match_pos.size()-1].data[2] = Matches_States.size() - 1;
*/        /*
        for(i=0; i<Res_Match_pos.size(); i++){
            cout<<"Res_pos: "<<Res_Match_pos[i].data[0]<<"\t"<<Res_Match_pos[i].data[1]<<"\t"<<Res_Match_pos[i].data[2]<<endl;
        }
        */
/*        for(j=0; j<Threshold.size(); j++){
            //cout<<Threshold[j].data[0]<<"\t"<<Threshold[j].data[1]<<"\t"<<Threshold[j].data[2]<<endl;
            var1 = Threshold[j].data[2];
            i = 0;
            if(var1 < 10){
                i = 4;
            }
            else if(var1 < 20){
                i = 3;
            }
            else if(var1 < 30){
                i = 2;
            }
            else if(var1 < 50){
                i = 1;
            }
            if(i > 0){
                m = 99999; n = 99999;
                for(k=0; k<Patch_Candidate_States.size(); k++){
                    var2 = Patch_Candidate_States[k].first;
                    if(var2 == Threshold[j].data[0]){
                        m = k;
                    }
                    else if(var2 == Threshold[j].data[1]){
                        n = k;
                    }
                    if(m != 99999 && n != 99999){
                        break;
                    }
                }
                //cout<<"m = "<<m<<" "<<"n = "<<n<<endl;
                if(m != 99999){
                    MaxNum[m] += i;
                }
                if(n != 99999){
                    MaxNum[n] += i;
                }
            }
        }
        //cout<<"test1"<<endl;
        cout<<"size of patch candidates: "<<Patch_Candidate_States.size()<<endl;
        cout<<"size of Max clusters: "<<MaxNum.size()<<endl;
*/
        /*
        for(i=0; i<MaxNum.size(); i++){
             cout<<Patch_Candidate_States[i].first<<'\t'<<MaxNum[i]<<endl;
        }
        */

/*        counter = Patch_Candidate_States.size() - 1;
        //counter = 36;

    if(input_sign == 0){
        //a divide-conquer method to find out all acceptable combinations of clusters
        Combinations = CombinationSearch<unsigned, float>(Patch_Candidate_States, Patches_Pair, Matches_States, Res_Match_pos, MaxNum, 0, 0, counter,counter,1,mu);
        //Combinations = CombinationSearch<unsigned, float>(Patch_Candidate_States, Patches_Pair, Matches_States, Res_Match_pos, MaxNum, 0, 115, 116, 2, 1,mu);

        ofstream out_file("solutions.txt");
        if(!out_file){
            cout<<"File open failed: cannot write to file: solutions.txt";
            return EXIT_FAILURE;
        }

        if(*NumCombs > Combinations.size()){
            *NumCombs = Combinations.size();
        }

        cout<<"combinations: \n";
        for(i=0; i<*NumCombs; i++){
            //cout<<"solution "<<i<<" "<<Combinations[i].second<<endl;
            out_file<<"$solution "<<i<<" "<<Combinations[i].second<<"\n";
            for(j=0; j<Combinations[i].first.size(); j++){
                //cout<<(Combinations[i].first)[j].data[1]<<"\t"<<(Combinations[i].first)[j].data[0]<<"\t"<<(Combinations[i].first)[j].data[2]<<"\t"<<endl;
                out_file<<(Combinations[i].first)[j].data[1]<<"\t"<<(Combinations[i].first)[j].data[0]<<"\t"<<(Combinations[i].first)[j].data[2]<<"\n";
            }
        }	

        out_file<<"END"<<"\n";
        out_file.close();
        strcpy(input_file, "solutions.txt");
	
	
    }//if(input_sign == 0)
    //else{
        //pre-generated solutions
        in_file.open(input_file);
        if(!in_file){
            cout<<"File open failed: cannot find "<<input_file<<endl;
            return EXIT_FAILURE;
        }
	    cout<<"Reading combinations from file: "<<input_file<<endl;
	Combinations.clear();
        tempT.clear();
        float prev_strength = 0.0;
        strength = 0.0;
        while(getline(in_file,temp_str1)){
            //cout<<temp_str1<<endl;
            if(temp_str1[0] == '$' && temp_str1[1] == 's'){ //solutions
                istringstream iss1(temp_str1);
                iss1 >> temp_str2;
                iss1 >> var1;
                iss1 >> fff;
                if(var1 == 0){
                    prev_strength = fff;
                }
                if(tempT.size() != 0){
                    Combinations.push_back(make_pair(tempT,prev_strength));
                    tempT.clear();
                    prev_strength = fff;
                }
                
            }
            else if(temp_str1[0] != ' ' && temp_str1[0] != 'E'){
                //cout<<temp_str1<<endl;
		istringstream sin(temp_str1);
		sin>>_triple.data[0];
		sin>>_triple.data[1];
		sin>>_triple.data[2];
                //cout<<_triple<<endl;
                tempT.push_back(_triple);
            }
            else if(temp_str1[0] == 'E'){
                Combinations.push_back(make_pair(tempT,prev_strength));
            }
        }
        in_file.close();
        
        if(*NumCombs != Combinations.size()){
            *NumCombs = Combinations.size();
        }

    //}//else
        
    //convert Combinations to C format and store in the file solutions.txt
    *comb = (Combination*) calloc (*NumCombs,sizeof(Combination));
    //cout<<"combinations: \n";
    for(i=0; i<*NumCombs; i++){
	(*comb)[i].strength = Combinations[i].second;
	(*comb)[i].sol = (ClusterT*) calloc (NumBPRes,sizeof(ClusterT));
	//cout<<"solution: "<<i<<"\t"<<comb[i].strength<<endl;
        for(j=0; j<Combinations[i].first.size(); j++){
	    assert(NumBPRes == Combinations[i].first.size());
	    ((*comb)[i].sol)[j].NumRes = (int)(Combinations[i].first)[j].data[1];
	    ((*comb)[i].sol)[j].NumFile = (int)(Combinations[i].first)[j].data[0];
	    ((*comb)[i].sol)[j].NumClt = (int)(Combinations[i].first)[j].data[2];
	    //cout<<(comb[i].sol)[j].NumRes<<"\t"<<(comb[i].sol)[j].NumFile <<"\t"<<(comb[i].sol)[j].NumClt<<endl;
        } 
    }
    cout<<"combinations intialization finished"<<endl;

    *sol_cl = (int **)calloc(NumRes, sizeof(int *));
    int **psol_cl = *(sol_cl);
    ClusterT *pCl;
    for(i=0; i<NumRes; i++){
	*(psol_cl+i) = (int *)calloc(MAXINDEX, sizeof(int));
    }
    for(i=0; i<*NumCombs; i++){
	for(j=0; j<NumRes; j++){
	    pCl = ((*comb)+i)->sol+j;
	    index = pCl->NumFile*100+pCl->NumClt;
	    psol_cl[j][index] = 1; 
	}
    }

	//cout<<"For test "<<"sign 21 8 10 = "<<psol_cl[21][810]<<endl;

   //read cluster files 
    string tail = ".pdb";
    string head = "clusters_for_residue_";
    PairVec<pair<unsigned, string>, vector<string> > Res_info;
    PairVec<unsigned, vector<Triple<unsigned> > > Res_cluster; 
    PairVec<unsigned, vector<Triple<unsigned> > > Res_cl_total; 
    PairVec<unsigned, vector<Triple<float> > > Prot_Cord;

    
    for(i=0; i<NumRes; i++){
        stringstream ss;
        ss << i;
        temp_str2 = ss.str();
        temp_str4 = head + temp_str2 + tail;
        //cout<<"i = "<<i<<"\t"<<"temp_str4 = "<<temp_str4<<endl;
        in_file.open(temp_str4.c_str());
        if(!in_file){
            cout<<"File open failed: cannot find "<<temp_str4<<endl;
            return EXIT_FAILURE;
        }
        flag = false;
        tempS.clear();
        tempU.clear();
	tempU2.clear();
        tempT.clear();
        while(getline(in_file,temp_str1)){
            //cout<<temp_str1<<" "<<temp_str1.length()<<endl;
            if(temp_str1[0] == 'M' && temp_str1[1] == 'O'){
                if(temp_str1.length()  < 22){
		    //cout<<"line = "<<temp_str1<<endl;
		    if(temp_str1[9] == ' '){
                        _tripleU.data[0] = atoi(temp_str1.substr(10,1).c_str());
		    }
		    else{
                        _tripleU.data[0] = atoi(temp_str1.substr(9,2).c_str());
		    }
		    if(temp_str1[15] == ' '){
                        _tripleU.data[1] = atoi(temp_str1.substr(16,1).c_str());
		    }
		    else{
			_tripleU.data[1] = atoi(temp_str1.substr(15,2).c_str());
		    }
                    //center of each cluster
                    _tripleU.data[2] = 999;
                }
                else{
		    if(temp_str1[9] == ' '){
                        _tripleU.data[0] = atoi(temp_str1.substr(10,1).c_str());
		    }
		    else{
                        _tripleU.data[0] = atoi(temp_str1.substr(9,2).c_str());
		    }
		    if(temp_str1[15] == ' '){
                        _tripleU.data[1] = atoi(temp_str1.substr(16,1).c_str());
		    }
		    else{
			_tripleU.data[1] = atoi(temp_str1.substr(15,2).c_str());
		    };
		    if(temp_str1[21] == ' '){
                        _tripleU.data[2] = atoi(temp_str1.substr(22,1).c_str());
		    }
		    else{
                        _tripleU.data[2] = atoi(temp_str1.substr(21,2).c_str());
		    }

                }

		index = _tripleU.data[0]*100+_tripleU.data[1];
		if(psol_cl[i][index] == 1){
                	tempU.push_back(_tripleU);
		}
		tempU2.push_back(_tripleU);
                //cout<<"ter1"<<endl;

            }
            else if(temp_str1[0] == 'E' && temp_str1[1] == 'N'){
                if(flag == false){
                    flag = true;
                    temp_str3 = temp_str4.substr(17,3);
                    pus = make_pair(i,temp_str3);
                    //cout<<"temp_str3 = "<<temp_str3<<endl;
                    Res_info.PushBack(pus, tempS);
                }
            }
            else if(temp_str1[0] == 'A' && temp_str1[1] == 'T'){
                temp_str4 = temp_str1;
                if(flag == false){
                    if(temp_str1[12] == ' ')
                        temp_str2 = temp_str1.substr(13,3);
                    else
                        temp_str2 = temp_str1.substr(12,4);
                    //temp_str3 = temp_str1.substr(17,3);
		    temp_str2.erase(temp_str2.find_last_not_of(" ") + 1);
                    tempS.push_back(temp_str2);
                }
                _triple.data[0] = atof(temp_str1.substr(31,7).c_str());
                _triple.data[1] = atof(temp_str1.substr(39,7).c_str());
                _triple.data[2] = atof(temp_str1.substr(47,7).c_str());
		if(psol_cl[i][index] == 1){
                	tempT.push_back(_triple);
		}
            }
        }//while(getline())
        Res_cluster.PushBack(i,tempU);
	Res_cl_total.PushBack(i,tempU2);
        Prot_Cord.PushBack(i,tempT);
        in_file.close();
    }

    //convert to C format 
    //initialize protein structure 
    int num_atoms = 0;
    for(i=0; i<NumRes; i++){
        num_atoms += Res_info[i].second.size();
    }
    cout<<"num_atoms: "<<num_atoms<<endl;
    prot->num_atoms = num_atoms;
    prot->atomlist = (Patom*) calloc (num_atoms,sizeof(Patom));

    Patom * pt;
    pt = prot->atomlist;
    *AR = (AtomR*) calloc (num_atoms,sizeof(AtomR));
    *resca = (ResCA*) calloc (*NumTotal,sizeof(ResCA));
    
    *cord = (float**) calloc (*NumTotal,sizeof(float *));
    float *pf;
    for(i=0; i<NumRes; i++){
        (*cord)[i] = (float*) calloc (3*Prot_Cord[i].second.size(),sizeof(float));
        pf = (*cord)[i];
        for(j=0; j<Prot_Cord[i].second.size(); j++){
            pf[3*j] = (Prot_Cord[i].second)[j].data[0];
            pf[3*j+1] = (Prot_Cord[i].second)[j].data[1];
            pf[3*j+2] = (Prot_Cord[i].second)[j].data[2];
        }
    }

    *clusters= (int**) calloc (*NumTotal,sizeof(int *));
    int *pi;
    for(i=0; i<NumRes; i++){
        (*clusters)[i] = (int *) calloc (3*Res_cluster[i].second.size() ,sizeof(int));
        pi = (*clusters)[i];
        for(j=0; j<Res_cluster[i].second.size(); j++){
            //if(i == 5){
                //cout<<"j = "<<j<<"\t"<<(Res_cluster[i].second)[j]<<endl;
            //}
            pi[3*j] = (Res_cluster[i].second)[j].data[0];
            pi[3*j+1] = (Res_cluster[i].second)[j].data[1];
            pi[3*j+2] = (Res_cluster[i].second)[j].data[2];
        }
    }

    
    k = 0;
    for(i=0; i<NumRes; i++){
        (*resca)[i].AtomNum = Res_info[i].second.size();
	(*resca)[i].atom_begin = k;
        for(j=0; j<Res_info[i].second.size(); j++){
            //cout<<"k = "<<k<<endl;
            pt[k].Num = k;
            strcpy(pt[k].name, (Res_info[i].second)[j].c_str());
	    strcpy(pt[k].resname, Res_info[i].first.second.c_str());
            (*AR)[k].NumRes = i;
            (*AR)[k].NumRespos = j;
	    ++k;
            //PLP type 
        }
	//(*resca)[i].atom_end = k - 1;
    }


    for(i=0; i<NumRes; i++){
        counter = 0;
        (*resca)[i].StateNum = Res_cluster[i].second.size();
        for(j=0; j<Res_cluster[i].second.size(); j++){
            if(999 == (Res_cluster[i].second)[j].data[2]){
                ++counter;
            }
        }
        (*resca)[i].ClNum = counter;
    }

    //Total number of clusters for residue not only the clusters appearing in the combinations
    vector<unsigned> Res_cl_num;
    for(i=0; i<NumRes; i++){
        counter = 0;
        for(j=0; j<Res_cl_total[i].second.size(); j++){
            if(999 == (Res_cl_total[i].second)[j].data[2]){
                ++counter;
            }
        }
        Res_cl_num.push_back(counter);
    }

    *cluster_pos = (int**) calloc (*NumTotal,sizeof(int *));

    for(i=0; i<NumRes; i++){
	(*cluster_pos)[i] = (int *) calloc (3*(*resca)[i].ClNum, sizeof(int));
        counter = 0;

        for(j=0; j<Res_cluster[i].second.size(); j++){
            if(999 == (Res_cluster[i].second)[j].data[2]){
		(*cluster_pos)[i][3*counter] = (Res_cluster[i].second)[j].data[0];
		(*cluster_pos)[i][3*counter+1] =(Res_cluster[i].second)[j].data[1];
		(*cluster_pos)[i][3*counter+2] = j; 
                ++counter;
            }

        }
    }
*/    /*
    for(i=0; i<NumRes; i++){
        cout<<"Num of clusters for res: "<<i<<"\t"<<Res_cl_num[i]<<endl;
    }
    */

    //Read Gridpoints and scores on gridpoints
/*    cout<<"Reading Grid scores..."<<endl;
    in_file.open("LimocScore.txt");
    if(!in_file){
        cout<<"File open failed: cannot find LimocScore.txt";
        return EXIT_FAILURE;
    }

    
    int cur_pos = -1;
    flag = false;
    unsigned index_cur_res, index_total_res;
    unsigned cur_size = 100000;
    vector<IRScore> vecIRS;
    IRScore tempIRS;
    IRScore *pIRS;
    *gridscore = (GridScore *)calloc(cur_size, sizeof(GridScore));
    GridScore *pg;
    while(getline(in_file,temp_str1)){
	istringstream iss(temp_str1);
        if(temp_str1[0] == '$'){
                //read a new grid point
                ++cur_pos;
                //need to reallocate memory
                if(cur_pos > cur_size){
 			cur_size += 100000;
			pg = (GridScore*) realloc (*gridscore, cur_size * sizeof(GridScore));
			if(pg == NULL){
				cout<<"bad memory allocation when resizing gridarray"<<endl;
				return EXIT_FAILURE;
			}
			else{
				*gridscore = pg;
			}
		}
		pg = *gridscore;
                iss >> temp_str2;
                iss >> cur_A;
		pg[cur_pos].NumGP = cur_A;
                iss >> fff;
		pg[cur_pos].GP_Cord[0] = fff;
                iss >> fff;
		pg[cur_pos].GP_Cord[1] = fff;
                iss >> fff;
		pg[cur_pos].GP_Cord[2] = fff;
                iss >> cur_B;
		pg[cur_pos].ERNum = cur_B;
		if(cur_B != 0){
			pg[cur_pos].ERlist = (int *)calloc(cur_B, sizeof(int));
			if(pg[cur_pos].ERlist == NULL){
				cout<<"bad memory allocation"<<endl;
				return EXIT_FAILURE;
			}	
		}
                else{
			pg[cur_pos].ERlist = NULL;
		}
                iss >> cur_D;
		pg[cur_pos].IRNum = cur_D;
		if(cur_D != 0){
			pg[cur_pos].IRlist = (int *)calloc(cur_D, sizeof(int));
			if(pg[cur_pos].IRlist == NULL){
				cout<<"bad memory allocation"<<endl;
				return EXIT_FAILURE;
			}	
		}
                else{
			pg[cur_pos].IRlist = NULL;
		}
        }//if(temp_str1[0] == '$')
	//explicit residue
	else if(temp_str1[0] == 'e'){
		iss >> temp_str2;
		for(i=0; iss >> cur_C; i++){
		    (pg[cur_pos].ERlist)[i] = cur_C;
		}
	}
        //average residue
	else if(temp_str1[0] == 'a'){
		vecIRS.clear();
		counter = 0;
		iss >> temp_str2;
		iss >> cur_E;
		(pg[cur_pos].IRlist)[counter] = cur_E;
		//num of states for current residue
		index_total_res = Res_cl_num[cur_E];
		index_cur_res = 1;

		tempIRS.cluster.NumRes = cur_E;
		iss >> cur_E;
		tempIRS.cluster.NumFile = cur_E;
		iss >> cur_E;
		tempIRS.cluster.NumClt = cur_E;

		iss >> fff;
		(tempIRS.Score)[0] = fff;
		iss >> fff;
		(tempIRS.Score)[1] = fff;
		iss >> fff;
		(tempIRS.Score)[2] = fff;
		iss >> fff;
		(tempIRS.Score)[3] = fff;

		index = tempIRS.cluster.NumFile*100 + tempIRS.cluster.NumClt;
		if(psol_cl[tempIRS.cluster.NumRes][index] == 1){
			vecIRS.push_back(tempIRS);
		}
	}
	//continue reading average residue
	else if(temp_str1 != ""){
		++index_cur_res;
		if(index_cur_res == index_total_res + 1){
			iss >> cur_E;
			++counter;
			index_total_res = Res_cl_num[cur_E];
			index_cur_res = 1;
			(pg[cur_pos].IRlist)[counter] = cur_E;
		}
		else{
			cur_E = (pg[cur_pos].IRlist)[counter];
		}
		tempIRS.cluster.NumRes = cur_E;
		iss >> cur_E;
		tempIRS.cluster.NumFile = cur_E;
		iss >> cur_E;
		tempIRS.cluster.NumClt = cur_E;
		iss >> fff;
		(tempIRS.Score)[0] = fff;
		iss >> fff;
		(tempIRS.Score)[1] = fff;
		iss >> fff;
		(tempIRS.Score)[2] = fff;
		iss >> fff;
		(tempIRS.Score)[3] = fff;

		index = tempIRS.cluster.NumFile*100 + tempIRS.cluster.NumClt;
		if(psol_cl[tempIRS.cluster.NumRes][index] == 1){
			vecIRS.push_back(tempIRS);
		}

		if(counter == pg[cur_pos].IRNum - 1  && index_cur_res == index_total_res){
			pg[cur_pos].irs = (IRScore *) calloc(vecIRS.size(), sizeof(IRScore));
			if(pg[cur_pos].irs == NULL){
				cout<<"bad memory allocation"<<endl;
				return EXIT_FAILURE;
			}
			pIRS = pg[cur_pos].irs;
			for(i=0; i<vecIRS.size(); i++){
				pIRS[i].cluster.NumRes = vecIRS[i].cluster.NumRes;
				pIRS[i].cluster.NumFile = vecIRS[i].cluster.NumFile;
				pIRS[i].cluster.NumClt = vecIRS[i].cluster.NumClt;
				(pIRS[i].Score)[0] = (vecIRS[i].Score)[0];
				(pIRS[i].Score)[1] = (vecIRS[i].Score)[1];
				(pIRS[i].Score)[2] = (vecIRS[i].Score)[2];
				(pIRS[i].Score)[3] = (vecIRS[i].Score)[3];	
			}
			pg[cur_pos].NumIRstates = vecIRS.size();
		}
	}
    }//while
    *NumGrid = cur_pos + 1;

    //Assign PLP type 
    cout<<"Assigning PLP type..."<<endl;
    AssignPLPtype(prot);


    return 1;
}
*/

/*=================================================================================================
GetClosesGridPoints
=================================================================================================*/
/*
int GetClosestGridPoints(const float * const xtmp, const float * const min, const float &delta, const int * const numGP, int *points)
{
	int	gx, gy, gz, closestGP, G_num;
	float	tx, ty, tz, alpha, beta, gama;
	float small = 0.00001;
	int numGP1 = numGP[2];
	int numGP2 = numGP[1] * numGP[2];

	alpha = modff( (xtmp[0] - min[0])/delta+small, &tx);
	beta  = modff( (xtmp[1] - min[1])/delta+small, &ty);
	gama = modff( (xtmp[2]  - min[2])/delta+small, &tz);
	gx = (int) (tx+small);
	gy = (int) (ty+small);
	gz = (int) (tz+small);
	//cout<<"grid: "<<gx<<"\t"<<gy<<"\t"<<gz<<endl;
	//cout<<numGP1<<"\t"<<numGP2<<endl;
	if (gx < 0 || gx >= numGP[0]-1 || gy < 0 || gy >= numGP[1]-1 || gz < 0 || gz >= numGP[2]-1)
	{
		G_num = -1;
                int test;
                for(test=0; test<8; test++){
                    points[test] = -1;
                }
	}
	else
	{
		if (alpha < small)
		{
			alpha = small;
		}
		if (beta < small)
		{
			beta = small;
		}
		if (gama < small)
		{
			gama = small;
		}
		closestGP = (int)(2*alpha) + 2* (int)(2*beta) + 4* (int)(2*gama);

		points[0] = gx*numGP2 + gy*numGP1 + gz;
		points[1] = (gx+1)*numGP2 + gy*numGP1 + gz;
		points[2] = gx*numGP2 + (gy+1)*numGP1 + gz;
		points[3] = (gx+1)*numGP2 + (gy+1)*numGP1 + gz;
		points[4] = gx*numGP2 + gy*numGP1 + (gz+1);
		points[5] = (gx+1)*numGP2 + gy*numGP1 + (gz+1);
		points[6] = gx*numGP2 + (gy+1)*numGP1 + (gz+1);
		points[7] = (gx+1)*numGP2 + (gy+1)*numGP1 + (gz+1);


		switch(closestGP)
		{
			case 0:
				G_num = points[0];
				break;

			case 1:
				G_num = points[1];
				break;

			case 2:
				G_num = points[2];
				break;

			case 3:
				G_num = points[3];
				break;

			case 4:
				G_num = points[4];
				break;

			case 5:
				G_num = points[5];
				break;

			case 6:
				G_num = points[6];
				break;

			case 7:
				G_num = points[7];
				break;

		}

	}
	return G_num;
}
*/
/*=================================================================================================
FindGridScore
=================================================================================================*/
/*
int FindGridScore(const GridScore * const gridscore, const int &NumGridScore, int &curGP){
	int start, end, mid;
	start = 0; end = NumGridScore - 1;
	int tmp;
        if(curGP == -1){
            return -1;
        }
	while(start < end - 1){
		mid = (start + end) / 2;
		tmp = gridscore[mid].NumGP;
		if( tmp == curGP){
			return mid;
		}
		else if(tmp > curGP){
			end = mid;
		}
		else{
			start = mid;
		}
	}
	tmp = gridscore[start].NumGP;
	if( tmp == curGP){
		return start;
	}
	tmp = gridscore[end].NumGP;
	if( tmp == curGP){
		return end;
	}
	return -1;

}
*/
/*=================================================================================================
FindLigScore
=================================================================================================*/
/*
int FindLigScore(const vector<pair<int, float> >& vec, int &hashkey){
	int start, end, mid;
	start = 0; end = vec.size() - 1;
	int tmp;
	while(start < end - 1){
		mid = (start + end) / 2;
		tmp = vec[mid].first;
		if( tmp == hashkey){
			return mid;
		}
		else if(tmp > hashkey){
			end = mid;
		}
		else{
			start = mid;
		}
	}
	tmp = vec[start].first;
	if( tmp == hashkey){
		return start;
	}
	tmp = vec[end].first;
	if( tmp == hashkey){
		return end;
	}
	return -1;

}
*/
/*=================================================================================================
FindClusterNum
=================================================================================================*/
/*
ClusterT* FindClusterNum(const int &resid, const Combination * const Comb, const int &NumRes){
	int start, end, mid;
	start = 0; end = NumRes - 1;
	int tmp;
	while(start < end - 1){
		mid = (start + end) / 2;
		tmp = Comb->sol[mid].NumRes;
		if( tmp == resid){
			return Comb->sol+mid;
		}
		else if(tmp > resid){
			end = mid;
		}
		else{
			start = mid;
		}
	}
	tmp = Comb->sol[start].NumRes;
	if( tmp == resid){
		return Comb->sol+start;
	}
	tmp = Comb->sol[end].NumRes;
	if( tmp == resid){
		return Comb->sol+end;
	}
	return NULL;
}
*/
/*=================================================================================================
FindIRNum
=================================================================================================*/
/*
int FindIRNum(const ClusterT *cl, const IRScore* irs, const int &Total){
	int start, end, mid;
	start = 0; end = Total - 1;
	int tmp;
	int key = 10000*cl->NumRes+ 100*cl->NumFile + cl->NumClt;
	while(start < end - 1){
		mid = (start + end) / 2;
		tmp = 10000*irs[mid].cluster.NumRes+ 100*irs[mid].cluster.NumFile + irs[mid].cluster.NumClt;
		if( tmp == key){
			return mid;
		}
		else if(tmp > key){
			end = mid;
		}
		else{
			start = mid;
		}
	}
	tmp = 10000*irs[start].cluster.NumRes+ 100*irs[start].cluster.NumFile + irs[start].cluster.NumClt;
	if( tmp == key){
		return start;
	}
	tmp = 10000*irs[end].cluster.NumRes+ 100*irs[end].cluster.NumFile + irs[end].cluster.NumClt;
	if( tmp == key){
		return end;
	}
	return -1;
}
*/
/*=================================================================================================
FindClusterCord
=================================================================================================*/
/*
inline float *FindClusterCord(int ** const clusters, int ** const cluster_pos, const ClusterT * const pCl, const int &NumRes, const ResCA * const resca, int &atom_index, float ** const prot_cord){
    int i;
    int resid = pCl->NumRes;
    int fileid = pCl->NumFile;
    int cltid = pCl->NumClt;
    int index;
    for(i=0; i<resca[resid].ClNum; i++){
	if(cluster_pos[resid][3*i] == fileid && cluster_pos[resid][3*i+1] == cltid){
		index = cluster_pos[resid][3*i+2];
		break;
	}
    }

    return prot_cord[resid] + 3 * resca[resid].AtomNum * index + 3 * atom_index;	
}
*/
/*=================================================================================================
FindClusterMemberCord
=================================================================================================*/
/*
inline float *FindClusterMemberCord(int ** const clusters, int ** const cluster_pos, const ClusterT * const pCl, const int &NumRes, const ResCA * const resca, int &atom_index, float ** const prot_cord, const int n){
    int i;
    int resid = pCl->NumRes;
    int fileid = pCl->NumFile;
    int cltid = pCl->NumClt;
    int index;
    for(i=0; i<resca[resid].ClNum; i++){
	if(cluster_pos[resid][3*i] == fileid && cluster_pos[resid][3*i+1] == cltid){
		index = cluster_pos[resid][3*i+2];
		break;
	}
    }
    //cout<<"resid: "<<resid<<"\t("<<fileid<<"\t"<<cltid<<"\t"<<")"<<"index = "<<index<<endl;

    return prot_cord[resid] + 3 * resca[resid].AtomNum * (index + n + 1) + 3 * atom_index;	
}
*/
/*=================================================================================================
FindCombCord
=================================================================================================*/
/*
int FindCombCord(int ** const clusters, const Combination * const Comb, const int NumRes, const ResCA * const resca, float ** const prot_cord, float * cur_p){
	int i,j, k;
    	int resid;
    	int fileid;
    	int cltid;
	ClusterT * pCl = Comb->sol;
	int index = 0;
	float *temp;
	for(i=0; i<NumRes; i++){
    		resid = pCl[i].NumRes;
    		fileid = pCl[i].NumFile;
    		cltid = pCl[i].NumClt;	
    		for(j=0; j<resca[resid].StateNum; j++){
        		if(clusters[resid][3*j] == fileid && clusters[resid][3*j+1] == cltid && clusters[resid][3*j+2] == 999){
           			break;	
			}
    		}
		for(k=0; k<resca[resid].AtomNum; k++){
			temp = prot_cord[resid] + 3 * resca[resid].AtomNum * j + 3 * k;
			cur_p[3 * index ] = temp[0];
			cur_p[3 * index + 1] = temp[1];
			cur_p[3 * index + 2] = temp[2];
			++index;
		}	
	}


	return 1;
}
*/
/*=================================================================================================
PLP_Limoc
=================================================================================================*/
/*
float PLP_Limoc(float *curLigcord, int *DA, const ClusterT* const pCl, const int &NumRes, const protein *const prot, float ** const prot_cord, int ** const clusters, int ** const cluster_pos, const ResCA * const resca, bool flag){
	int i,j,k;
	int resid, atomid;
	float *xtmp;
	float PLPscore = 0.0, dist, dist2;
	float energy;
	int PLPtype;
	char cc;

	int	st = -1;
	float	A[2], B[2], C[2], D[2], E[2], F[2];
	float	gA[2], gB[2], gC[2], gD[2];
	A[0] = 3.4;		A[1] = 2.3;
	B[0] = 3.6;		B[1] = 2.6;
	C[0] = 4.5;		C[1] = 3.1;
	D[0] = 5.5;		D[1] = 3.4;
	E[0] = -0.4;		E[1] = -2.0;
	F[0] = 20.0;		F[1] = 20.0;
	gA[0] = -5.882352941;    //-20.0/3.4	
	gA[1] = -8.695652174;    //-20.0/2.3
	gB[0] = -2.0;    	 //-0.4/0.2
	gB[1] = -6.666666667;	 //-2.0/0.3
	gC[0] = 0.0;		
	gC[1] = 0.0;
	gD[0] = 0.4;             //0.4/1.0	
	gD[1] = 6.6666666667;    //2.0/0.3


	//cout<<"Don: "<<DA[0]<<"\t"<<DA[1]<<endl;
	//calculate explicit scores

		//find the cluster num with current combination 
		//pCl = FindClusterNum(resid, Comb, NumRes);
		if(pCl == NULL){
			cout<<"wrong residue num: "<<resid<<endl;
			exit(0);
		}
		resid = pCl->NumRes;
		//cout<<"\tExplicit Res: "<<pCl->NumRes<<"\t"<<pCl->NumFile<<"\t"<<pCl->NumClt<<endl;
		//if(flag == true){
			//cout<<"\tLig atom cord: "<<curLigcord[0]<<"\t"<<curLigcord[1]<<"\t"<<curLigcord[2]<<endl;
		//}
		for(j=0; j<resca[resid].AtomNum; j++){
			atomid = resca[resid].atom_begin + j;
			cc = prot->atomlist[atomid].name[0];
			if(cc != 'H' && cc != '1' && cc != '2' && cc != '3'){
			    //find the coordinates for current residue/cluster/atom
			    xtmp = FindClusterCord(clusters, cluster_pos, pCl, NumRes, resca, j, prot_cord);

			    PLPtype = prot->atomlist[atomid].plptype;

			    //if(flag == true){
			    	//cout<<"\t\tatom: "<<j<<"\t"<<xtmp[0]<<"\t"<<xtmp[1]<<"\t"<<xtmp[2]<<"\tPLPtype: "<<PLPtype<<"\t";
			    //}
			    if((!DA[0] && !DA[1]) || (PLPtype == 1) || (DA[0] && !DA[1] && PLPtype == 2) || (DA[1] && !DA[0] && PLPtype == 3) ){
				st = 0;
			    }
			    else{
				st = 1;
			    }
			    DistSq(curLigcord, xtmp, &dist2);
			    dist = sqrt(dist2);
			    if(dist < D[st]){
					
				if(dist <= A[st]){
					energy = gA[st]*dist + F[st];
				}
				else if(dist <= B[st]){
					energy = gB[st]*dist - gB[st]*A[st];
				}
				else if(dist <= C[st]){
					energy = E[st];
				}
				else{
					energy = gD[st]*dist - gD[st]*D[st];
				}
					
				PLPscore += energy;
			    }
			    //if(flag == true){
			    	//cout<<"\t\tenergy: "<<energy<<endl;
			    //}	
			}		
		}//for(j=0; j<resca[resid].AtomNum; j++)
	

	//cout<<"\tScore: "<<PLPscore<<endl;
        return PLPscore;

}
*/
/*=================================================================================================
PLP_Limoc_lv2
=================================================================================================*/
/*
float PLP_Limoc_lv2(float *curLigcord, int *DA, const ClusterT* const pCl, const int &NumRes, const protein *const prot, float ** const prot_cord, int ** const clusters, int ** const cluster_pos, const ResCA * const resca, bool flag, const int n){
	int i,j,k;
	int resid, atomid;
	float *xtmp;
	float PLPscore = 0.0, dist, dist2;
	float energy;
	int PLPtype;
	char cc;

	int	st = -1;
	float	A[2], B[2], C[2], D[2], E[2], F[2];
	float	gA[2], gB[2], gC[2], gD[2];
	A[0] = 3.4;		A[1] = 2.3;
	B[0] = 3.6;		B[1] = 2.6;
	C[0] = 4.5;		C[1] = 3.1;
	D[0] = 5.5;		D[1] = 3.4;
	E[0] = -0.4;		E[1] = -2.0;
	F[0] = 20.0;		F[1] = 20.0;
	gA[0] = -5.882352941;    //-20.0/3.4	
	gA[1] = -8.695652174;    //-20.0/2.3
	gB[0] = -2.0;    	 //-0.4/0.2
	gB[1] = -6.666666667;	 //-2.0/0.3
	gC[0] = 0.0;		
	gC[1] = 0.0;
	gD[0] = 0.4;             //0.4/1.0	
	gD[1] = 6.6666666667;    //2.0/0.3


	//cout<<"Don: "<<DA[0]<<"\t"<<DA[1]<<endl;
	//calculate explicit scores

		//find the cluster num with current combination 
		//pCl = FindClusterNum(resid, Comb, NumRes);
		if(pCl == NULL){
			cout<<"wrong residue num: "<<resid<<endl;
			exit(0);
		}
		resid = pCl->NumRes;
		//cout<<"\tExplicit Res: "<<pCl->NumRes<<"\t"<<pCl->NumFile<<"\t"<<pCl->NumClt<<endl;
		//if(flag == true){
			//cout<<"\tLig atom cord: "<<curLigcord[0]<<"\t"<<curLigcord[1]<<"\t"<<curLigcord[2]<<endl;
		//}
		//cout<<"step1.501"<<"resid: "<<resid<<"\tNumatoms: "<<resca[resid].AtomNum<<endl;
    		//cout<<pCl->NumRes<<"\t"<<pCl->NumFile<<"\t"<<pCl->NumClt<<endl;
		for(j=0; j<resca[resid].AtomNum; j++){
			atomid = resca[resid].atom_begin + j;
			//cout<<"step1.502"<<endl;
			cc = prot->atomlist[atomid].name[0];
			//cout<<"step1.503"<<endl;
			if(cc != 'H' && cc != '1' && cc != '2' && cc != '3'){
			    //find the coordinates for current residue/cluster/atom
			    //cout<<"step1.51 "<<j<<"\t"<<n<<endl;
			    xtmp = FindClusterMemberCord(clusters, cluster_pos, pCl, NumRes, resca, j, prot_cord, n);
			    //cout<<"step1.52"<<endl;

			    PLPtype = prot->atomlist[atomid].plptype;

			    //if(flag == true){
			    	//cout<<"\t\tatom: "<<j<<"\t"<<xtmp[0]<<"\t"<<xtmp[1]<<"\t"<<xtmp[2]<<"\tPLPtype: "<<PLPtype<<"\t";
			    //}
			    if((!DA[0] && !DA[1]) || (PLPtype == 1) || (DA[0] && !DA[1] && PLPtype == 2) || (DA[1] && !DA[0] && PLPtype == 3) ){
				st = 0;
			    }
			    else{
				st = 1;
			    }
			    DistSq(curLigcord, xtmp, &dist2);
			    dist = sqrt(dist2);
			    if(dist < D[st]){
					
				if(dist <= A[st]){
					energy = gA[st]*dist + F[st];
				}
				else if(dist <= B[st]){
					energy = gB[st]*dist - gB[st]*A[st];
				}
				else if(dist <= C[st]){
					energy = E[st];
				}
				else{
					energy = gD[st]*dist - gD[st]*D[st];
				}
					
				PLPscore += energy;
			    }
			    //if(flag == true){
			    	//cout<<"\t\tenergy: "<<energy<<endl;
			    //}	
			}		
		}//for(j=0; j<resca[resid].AtomNum; j++)
	

	//cout<<"\tScore: "<<PLPscore<<endl;
        return PLPscore;

}
*/
/*=================================================================================================
trilinear
=================================================================================================*/
/*
inline float trilinear(const float * const avg, const float * const lig_cord, const float &delta, const float * const start){
	float a, b, c;
	a = lig_cord[0] - start[0];
	b = lig_cord[1] - start[1];
	c = lig_cord[2] - start[2];
	if(a<0 || b<0 || c<0 || a > delta || b > delta || c > delta){
		cout<<"Finding wrong grid points to interpolate"<<endl;
		cout<<"a = "<<a<<"\t"<<"b = "<<b<<"\t"<<"c = "<<c<<endl;
		cout<<"ligand cord: "<<lig_cord[0]<<"\t"<<lig_cord[1]<<"\t"<<lig_cord[2]<<endl;
		exit(0);
	}
	float g01, g23, g45, g67, g0123, g4567;
	g01 = avg[0] * (delta - a) + avg[1] * a;
	g23 = avg[2] * (delta - a) + avg[3] * a;
	g45 = avg[4] * (delta - a) + avg[5] * a;
	g67 = avg[6] * (delta - a) + avg[7] * a;
	g0123 = g01 * (delta - b) + g23 * b;;
	g4567 = g45 * (delta - b) + g67 * b;
	return (g0123 * (delta - c) + g4567 * c) / (delta*delta*delta);
}
*/
//YY add
/*===================================================================
getCentroid
=====================================================================*/
/*
float *getCentroid(float *xa,ligand *lig) 
{
    static float centroid[3];
    int j=0;

    float sum_x, sum_y, sum_z;
    sum_x = 0;
    sum_y = 0;
    sum_z = 0;

    int i;
    float	mass = 0.0;
    for (i = 0; i<lig->numatoms; i++ ) 
    {
    	if ( lig->atm[i].element == eB ) 
        {
            sum_x += xa[3*i] * mB; 
            sum_y += xa[3*i + 1] * mB;
            sum_z += xa[3*i + 2] * mB;
 	    j++;
	    mass += mB;
        }
    	if ( lig->atm[i].element == eC ) 
        {
            sum_x += xa[3*i] * mC;
            sum_y += xa[3*i + 1] * mC;
            sum_z += xa[3*i + 2] * mC;
 	    j++;
            mass += mC;
        }
    	if ( lig->atm[i].element == eN ) 
        {
            sum_x += xa[3*i] * mN;
            sum_y += xa[3*i + 1] * mN;
            sum_z += xa[3*i + 2] * mN;
 	    j++;
	    mass += mN;
        }
    	if ( lig->atm[i].element == eO ) 
        {
            sum_x += xa[3*i] * mO;
            sum_y += xa[3*i + 1] * mO;
            sum_z += xa[3*i + 2] * mO;
 	    j++;
		mass += mO;
        }
    	if ( lig->atm[i].element == eF ) 
        {
            sum_x += xa[3*i] * mF;
            sum_y += xa[3*i + 1] * mF;
            sum_z += xa[3*i + 2] * mF;
 	    j++;
		mass += mF;
        }
    	if ( lig->atm[i].element == eP ) 
        {
            sum_x += xa[3*i] * mP;
            sum_y += xa[3*i + 1] * mP;
            sum_z += xa[3*i + 2] * mP;
 	    j++;
		mass += mP;
        }
    	if ( lig->atm[i].element == eS ) 
        {
            sum_x += xa[3*i] * mS;
            sum_y += xa[3*i + 1] * mS;
            sum_z += xa[3*i + 2] * mS;
 	    j++;
		mass += mS;
        }
    	if ( lig->atm[i].element == eCL ) 
        {
            sum_x += xa[3*i] * mCL;
            sum_y += xa[3*i + 1] * mCL;
            sum_z += xa[3*i + 2] * mCL;
 	    j++;
		mass += mCL;
        }
    	if ( lig->atm[i].element == eBR ) 
        {
            sum_x += xa[3*i] * mBR;
            sum_y += xa[3*i + 1] * mBR;
            sum_z += xa[3*i + 2] * mBR;
 	    j++;
		mass += mBR;
        }
    	if ( lig->atm[i].element == eI ) 
        {
            sum_x += xa[3*i] * mI;
            sum_y += xa[3*i + 1] * mI;
            sum_z += xa[3*i + 2] * mI;
 	    j++;
		mass += mI;
        }
    	if ( lig->atm[i].element == eMe1 ) 
        {
            sum_x += xa[3*i] * mMe1;
            sum_y += xa[3*i + 1] * mMe1;
            sum_z += xa[3*i + 2] * mMe1;
 	    j++;
		mass += mMe1;
        }
    	if ( lig->atm[i].element == eMe2 ) 
        {
            sum_x += xa[3*i] * mMe2;
            sum_y += xa[3*i + 1] * mMe2;
            sum_z += xa[3*i + 2] * mMe2;
 	    j++;
		mass += mMe2;
        }
    	if ( lig->atm[i].element == eMe3 ) 
        {
            sum_x += xa[3*i] * mMe3;
            sum_y += xa[3*i + 1] * mMe3;
            sum_z += xa[3*i + 2] * mMe3;
 	    j++;
		mass += mMe3;
        }
    }   
    if (j != lig->numhatms)
    {
	printf("Number of heavy atoms does not match\n");
	exit(0);
    }
    centroid[0] = sum_x / mass; 
    centroid[1] = sum_y / mass; 
    centroid[2] = sum_z / mass; 
    
    return centroid;
}
*/

/*=================================================================================================
LimocScore
=================================================================================================*/
/*
float LimocScore(const float * const Gridmin, const float Griddelta, const int * const numGP, const protein * const prot_Limoc, float ** const prot_cord, int ** const clusters, int ** const cluster_pos, const ResCA * const resca, const GridScore * const gridscore, const int NumGridScore, const Combination * const Comb, const int NumCombs, const int NumRes, const ligand * const lig, const float * const lig_cord, int **sol_cl, int *CurComb){


	float *energy;
	energy = (float *)calloc(NumCombs, sizeof(float));
	if(energy == NULL){
		cout<<"Bad memory allocation"<<endl;
		exit(0);
	}
	
	int i, j, k, m, n, index;
	
	int closestGP;
	int resid;
	int GP8[8];
	//vector<unsigned> ExpRes;
	//vector<unsigned> AvgRes;
	unsigned ExpRes[100];
	unsigned AvgRes[100];
	unsigned Total_exp[200];
	unsigned Total_avg[200];
	//int exp_size, avg_size;	
	unsigned NumExp, NumAvg, NumTotal_exp, NumTotal_avg;
	//vector<unsigned> tempV;
	unsigned tmp;
	unsigned tempV[100];
	int GridIndex[8];
	float curLigcord[3];
	int indexGrid, indexGridClosest;
	int *pi;
    	IRScore *pIRS;
	ClusterT *pCl;
	bool flag;
	int count;
	int DA[2];
	float avgE[8];

	//vector<float > energy(NumCombs,0.0);
	bool flag_test = false;
	int sum;
	//need to be changed to a larger number if having more than 10 files 
	//float Score_exp[100][1000];
	float Score_avg[8][MAXINDEX];
	//float Score_tri_avg[100][1000];
	float **Score_exp, **Score_tri_avg;
	Score_exp = (float **)calloc(NumRes, sizeof(float*));
	if(Score_exp == NULL){
		cout<<"Bad memory allocation"<<endl;
		exit(0);
	}
	for(i=0; i<NumRes; i++){
		Score_exp[i] = (float *)calloc(MAXINDEX, sizeof(float));
	}	
	Score_tri_avg = (float **)calloc(NumRes, sizeof(float*));
	if(Score_tri_avg == NULL){
		cout<<"Bad memory allocation"<<endl;
		exit(0);
	}
	for(i=0; i<NumRes; i++){
		Score_tri_avg[i] = (float *)calloc(MAXINDEX, sizeof(float));
	}
	float tmp_energy;
	ClusterT tmp_cluster;
	float start[3], dis[3], dis2[3];
	float g01, g23, g45, g67, g0123, g4567;
	//float recip_cubic_delta = 1.00 / Griddelta*Griddelta*Griddelta;
	int tmp_res;
	float Punish = 0.0;

	

	NumTotal_exp = 0;
	NumTotal_avg = 0;
	//go through each ligand heavy atom 
	for(i = 0, count = 0; i < lig->numatoms; i++){
		if(lig->atm[i].element != eH){
			//get ligand atom coordinates 
			
			sum = 3*i;
			curLigcord[0] = lig_cord[sum];
			curLigcord[1] = lig_cord[sum+1];
			curLigcord[2] = lig_cord[sum+2];
			//get closest eight grid points for current ligand atom 
			//cout<<Gridmin[0]<<"\t"<<Griddelta<<"\t"<<numGP[0]<<endl;
			closestGP = GetClosestGridPoints(curLigcord, Gridmin, Griddelta, numGP, GP8);
			
			//cout<<"**************************************"<<endl;
			//cout<<"atom: "<<count<<"\nnearest grids: "<<endl;
			//for(k=0; k<8; k++){
				//cout<<GP8[k]<<"\t";
			//}
			//cout<<endl;
			//cout<<"closest grid: "<<closestGP<<endl;

			//find the explicit and average residues
			//ExpRes.clear();
			//AvgRes.clear();
			NumExp = 0;
			NumAvg = 0;
			
			indexGridClosest = FindGridScore(gridscore, NumGridScore, closestGP);
			
			
			if(indexGridClosest != -1){
				pi = gridscore[indexGridClosest].ERlist;
				for(m=0; m<gridscore[indexGridClosest].ERNum; m++){
					//ExpRes.push_back(pi[m]);
					ExpRes[NumExp] = pi[m];
					NumExp++;	
				}
			}
			else{
				for(n=0; n<8; n++){
					indexGridClosest = FindGridScore(gridscore, NumGridScore, GP8[n]);
					if(indexGridClosest != -1){
						closestGP = GP8[n];
						break;
					}
				}
				if(n != 8){
					pi = gridscore[indexGridClosest].ERlist;
					for(m=0; m<gridscore[indexGridClosest].ERNum; m++){
						//ExpRes.push_back(pi[m]);
						ExpRes[NumExp] = pi[m];
						NumExp++;
					}					
				}
				else{
					Punish += 20.0;
					continue;
				}
			}
			
			
			for(k=0; k<8; k++){
				if(GP8[k] == closestGP){
					indexGrid = indexGridClosest;
				}
				else{
					indexGrid = FindGridScore(gridscore, NumGridScore, GP8[k]);
				}
				GridIndex[k] = indexGrid;
				if(indexGrid != -1){
					pi = gridscore[indexGrid].IRlist;
					for(m=0; m<gridscore[indexGrid].IRNum; m++){
						flag = true;
						for(n=0; n<NumAvg; n++){
							if(AvgRes[n] == pi[m]){
								flag = false;
								break;
							}
						}
						if(flag == true){
							//AvgRes.push_back(pi[m]);
							AvgRes[NumAvg] = pi[m];
							NumAvg++;
						}
					}	
				}	
				
			}//for(k=0; k<8)
			//remove explicit residues from average residue list 
			//tempV.clear();
			
			tmp = 0;
			for(m=0; m<NumAvg; m++){
				flag = true;
				for(n=0; n<NumExp; n++){
					if(AvgRes[m] == ExpRes[n]){
						flag = false;
						break;
					}
				}
				if(flag == true){
					tempV[tmp] = AvgRes[m];
					tmp++;
				}
			}

			NumAvg = tmp;
			//LigGridRes.push_back(make_pair(ExpRes, tempV));

			if(NumTotal_exp != 0){
				for(k=0; k<NumExp; k++){
					flag = true;
					for(j=0; j<NumTotal_exp; j++){
						if(Total_exp[j] == ExpRes[k]){
							flag = false;
							break;
						}
					}
					if(flag == true){
						Total_exp[NumTotal_exp] = ExpRes[k];
						NumTotal_exp++;
					}
				}
			}
			else{
				for(k=0; k<NumExp; k++){
					Total_exp[NumTotal_exp] = ExpRes[k];
					NumTotal_exp++;
				}
			}

			if(NumTotal_avg != 0){
				for(k=0; k<NumAvg; k++){
					flag = true;
					for(j=0; j<NumTotal_avg; j++){
						if(Total_avg[j] == tempV[k]){
							flag = false;
							break;
						}
					}
					if(flag == true){
						Total_avg[NumTotal_avg] = tempV[k];
						NumTotal_avg++;
					}
				}
			}
			else{
				for(k=0; k<NumAvg; k++){
					Total_avg[NumTotal_avg] = tempV[k];
					NumTotal_avg++;
				}
			}

			
			DA[0] = lig->atm[i].Don;
			DA[1] = lig->atm[i].Acc;
			//cout<<"DA: "<<DA[0]<<"\t"<<DA[1]<<endl;

			//exp_size = NumExp;
			for(m=0; m<NumExp ; m++){
				resid = ExpRes[m];
				//cout<<"exp resid: "<<resid<<"\t";
				for(n=0; n<resca[resid].ClNum; n++){
					sum = 3*n;
					tmp_cluster.NumRes = resid;
					tmp_cluster.NumFile = cluster_pos[resid][sum];
					tmp_cluster.NumClt = cluster_pos[resid][sum+1];	
					index = tmp_cluster.NumFile*100 + tmp_cluster.NumClt;
					//cout<<"cluster: "<<tmp_cluster.NumFile<<"\t"<<tmp_cluster.NumClt<<endl;

					//if(sol_cl[resid][index] == 1){
					tmp_energy = PLP_Limoc(curLigcord, DA, &tmp_cluster, NumRes, prot_Limoc, 	prot_cord, clusters, cluster_pos, resca,true);
					//tmp_energy = 0.0;
					Score_exp[resid][index] += tmp_energy;
					//}

				}

				//pCl = Comb->sol+resid;
				//index = pCl->NumFile*100 + pCl->NumClt;
				//cout<<"Energy: "<<Score_exp[resid][index]<<"\t";
			}
			//cout<<endl;

			//calculate average energy

			for(k=0; k<8; k++){
				if(GridIndex[k] != -1){
					if(k == 0){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0];
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1];
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2];
						break;
					}
					else if(k == 1){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0] - Griddelta;
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1];
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2];
						break;
					}
					else if(k == 2){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0];
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1] - Griddelta;
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2];
						break;
					}
					else if(k == 3){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0] - Griddelta;
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1] - Griddelta;
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2];
						break;
					}
					else if(k == 4){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0];
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1];
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2] - Griddelta;
						break;
					}
					else if(k == 5){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0] - Griddelta;
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1];
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2] - Griddelta;
						break;
					}
					else if(k == 6){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0];
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1] - Griddelta;
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2] - Griddelta;
						break;
					}	
					else if(k == 7){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0] - Griddelta;
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1] - Griddelta;
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2] - Griddelta;
						break;
					}						
				}
			}//for(k=0; k<8; k++)
			dis[0] = (curLigcord[0] - start[0]) / Griddelta;
			dis[1] = (curLigcord[1] - start[1]) / Griddelta;
			dis[2] = (curLigcord[2] - start[2]) / Griddelta;
			if(dis[0] < -0.1 || dis[1] < -0.1 || dis[2] < -0.1 || dis[0] > 1.1 || dis[1] > 1.1 || dis[2] > 1.1){
				cout<<"wrong grid points"<<endl;
				cout<<"Ligand atom cord: "<<curLigcord[0]<<"\t"<<curLigcord[1]<<"\t"<<curLigcord[2]<<endl;
				for(k=0; k<8; k++){
					cout<<GP8[k]<<"\t";
				}	
				cout<<endl;		
				for(k=0; k<8; k++){
					cout<<GridIndex[k]<<"\t";
				}
				cout<<endl;
				cout<<"a = "<<dis[0]<<"\t"<<"b = "<<dis[1]<<"\t"<<"c = "<<dis[2]<<endl;
				exit(0);
			}	
			//dis2[0] = 1.00 - dis[0];
			//dis2[1] = 1.00 - dis[1];
			//dis2[2] = 1.00 - dis[2];
			 
			//cout<<dis[0]<<"\t"<<dis[1]<<"\t"<<dis[2]<<endl;

			//avg_size = NumAvg;
			for(m=0; m<NumAvg; m++){
				resid = tempV[m];
				//cout<<"avg resid: "<<resid<<"\t";
    				int total_states = resca[resid].ClNum;
				for(k=0; k<8; k++){
					flag = false;
					//if can identify a grid point in the LimocScore.txt
					if(GridIndex[k] != -1){
						sum = 0;
						for(n=0; n<(gridscore + GridIndex[k])->IRNum; n++){
							tmp_res = ((gridscore + GridIndex[k])->IRlist)[n];
							if(tmp_res == resid){
								flag = true;
								break;
							}
							sum += resca[tmp_res].ClNum;
						}
					}
					if(flag != true){
						//avgE[k] = 0.0;
						//for(n=0; n<1000; n++){
							//Score_avg[k][n] = 0.0;
						//}
						for(n=0; n<total_states; n++){
							index = cluster_pos[resid][3*n+0]*100 + cluster_pos[resid][3*n+1];
							if(sol_cl[resid][index] == 1){
								Score_avg[k][index] = 0.0;
							}
						}
					}
					else{
						//cout<<"\t"<<gridscore[GridIndex[k]].NumGP<<endl;
						//find the score on grid 
						pIRS = (gridscore + GridIndex[k])->irs + sum;
						for(n=0; n<total_states; n++){
							//support maxmium 10 Files and 100 clusters in each file
							//index = pIRS[n].cluster.NumFile*100+pIRS[n].cluster.NumClt;
							index = cluster_pos[resid][3*n+0]*100 + cluster_pos[resid][3*n+1];
							if(sol_cl[resid][index] == 1){
								if(DA[0] && !DA[1]){
									Score_avg[k][index] = (pIRS[n].Score)[0];
								}
								else if(!DA[0] && DA[1]){
									Score_avg[k][index] = (pIRS[n].Score)[1];
								}
								else if(DA[0] && DA[1]){
									Score_avg[k][index] = (pIRS[n].Score)[2];
								}
								else{
									Score_avg[k][index] = (pIRS[n].Score)[3];
								}
							}
							else{
								cout<<"test--test"<<endl;
							}
							//cout<<"\t\t"<<"index: "<<index<<"\t"<<Score_avg[k][index]<<endl;
						}
					}//else 

				}//for(k=0; k<8)


				for(n=0; n<total_states; n++){
					index = cluster_pos[resid][3*n+0]*100 + cluster_pos[resid][3*n+1];
					//if(sol_cl[resid][index] == 1){
					if(Score_avg[0][index] != 0 || Score_avg[1][index] != 0 || Score_avg[2][index] != 0 || Score_avg[3][index] != 0 || Score_avg[4][index] != 0 || Score_avg[5][index] != 0 || Score_avg[6][index] != 0 || Score_avg[7][index] != 0){
						
						g01 = Score_avg[0][index] + (Score_avg[1][index] - Score_avg[0][index]) * dis[0];
						g23 = Score_avg[2][index] + (Score_avg[3][index] - Score_avg[2][index]) * dis[0];
						g45 = Score_avg[4][index] + (Score_avg[5][index] - Score_avg[4][index]) * dis[0];
						g67 = Score_avg[6][index] + (Score_avg[7][index] - Score_avg[6][index]) * dis[0];
						g0123 = g01 + (g23 - g01) * dis[1];;
						g4567 = g45 + (g67 - g45) * dis[1];
						Score_tri_avg[resid][index] += (g0123 + (g4567 - g0123) * dis[2]);
					}
					//}
				}

				//pCl = Comb->sol+resid;
				//index = pCl->NumFile*100 + pCl->NumClt;
				//cout<<"Energy: "<<Score_tri_avg[resid][index]<<"\t";
			}//for(m=0; m<avg_size; m++)
			//cout<<endl;
			count++;
			
		}//if(lig->atm[i].element != eH)
	}// for(i = 0; i < lig->numatoms; i++)


	for(j=0; j<NumCombs; j++){
		for(m=0; m<NumTotal_exp; m++){
			resid = Total_exp[m];		
			pCl = (Comb+j)->sol+resid;
		//cout<<"\t"<<pCl->NumRes<<"\t"<<pCl->NumFile<<"\t"<<pCl->NumClt<<endl;
			index = pCl->NumFile*100 + pCl->NumClt;
		//explicit score
		//index = 100;

			//cout<<"resid: "<<resid<<"\tCluster: "<<pCl->NumFile<<"\t"<<pCl->NumClt;
			//cout<<"\tEnergy: "<<Score_exp[resid][index]<<endl;
			energy[j] += Score_exp[resid][index];
			//if(j == 745){
				//cout<<"resid: "<<resid<<"\tE"<<" cluster: ("<<pCl->NumFile<< "\t"<<pCl->NumClt<<")\t"<<Score_exp[resid][index]<<"\tsum: "<<energy[j]<<endl;
			//}
		}
	}//for(j=0; j<NumCombs; j++)

	for(j=0; j<NumCombs; j++){
		for(m=0; m<NumTotal_avg ; m++){
			resid = Total_avg[m];		
			pCl = (Comb+j)->sol+resid;
		//cout<<"\t"<<pCl->NumRes<<"\t"<<pCl->NumFile<<"\t"<<pCl->NumClt<<endl;
			index = pCl->NumFile*100 + pCl->NumClt;
		//explicit score
		//index = 100;

			//cout<<"resid: "<<resid<<"\tCluster: "<<pCl->NumFile<<"\t"<<pCl->NumClt;
			//cout<<"\tEnergy: "<<Score_tri_avg[resid][index]<<endl;
			energy[j] += Score_tri_avg[resid][index];
			//if(j == 745){
				//cout<<"resid: "<<resid<<"\tA"<<" cluster: ("<<pCl->NumFile<< "\t"<<pCl->NumClt<<")\t"<<Score_tri_avg[resid][index]<<"\tsum: "<<energy[j]<<endl;
			//}
		}
	}//for(j=0; j<NumCombs; j++)

      
	if(Score_exp != NULL){
		for(i=0; i<NumRes; i++){
			free(Score_exp[i]);
		}
		free(Score_exp);
	}
	if(Score_tri_avg != NULL){
		for(i=0; i<NumRes; i++){
			free(Score_tri_avg[i]);
		}
		free(Score_tri_avg);
	}

	
	float minE = energy[0];
	int minNum = 0;
		//cout<<"comb Energy: "<<0<<"\t"<<energy[0]<<endl;
	for(i=1; i<NumCombs; i++){
                //if(i%1000 == 0){
		    //cout<<"comb Energy: "<<i<<"\t"<<energy[i]<<endl;
                //}
		if(minE > energy[i]){
			minE = energy[i];
			minNum = i;
		}
	}


	*CurComb = minNum;
	//sort(energy.begin(), energy.end(), Compare_Pair1<int, float>);
	free(energy);
	return (minE+Punish)/SCALAR_FACTOR;
}
*/

/*=================================================================================================
LimocScore_lv2
=================================================================================================*/
/*
float LimocScore_lv2(const float * const Gridmin, const float Griddelta, const int * const numGP, const protein * const prot_Limoc, float ** const prot_cord, int ** const clusters, int ** const cluster_pos, const ResCA * const resca, const GridScore * const gridscore, const int NumGridScore, const Combination * const Comb, const int NumCombs, const int NumRes, const ligand * const lig, const float * const lig_cord, int **sol_cl, const int CombIndex){


	float energy = 0.0;
	//cout<<"Comb used is: "<<CombIndex<<endl;
	
	int i, j, k, m, n, index;
	
	int closestGP;
	int resid;
	int GP8[8];
	//vector<unsigned> ExpRes;
	//vector<unsigned> AvgRes;
	unsigned ExpRes[100];
	unsigned AvgRes[100];
	unsigned Total_exp[200];
	unsigned Total_avg[200];
	//int exp_size, avg_size;	
	unsigned NumExp, NumAvg, NumTotal_exp, NumTotal_avg;
	//vector<unsigned> tempV;
	unsigned tmp;
	unsigned tempV[100];
	int GridIndex[8];
	float curLigcord[3];
	int indexGrid, indexGridClosest;
	int *pi;
    	IRScore *pIRS;
	ClusterT *pCl;
	bool flag;
	int count;
	int DA[2];
	float avgE[8];

	//vector<float > energy(NumCombs,0.0);
	bool flag_test = false;
	int sum;
	//need to be changed to a larger number if having more than 10 files 
	//float Score_exp[100][1000];
	float Score_avg[8];
	//float Score_tri_avg[100][1000];
	float **Score_exp, *Score_tri_avg;
	int *MemNum;
	Score_exp = (float **)calloc(NumRes, sizeof(float*));
	if(Score_exp == NULL){
		cout<<"Bad memory allocation"<<endl;
		exit(0);
	}
	for(i=0; i<NumRes; i++){
		Score_exp[i] = (float *)calloc(100, sizeof(float));
	}	
	Score_tri_avg = (float *)calloc(NumRes, sizeof(float));
	if(Score_tri_avg == NULL){
		cout<<"Bad memory allocation"<<endl;
		exit(0);
	}
	MemNum = (int *)calloc(NumRes, sizeof(int));
	if(MemNum == NULL){
		cout<<"Bad memory allocation"<<endl;
		exit(0);
	}
	float tmp_energy;
	ClusterT tmp_cluster;
	float start[3], dis[3], dis2[3];
	float g01, g23, g45, g67, g0123, g4567;
	//float recip_cubic_delta = 1.00 / Griddelta*Griddelta*Griddelta;
	int tmp_res;
	float Punish = 0.0;
	int StartIndex, EndIndex;

	

	NumTotal_exp = 0;
	NumTotal_avg = 0;
	//go through each ligand heavy atom 
	for(i = 0, count = 0; i < lig->numatoms; i++){
		if(lig->atm[i].element != eH){
			//get ligand atom coordinates 
			
			sum = 3*i;
			curLigcord[0] = lig_cord[sum];
			curLigcord[1] = lig_cord[sum+1];
			curLigcord[2] = lig_cord[sum+2];
			//get closest eight grid points for current ligand atom 
			//cout<<Gridmin[0]<<"\t"<<Griddelta<<"\t"<<numGP[0]<<endl;
			closestGP = GetClosestGridPoints(curLigcord, Gridmin, Griddelta, numGP, GP8);
			
			//cout<<"**************************************"<<endl;
			//cout<<"atom: "<<count<<"\nnearest grids: "<<endl;
			//for(k=0; k<8; k++){
				//cout<<GP8[k]<<"\t";
			//}
			//cout<<endl;
			//cout<<"closest grid: "<<closestGP<<endl;

			//find the explicit and average residues
			//ExpRes.clear();
			//AvgRes.clear();
			NumExp = 0;
			NumAvg = 0;
			
			indexGridClosest = FindGridScore(gridscore, NumGridScore, closestGP);
			
			
			if(indexGridClosest != -1){
				pi = gridscore[indexGridClosest].ERlist;
				for(m=0; m<gridscore[indexGridClosest].ERNum; m++){
					//ExpRes.push_back(pi[m]);
					ExpRes[NumExp] = pi[m];
					NumExp++;	
				}
			}
			else{
				for(n=0; n<8; n++){
					indexGridClosest = FindGridScore(gridscore, NumGridScore, GP8[n]);
					if(indexGridClosest != -1){
						closestGP = GP8[n];
						break;
					}
				}
				if(n != 8){
					pi = gridscore[indexGridClosest].ERlist;
					for(m=0; m<gridscore[indexGridClosest].ERNum; m++){
						//ExpRes.push_back(pi[m]);
						ExpRes[NumExp] = pi[m];
						NumExp++;
					}					
				}
				else{
					//couldn't find even one near grid point
					Punish += 20.0;
					continue;
				}
			}
			
			
			for(k=0; k<8; k++){
				if(GP8[k] == closestGP){
					indexGrid = indexGridClosest;
				}
				else{
					indexGrid = FindGridScore(gridscore, NumGridScore, GP8[k]);
				}
				GridIndex[k] = indexGrid;
				if(indexGrid != -1){
					pi = gridscore[indexGrid].IRlist;
					for(m=0; m<gridscore[indexGrid].IRNum; m++){
						flag = true;
						for(n=0; n<NumAvg; n++){
							if(AvgRes[n] == pi[m]){
								flag = false;
								break;
							}
						}
						if(flag == true){
							//AvgRes.push_back(pi[m]);
							AvgRes[NumAvg] = pi[m];
							NumAvg++;
						}
					}	
				}	
				
			}//for(k=0; k<8)
			//remove explicit residues from average residue list 
			//tempV.clear();
			
			tmp = 0;
			for(m=0; m<NumAvg; m++){
				flag = true;
				for(n=0; n<NumExp; n++){
					if(AvgRes[m] == ExpRes[n]){
						flag = false;
						break;
					}
				}
				if(flag == true){
					tempV[tmp] = AvgRes[m];
					tmp++;
				}
			}

			NumAvg = tmp;
			//LigGridRes.push_back(make_pair(ExpRes, tempV));

			if(NumTotal_exp != 0){
				for(k=0; k<NumExp; k++){
					flag = true;
					for(j=0; j<NumTotal_exp; j++){
						if(Total_exp[j] == ExpRes[k]){
							flag = false;
							break;
						}
					}
					if(flag == true){
						Total_exp[NumTotal_exp] = ExpRes[k];
						NumTotal_exp++;
					}
				}
			}
			else{
				for(k=0; k<NumExp; k++){
					Total_exp[NumTotal_exp] = ExpRes[k];
					NumTotal_exp++;
				}
			}

			if(NumTotal_avg != 0){
				for(k=0; k<NumAvg; k++){
					flag = true;
					for(j=0; j<NumTotal_avg; j++){
						if(Total_avg[j] == tempV[k]){
							flag = false;
							break;
						}
					}
					if(flag == true){
						Total_avg[NumTotal_avg] = tempV[k];
						NumTotal_avg++;
					}
				}
			}
			else{
				for(k=0; k<NumAvg; k++){
					Total_avg[NumTotal_avg] = tempV[k];
					NumTotal_avg++;
				}
			}

			
			DA[0] = lig->atm[i].Don;
			DA[1] = lig->atm[i].Acc;
			//cout<<"DA: "<<DA[0]<<"\t"<<DA[1]<<endl;

			//cout<<"step1"<<endl;

			//exp_size = NumExp;
			for(m=0; m<NumExp ; m++){
				resid = ExpRes[m];
				//cout<<"exp resid: "<<resid<<"\t";
				pCl = (Comb+CombIndex)->sol+resid;
				tmp_cluster.NumRes = resid;
				tmp_cluster.NumFile = pCl->NumFile;
				tmp_cluster.NumClt = pCl->NumClt;
				flag = false;


				//cout<<"step1.49"<<"\t"<<resid<<"\t"<<tmp_cluster.NumFile<<"\t"<<tmp_cluster.NumClt<<endl;

				if(MemNum[resid] == 0){
					for(n=0; n<resca[resid].ClNum; n++){
						sum = 3*n;
						if(tmp_cluster.NumFile == cluster_pos[resid][sum] &&
						tmp_cluster.NumClt == cluster_pos[resid][sum+1]){
							flag = true;
							StartIndex = cluster_pos[resid][sum+2];
							if(n == resca[resid].ClNum - 1){
								int sts;
								for(sts=resca[resid].StateNum-1; ; sts--){
									if(clusters[resid][3*sts+2] == 999){
										EndIndex = StartIndex + resca[resid].StateNum - sts;
										break; 
									}
								}
							}
						}	
						else{
							if(flag == true){
								EndIndex = cluster_pos[resid][sum+2];
								break;
							}
						}
					//}
					}
					MemNum[resid] = EndIndex-StartIndex-1;
					//if(resid == 66){
						//cout<<"testtesttest: "<<MemNum[resid]<<endl;
					//}
				}

				//cout<<"step1.5"<<"\t"<<resid<<"\t"<<MemNum[resid]<<endl;

				for(n=0; n<MemNum[resid]; n++){
					//if(sol_cl[resid][index] == 1){
					//cout<<"step1.50"<<endl;
					tmp_energy = PLP_Limoc_lv2(curLigcord, DA, &tmp_cluster, NumRes, prot_Limoc, 	prot_cord, clusters, cluster_pos, resca, true, n);
					//cout<<"step1.6"<<endl;
					//tmp_energy = 0.0;
					Score_exp[resid][n] += tmp_energy;
					//cout<<"n = "<<n<<endl;
				}

				//pCl = Comb->sol+resid;
				//index = pCl->NumFile*100 + pCl->NumClt;
				//cout<<"Energy: "<<Score_exp[resid][index]<<"\t";
				//cout<<"step1.6"<<endl;
			}
			//cout<<endl;

			//cout<<"step2"<<endl;

			//calculate average energy

			for(k=0; k<8; k++){
				if(GridIndex[k] != -1){
					if(k == 0){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0];
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1];
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2];
						break;
					}
					else if(k == 1){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0] - Griddelta;
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1];
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2];
						break;
					}
					else if(k == 2){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0];
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1] - Griddelta;
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2];
						break;
					}
					else if(k == 3){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0] - Griddelta;
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1] - Griddelta;
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2];
						break;
					}
					else if(k == 4){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0];
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1];
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2] - Griddelta;
						break;
					}
					else if(k == 5){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0] - Griddelta;
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1];
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2] - Griddelta;
						break;
					}
					else if(k == 6){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0];
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1] - Griddelta;
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2] - Griddelta;
						break;
					}	
					else if(k == 7){
						start[0] = (gridscore+GridIndex[k])->GP_Cord[0] - Griddelta;
						start[1] = (gridscore+GridIndex[k])->GP_Cord[1] - Griddelta;
						start[2] = (gridscore+GridIndex[k])->GP_Cord[2] - Griddelta;
						break;
					}						
				}
			}//for(k=0; k<8; k++)
			dis[0] = (curLigcord[0] - start[0]) / Griddelta;
			dis[1] = (curLigcord[1] - start[1]) / Griddelta;
			dis[2] = (curLigcord[2] - start[2]) / Griddelta;
			if(dis[0] < -0.1 || dis[1] < -0.1 || dis[2] < -0.1 || dis[0] > 1.1 || dis[1] > 1.1 || dis[2] > 1.1){
				cout<<"wrong grid points"<<endl;
				cout<<"Ligand atom cord: "<<curLigcord[0]<<"\t"<<curLigcord[1]<<"\t"<<curLigcord[2]<<endl;
				for(k=0; k<8; k++){
					cout<<GP8[k]<<"\t";
				}	
				cout<<endl;		
				for(k=0; k<8; k++){
					cout<<GridIndex[k]<<"\t";
				}
				cout<<endl;
				cout<<"a = "<<dis[0]<<"\t"<<"b = "<<dis[1]<<"\t"<<"c = "<<dis[2]<<endl;
				exit(0);
			}	
			//dis2[0] = 1.00 - dis[0];
			//dis2[1] = 1.00 - dis[1];
			//dis2[2] = 1.00 - dis[2];
			 
			//cout<<dis[0]<<"\t"<<dis[1]<<"\t"<<dis[2]<<endl;

			//avg_size = NumAvg;
			for(m=0; m<NumAvg; m++){
				resid = tempV[m];
				//cout<<"avg resid: "<<resid<<"\t";
				pCl = (Comb+CombIndex)->sol+resid;
				tmp_cluster.NumRes = resid;
				tmp_cluster.NumFile = pCl->NumFile;
				tmp_cluster.NumClt = pCl->NumClt;
    				int total_states = resca[resid].ClNum;
				for(k=0; k<8; k++){
					flag = false;
					//if can identify a grid point in the LimocScore.txt
					if(GridIndex[k] != -1){
						sum = 0;
						for(n=0; n<(gridscore + GridIndex[k])->IRNum; n++){
							tmp_res = ((gridscore + GridIndex[k])->IRlist)[n];
							if(tmp_res == resid){
								flag = true;
								break;
							}
							sum += resca[tmp_res].ClNum;
						}
					}
					if(flag != true){
						//avgE[k] = 0.0;
						//for(n=0; n<1000; n++){
							//Score_avg[k][n] = 0.0;
						Score_avg[k] = 0.0;
					}
					else{
						//cout<<"\t"<<gridscore[GridIndex[k]].NumGP<<endl;
						//find the score on grid 
						pIRS = (gridscore + GridIndex[k])->irs + sum;
						for(n=0; n<total_states; n++){
							//support maxmium 10 Files and 100 clusters in each file
							//index = pIRS[n].cluster.NumFile*100+pIRS[n].cluster.NumClt;
							if(tmp_cluster.NumFile == cluster_pos[resid][3*n+0] && tmp_cluster.NumClt == cluster_pos[resid][3*n+1]){
								if(DA[0] && !DA[1]){
									Score_avg[k] = (pIRS[n].Score)[0];
								}
								else if(!DA[0] && DA[1]){
									Score_avg[k] = (pIRS[n].Score)[1];
								}
								else if(DA[0] && DA[1]){
									Score_avg[k] = (pIRS[n].Score)[2];
								}
								else{
									Score_avg[k] = (pIRS[n].Score)[3];
								}
								break;
							}
							//cout<<"\t\t"<<"index: "<<index<<"\t"<<Score_avg[k][index]<<endl;
						}
					}//else 

				}//for(k=0; k<8)


						
				g01 = Score_avg[0] + (Score_avg[1] - Score_avg[0]) * dis[0];
				g23 = Score_avg[2] + (Score_avg[3] - Score_avg[2]) * dis[0];
				g45 = Score_avg[4] + (Score_avg[5] - Score_avg[4]) * dis[0];
				g67 = Score_avg[6] + (Score_avg[7] - Score_avg[6]) * dis[0];
				g0123 = g01 + (g23 - g01) * dis[1];;
				g4567 = g45 + (g67 - g45) * dis[1];
				Score_tri_avg[resid] += (g0123 + (g4567 - g0123) * dis[2]);


				//pCl = Comb->sol+resid;
				//index = pCl->NumFile*100 + pCl->NumClt;
				//cout<<"Energy: "<<Score_tri_avg[resid][index]<<"\t";
			}//for(m=0; m<avg_size; m++)
			//cout<<endl;
			//cout<<"step2"<<endl;
			count++;
			
		}//if(lig->atm[i].element != eH)
	}// for(i = 0; i < lig->numatoms; i++)
	//cout<<"step3"<<endl;

	for(m=0; m<NumTotal_exp; m++){
		resid = Total_exp[m];
		tmp_energy = 999.0;  
		for(n=0; n< MemNum[resid]; n++){
			if(tmp_energy > Score_exp[resid][n]){     //!if all Scores are larger than +10.0, Just use mem 0
				tmp_energy = Score_exp[resid][n];
			}
		}
	
		
                if(tmp_energy > 10.0){  //sigle residue maximum energy: +10.0 
                    energy += 10.0;
                }
                else{
		    energy += tmp_energy;
                }
		//cout<<"resid: "<<resid<<"E\t"<<tmp_energy<<endl;
	}
        //cout <<"--------------------------------------------------------"<<endl;

	for(m=0; m<NumTotal_avg ; m++){
		resid = Total_avg[m];		

		//cout<<"resid: "<<resid<<"\tCluster: "<<pCl->NumFile<<"\t"<<pCl->NumClt;
		//cout<<"\tEnergy: "<<Score_tri_avg[resid][index]<<endl;
		energy += Score_tri_avg[resid];
		//cout<<"resid: "<<resid<<"A\t"<<Score_tri_avg[resid]<<endl;
	}

        
	if(Score_exp != NULL){
		for(i=0; i<NumRes; i++){
			free(Score_exp[i]);
		}
		free(Score_exp);
	}
	if(Score_tri_avg != NULL){
		free(Score_tri_avg);
	}
	if(MemNum != NULL){
		free(MemNum);
	}

	
	return (energy+Punish)/SCALAR_FACTOR;
}
*/

/*=================================================================================================
Initial_Limoc
=================================================================================================*/
/*
void Initial_Limoc(Limoc* limoc){
    limoc->prot_Limoc = (protein*) calloc (1,sizeof(protein));
    limoc->prot_cord = NULL;
    limoc->clusters = NULL;
    limoc->cluster_pos = NULL;
    limoc->NumRes = (int*) calloc (1,sizeof(int));
    limoc->AtomRes = NULL;
    limoc->resca = NULL;
    limoc->gridscore = NULL;
    limoc->NumGridpoints = (int*) calloc (1,sizeof(int));
    limoc->NumCombs = (int*) calloc (1,sizeof(int));
    *(limoc->NumCombs) = 10000;
    limoc->Comb = NULL;
    limoc->sol_cl = NULL;
}
*/
/*=================================================================================================
Free_Limoc
=================================================================================================*/
/*
void Free_Limoc(Limoc* limoc){
    int size = 0, ncomb = 0, ngrid = 0;
    if(limoc->NumRes != NULL){
    	size = *(limoc->NumRes);
    	free(limoc->NumRes);
    }
    if(limoc->NumCombs != NULL){
        ncomb = *(limoc->NumCombs);
	free(limoc->NumCombs);
    }
    if(limoc->NumGridpoints != NULL){
        ngrid = *(limoc->NumGridpoints);
	free(limoc->NumGridpoints);
    }
    int i;
    //free everything

    if(limoc->prot_cord != NULL){
    	for(i=0; i<size; i++){
        	free(limoc->prot_cord[i]);
    	}
    	free(limoc->prot_cord);
    }

    if(limoc->clusters != NULL){
    	for(i=0; i<size; i++){
        	free(limoc->clusters[i]);
    	}
    	free(limoc->clusters);
    }

    if(limoc->cluster_pos != NULL){
    	for(i=0; i<size; i++){
        	free(limoc->cluster_pos[i]);
    	}
    	free(limoc->cluster_pos);
    }

    if(limoc->AtomRes != NULL){
        free(limoc->AtomRes);
    }
    if(limoc->resca != NULL){
        free(limoc->resca);
    }

    if(limoc->prot_Limoc != NULL){
        if(limoc->prot_Limoc->atomlist != NULL){
            free(limoc->prot_Limoc->atomlist);
        }
        free(limoc->prot_Limoc);
    }
    if(limoc->Comb != NULL){
        for(i=0; i<ncomb; i++){
	    //cout<<"i = "<<i<<endl;
            free(limoc->Comb[i].sol);	
        }
        free(limoc->Comb);
    }

    if(limoc->sol_cl != NULL){
	for(i=0; i<size; i++){
		free(limoc->sol_cl[i]);
	}
	free(limoc->sol_cl);
    }

    if(limoc->gridscore != NULL){
        for(i=0; i<ngrid; i++){
	    //cout<<"i = "<<i<<endl;
            free(limoc->gridscore[i].ERlist);
            free(limoc->gridscore[i].IRlist);
            free(limoc->gridscore[i].irs);	
        }
        free(limoc->gridscore);
    }
}
*/
/*=================================================================================================
Limoc
=================================================================================================*/
/*
void LimocRefine(int argc, char *argv[], ligand *lig, float *min, float delta, int *numGP, int lig_nr){
    int i, j, index;
    ClusterT* pCl;
    cout<<"Limoc Start"<<endl;
    
    if(lig == NULL){
        cout<<"Please input valid ligand structure"<<endl;
	exit(0);
    }
    
    Limoc limoc;
    Initial_Limoc(&limoc);
    cout<<"After Initialization"<<endl;

 
    limoc.min[0] = min[0];
    limoc.min[1] = min[1];
    limoc.min[2] = min[2];
    limoc.delta = delta;
    limoc.numGP[0] = numGP[0];
    limoc.numGP[1] = numGP[1];
    limoc.numGP[2] = numGP[2];


    
    //limoc.min[0] = -14.222;
    //limoc.min[1] = -1.757;
    //limoc.min[2] = 26.74;
    //limoc.delta =0.4;
    //limoc.numGP[0] = 101;
    //limoc.numGP[1] = 96;
    //limoc.numGP[2] = 106;
    

    cout <<"Grid information:"<<endl;
    cout << limoc.min[0] << " " << limoc.min[1] << " "<<limoc.min[2] << " delta: " << delta << " numGP: "<<limoc.numGP[0]<<" "<<limoc.numGP[1]<<" "<<limoc.numGP[2]<<endl;

    i = LimocSearch(argc, argv, limoc.prot_Limoc, &limoc.prot_cord, &limoc.AtomRes, &limoc.resca, &limoc.clusters, &limoc.cluster_pos, limoc.NumRes, &limoc.gridscore, limoc.NumGridpoints, &limoc.Comb, limoc.NumCombs, 10, 0.2, &limoc.sol_cl);

    //size_t mem_used = memory_used();
    //cout << "Memory usage after LimocSearch: " << mem_used <<endl;

    cout<<"after LimocSearch"<<endl;

    //FindGridScore();
    
    //int GP[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    //float curlig[3] = {44.2823, 60.9727, 70.0562};


    //GetClosestGridPoints(curlig, limoc.min, limoc.delta, limoc.numGP, GP);

    //for(i = 0; i < 8; i++){
        //cout<<"Grid point: "<<GP[i]<<"\t";
        //j = FindGridScore(limoc.gridscore, *limoc.NumGridpoints, GP[i]);
        //cout<<"index: "<<j<<endl;
    //}
    
    

    //unsigned **sol_cl;

    MonteCarloLimoc(lig, limoc.prot_Limoc, limoc.min, limoc.delta, limoc.numGP, limoc.prot_cord, limoc.clusters, limoc.cluster_pos, limoc.resca, limoc.gridscore, *(limoc.NumGridpoints), limoc.Comb, *(limoc.NumCombs), *(limoc.NumRes), limoc.sol_cl, lig_nr);
    cout<<"Limoc Refinement finished"<<endl;

    //mem_used = memory_used();
    //cout << "Memory usage after MonteCarloLimoc: " << mem_used <<endl;
    
    Free_Limoc(&limoc);
    cout <<"Limoc freed"<<endl;

}
*/

//=================================================================================================
// find_clash
//=================================================================================================
//bool MonteCarloLigand::findClash()
//void MonteCarloLigand::findNeighborLeaf(CChain & m_chain, int pose, vector<int > res_index)
void MonteCarloLigand::findNeighborLeaf(CChain & m_chain, int pose)
{
  //vector<CLeaf * > res_list;
  //vector<int *> res_list;
  // Find C-alpha or center within 8Angs of any of the ligand center
  res_index.clear();
  int l = 0;
  double (*tmp)[3];
  for (int i = 0; i < m_chain.getLength(); i++)
    { 
      int dist = 0.0;
      
      for (int j = 0; j < LIGAND::m_liglist[0]->m_nGroups; j++)
       {
          //dist = sqrt (pow((m_chain.getLink(i)->getCenter()[0] - LIGAND::m_liglist[0]->m_lconformer[pose]->m_centers[j][0]),2) + pow((m_chain.getLink(i)->getCenter()[1] - LIGAND::m_liglist[0]->m_lconformer[pose]->m_centers[j][1]),2) + pow((m_chain.getLink(i)->getCenter()[2] - LIGAND::m_liglist[0]->m_lconformer[pose]->m_centers[j][2]),2) );
          dist = sqrt (pow((m_chain.getLink(i)->getCenter()[0] - cxs[j][0]),2) + pow((m_chain.getLink(i)->getCenter()[1] - cxs[j][1]),2) + pow((m_chain.getLink(i)->getCenter()[2] - cxs[j][2]),2));
          if (dist <= 8.0)
          {
             //cout << "Nearest residues "  << i << endl;
             //res_list.push_back(m_link[i]);
             res_index.push_back(i);
             //cout << res_index[l] << endl;
             l++;
             //exit(0);
          }
       }
    }
  //return CNode::findClash(res_list);
}

//=================================================================================================
// find_clash
//=================================================================================================
void MonteCarloLigand::undoLastMove(vector<LPOSE> & accept_coords)
{
  accept_coords.pop_back();   
}
/*
void MonteCarloLigand::findClash()
{
  assert(!isLeaf());

  int time = m_pSL->getTime();
  CNode * curr = goDown();

  // Find clashes inside the first child.
  if (curr->isAffected(time) && curr->findSelfClash())
    return true;

  CNode * next = curr->getNext();
  if (next)
    {
      // Find clashes between the two children.
      if (curr->findPairClash(next, curr->m_rotate,
                              curr->m_translate,
                              false))
        return true;

      // Find clashes inside the second child.
      if (next->isAffected(time) && next->findSelfClash())
        return true;
    }

  return false;
 
}
*/
