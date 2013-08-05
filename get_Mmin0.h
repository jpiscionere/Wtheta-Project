/* PROGRAM get_Mmin0
   
   --- get_Mmin0 ngal siglogM logM0 logM1 alpha BGCfile < Halocenters > logMmin
   --- Computes Mmin required to get a desired space density for a particular halo list and HOD.
   --- Version 0 uses the BGC halo format and the Zheng, Coil, Zehavi 2007 HOD model).

      * ngal    = desired space density of galaxies (in h^3/Mpc^3)
      * siglogM = width of cutoff at Mmin (scatter in the mass-luminosity relation)
      * logM0   = minimum halo mass that can contain a satellite galaxy (in units of Msun/h)
      * logM1   = halo mass for which <Nsat>=1 (in units of Msun/h)
      * alpha   = slope of the <Nsat> - M relation
      * BGCfile (BGC format)
      * < Halo center file (output of calc_group_stats.c)
      * > logMmin
*/
#ifndef _GET_MMIN0_H
#define _GET_MMIN0_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "/home/piscioja/NRC/include/nr.h"
#include "/home/piscioja/halobias/bgc_read_utils.c"

#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define cnint(x) ((x-floor(x)) < 0.5 ? floor(x) : ceil(x))
#define csign(x) (x < 0.0 ? -1 : 1)
#define PI (3.141592)
#define epsilon 1e-9

/*
int main()
{

double get_Mmin(double , double , double , double , double , char* , char*);
double logMmin;
char bgc_file[MAXLEN],halos_file[MAXLEN];

snprintf(bgc_file,MAXLEN,"/net/bender/data0/LasDamas/Consuelo/4004/fof_b0p2/Consuelo_4004_z0p054_fof_b0p2.0000.bgc");
snprintf(halos_file,MAXLEN,"/net/bender/data0/LasDamas/Consuelo/4004/Consuelo_4004_z0p054_fof_b0p2.com.halos");


logMmin=get_Mmin(0.3,0.1,11.0,12.5,0.1,bgc_file,halos_file);

fprintf(stderr,"get_Mmin0> logMmin = %5.2f\n",logMmin) ;

return 0;


}
*/


//function declaration
double get_Mmin0(double rhogal, double siglogM, double logM0, double logM1, double alpha, char *bgc_file, char *mf_file);


#endif




