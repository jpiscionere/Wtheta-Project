/*

covariance_fitting_sigma size_x size_y x0 dx y0 dy  dimension? < chisquared_matrix > 1-D Pdf distribution

*/




#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "/home/piscioja/NRC/include/nr.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#define sigma (5.67E-5)
#define NMAX (5000)
#define sqr(x) ((x)*(x))


int main(int argc, char *argv[])
{


int i,j,k,l,m,x,y,z,slope_min,amp_min;
double **chi_squared,*prob_sum,sum=0;
int size_x,size_y;
int dimension;
double x0,dx,y0,dy;
double chisquared_min=100000.0;

sscanf(argv[1],"%d",&size_x);
sscanf(argv[2],"%d",&size_y);
sscanf(argv[3],"%lf",&x0);
sscanf(argv[4],"%lf",&dx);
sscanf(argv[5],"%lf",&y0);
sscanf(argv[6],"%lf",&dy);
sscanf(argv[7],"%d",&dimension);



fprintf(stderr,"Dimension 1 runs from %lf to %lf\n",x0,x0 + size_x*dx);
fprintf(stderr,"Dimension 2 runs from %lf to %lf\n",y0,y0 + size_y*dy);
fprintf(stderr,"You are integrating over Dimension %d\n",dimension);
fprintf(stderr,"Your output will be the 1-D probablility distribution of Dimension %d\n",dimension);

chi_squared = malloc(size_x * sizeof(double *));

for (i = 0; i < size_y; i++)
{
	chi_squared[i] = malloc(size_x * sizeof(double *));
}


// everything is backwards.  x - actualy y dimen and y is actually x dimension

for(y=0;y<size_y;y++)
	for(x=0;x<size_x;x++)
		{

			fscanf(stdin,"%lf ",&chi_squared[y][x]);
			if(chi_squared[y][x] < chisquared_min)
				{
					chisquared_min = chi_squared[y][x];
					slope_min=x;
					amp_min=y;
				}
		}


fprintf(stderr,"Chisquared Min = %lf, for slope %lf and amp %lf\n",chisquared_min/18.0,y0+dy*slope_min,x0+dx*amp_min);



if(dimension==1) //integrating out the x dimension (amplitude)
{


double min_y = 100.0; //one bigger than the biggest possible
double max_y = -1000.0;
double value_y;
prob_sum=(double *) calloc(size_y,sizeof(double)) ;


//chi_squared[amp][slope]
//chi_squared[realx][realy]

i=0;
for(x=0;x<size_x;x++)
	for(y=0;y<size_y;y++)
	{


		prob_sum[x] += exp( - 1. * chi_squared[y][x] );
		sum += prob_sum[x];
		if( (chi_squared[y][x] - chisquared_min) <= 2.3)
		{
			value_y = y0 + dy*x;

			i++;
			if(value_y < min_y)
				min_y=value_y;
			if(value_y > max_y)
				max_y=value_y;
			
		}
		
	}












//fprintf(stderr,"1-sigma error on y dimension %lf to %lf for a range of %lf\n",max_y,min_y,max_y-min_y);
fprintf(stderr,"%lf %lf %lf\n",y0+dy*slope_min,max_y,min_y);
fprintf(stdout,"%lf %lf %lf\n",y0+dy*slope_min,max_y,min_y);


for(x=0;x<size_x;x++)
	fprintf(stdout,"%lf %e\n",y0 + dy*x,prob_sum[x]/sum);

free(prob_sum);

}

free(chi_squared);


return 0;
}



