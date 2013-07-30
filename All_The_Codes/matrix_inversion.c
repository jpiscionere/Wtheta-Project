/*

covariance_fitting3 blantoncorrectedfile uncorrectedfile randomsfile modelfiles >output

*/



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "/home/piscioja/NRC/include/nr.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "/home/piscioja/Clustering/WthetaPaper/Source/inv.h"

#define sigma (5.67E-5)
#define NMAX (500)
#define sqr(x) ((x)*(x))

int main(int argc, char *argv[])
{

int i,j,k,l;
double **covariance_matrix,trash;
double *covar;
int bins = 20;

	sscanf(argv[1],"%d",&bins);


//Allocating for covariance_matrix	
	covariance_matrix = malloc(bins * sizeof(double *));
    	covar=(double * ) calloc(bins*bins,sizeof(double));	

	for (i = 0; i < bins; i++)
     		{
			 covariance_matrix[i] = malloc(bins * sizeof(double *));
        	 }


	while(fscanf(stdin,"%lf %lf %lf",&trash,&trash,&covar[l])!=EOF)
                {

                        l++;
                }

        
	for(j=0;j<bins;j++)
       		for(i=0;i<bins;i++)
                	{
                        	covariance_matrix[j][i]=covar[k];
                       		 k++;
                	}


	 inv(covariance_matrix,bins);

         for(i=0;i<bins;i++){
                for(j=0;j<bins;j++)
                        {
                                fprintf(stdout,"%lf ",covariance_matrix[i][j]);

                        }
                fprintf(stdout,"\n");
        }



	return 0;
}














