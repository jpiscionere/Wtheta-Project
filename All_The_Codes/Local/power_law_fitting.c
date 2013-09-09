#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "chisquared_function.h"
#include "chisquared_function.c"
#include "chi2.h"


#define MAXLEN (100000)

int main(){


	double chi2(int lum_sample,double amp,double slope);	

	double chisquared;

	chisquared=chi2(3,152.81658748,2.05509656); 
	fprintf(stderr,"chisquared = %lf\n",chisquared);


return 0;

}



double chi2(int lum_sample,double amp, double slope)
{

	fprintf(stderr,"IN CHI2 FUNCTION\n");

	double *datatheta,*error,*datawtheta,*theory,**covariance_matrix;
	char	data_directory[MAXLEN],data_file[MAXLEN],covariance_file[MAXLEN];
	FILE	*fp1,*fp2;

	int bins2=100;
	
	int     i,j,k,l=0,diagonal_errors=0;
	double chisquared=100000;
	fprintf(stderr,"1\n");

	datatheta=(double *) calloc(bins2,sizeof(double)) ;
        error=(double *) calloc(bins2,sizeof(double)) ;
        datawtheta=(double *) calloc(bins2,sizeof(double)) ;
	theory=(double *) calloc(bins2,sizeof(double));	


	snprintf(data_directory,MAXLEN,"/home/piscioja/Clustering/WthetaPaper/Data");
        snprintf(data_file,MAXLEN,"%s/Wtheta/Wtheta_vollim_Mr%d_fib0.20rand.overlap.short",data_directory,lum_sample);



	fp1=fopen(data_file,"r") ;
        assert(fp1 != NULL );

	i=0;
	assert( i!=0);

        while(fscanf(fp1,"%lf %lf %lf",&datatheta[l],&datawtheta[l],&error[l])!=EOF)
       {
             l++ ;
       }

	bins2=l;

        fclose(fp1) ;
	


	covariance_matrix = malloc(bins2 * sizeof(double *));

        for (i = 0; i < bins2; i++)
                {
                         covariance_matrix[i] = malloc(bins2 * sizeof(double *));
                 }


	snprintf(covariance_file,MAXLEN,"%s/Wtheta/CovarWtheta_vollim_Mr%d_fib0.20rand.overlap.short.inv",data_directory,lum_sample);	
	fp2=fopen(covariance_file,"r");
        assert(fp2!= NULL);
        l=0;


	for(i=0;i<bins2;i++)
                        for(j=0;j<bins2;j++)
                        {
                                fscanf(fp2,"%lf ",&covariance_matrix[i][j]);

                        }

	fclose(fp2) ;
	for(i=0;i<bins2;i++) //loop over bins
                        {

                                theory[i] = amp*pow(datatheta[i]/datatheta[(int)bins2/2],slope); //theory is a power law
                        }



                  chisquared=chisquared_function(bins2,covariance_matrix, theory, datawtheta, error, diagonal_errors);


return chisquared;

}
