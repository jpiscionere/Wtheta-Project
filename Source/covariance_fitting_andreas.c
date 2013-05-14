/*

covariance_fitting datafile theoryfile invert? diagonalerrors? start_bin end_bin <covariance_matvariancechisquared matrix 
 
*/ 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "/home/piscioja/NRC/include/nr.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#define sigma (5.67E-5)
#define NMAX (500)
#define sqr(x) ((x)*(x))


void Printhelp(void)
{
  fprintf(stderr, "%s", "\n"
  "covariance_fitting_andreas datafile theoryfile number_density error_number_density theory_number_density invert_covar? diagonalerrors? start_bin end_bin < covariance_matrix > chisquared\n"
  "--data file in 6 column format\n"
  "--theory file in 1 column format\n"
  "--number density\n"
  "--error on number density\n"
  "--theory number density\n"
  "--invert_covar	= 0 : no\n"
  "                 	= 1 : yes\n"
  "--diagonalerrors 	= 0 : use covar for fit\n"
  "		    	= 1 : use diagonal errors for fit\n"
  "--total_bins     	= total number of bins, normal counting\n"
  "--start_bin	    	= first bin you care about, counting starts at zero\n"
  "--end_bin        	= last bin you care about, counting starts at zero\n"
  "< covariance matrix\n" 	
 );
}



int main(int argc, char *argv[])
{

   	int     i,j,l,input_format,diagonal_errors;
        FILE    *fp1,*fp2;
        char    *data_file,*theory_file;
        void    inv (double**,int),filter(double**,int,int);
	double 	chisquared_function (int , double**, double*, double*, double*, double);
        double  **covariance_matrix;
        int     N_file;
        double  *datar,*datawp;
        double  chisquared;  
        double  *intermediate_step,fjunk;
        int     iarg,bins2=14,start_bin,end_bin,bins_final;
        double  *theory,*error,*covar,*diagonal;
	double  num_den,err_num_den,theory_num_den;

	if(argc<9)
	 {
     		 Printhelp() ;
      		 return -1 ;
    	 }
  
	
        data_file = argv[1];
        theory_file = argv[2];
	sscanf(argv[3],"%lf",&num_den);
	sscanf(argv[4],"%lf",&err_num_den);
	sscanf(argv[5],"%lf",&theory_num_den);
        sscanf(argv[6],"%d",&diagonal_errors);
        sscanf(argv[7],"%d",&bins2);
        sscanf(argv[8],"%d",&start_bin);
        sscanf(argv[9],"%d",&end_bin);
                if(end_bin < 1)
                        end_bin=bins2;
                bins_final=end_bin - start_bin;
                if(bins2 != bins_final)
			{
                        	fprintf(stderr,"covariance > You have decided to filter out %d bins of the covariance matrix\n",bins2 - bins_final);
                		fprintf(stderr,"covariance > WARNING!!! Data file filtering not supported\n");
				fprintf(stderr,"covariance > I expect YOU to have filtered data & theory files\n");
			}	
		if(bins2 < bins_final)
                        {
                                fprintf(stderr,"Bins total less than filtered bins! Bailing Out!\n");
                                return -1;
                        }
        iarg=10;

        covar=(double * ) calloc(bins2*bins2,sizeof(double));
        diagonal=(double * ) calloc(bins2,sizeof(double));
        datar=(double *) calloc(bins2,sizeof(double)) ;
        error=(double *) calloc(bins2,sizeof(double)) ;
        datawp=(double *) calloc(bins2,sizeof(double)) ;
        intermediate_step=(double *) calloc(bins2,sizeof(double)) ; 
	theory=(double *) calloc(bins2,sizeof(double));

    	fp1=fopen(data_file,"r") ;
	assert(fp1 != NULL );

        l=0;

        while(fscanf(fp1,"%lf %lf %lf %lf %lf",&fjunk,&fjunk,&datar[l],&datawp[l],&error[l])!=EOF)
       {
             l++ ;
       }

        N_file = l ;

	if(N_file != bins_final)
		{
			fprintf(stderr,"covariance > BIN NUMBER MISMATCH!! You DIDN'T filter the DATA FILE!\n");
			fprintf(stderr,"covariance > BAILING OUT!\n");
			return -1;		
		}

        fclose(fp1) ;

	fp2=fopen(theory_file,"r") ;
	assert(fp2 != NULL );
        l=0;
	
        while(fscanf(fp2,"%lf %lf ",&fjunk,&theory[l])!=EOF)
		{
             		l++ ;
		}
	
        N_file = l ;

        if(N_file != bins_final)
                {
                        fprintf(stderr,"covariance > BIN NUMBER MISMATCH!! You DIDN'T filter the THEORY FILE!\n");
                        fprintf(stderr,"covariance > BAILING OUT!\n");
			return -1;                
                }
        fclose(fp2) ;	



	covariance_matrix = malloc(bins2 * sizeof(double *));
    	
	for (i = 0; i < bins2; i++)
     		{
			 covariance_matrix[i] = malloc(bins2 * sizeof(double *));
        	 }
	
	
	for(i=0;i<bins2;i++)
		for(j=0;j<bins2;j++)
			{
				fscanf(stdin,"%lf ",&covariance_matrix[i][j]);
				
			}



	// Filter Covar Matrix 
	if(bins2 != bins_final)
        	filter(covariance_matrix,start_bin,bins_final);
	bins2=bins_final;	

//End of filter
	
	if(input_format == 1)
		inv(covariance_matrix,bins2);	


	chisquared=chisquared_function(bins2,covariance_matrix, theory, datawp, error, diagonal_errors);
	chisquared += sqr((num_den-theory_num_den)/(err_num_den));

	fprintf(stderr,"%s> Chisquared = %lf \n",data_file,chisquared);
	fprintf(stdout,"%lf \n",chisquared);


	free(covar);
	free(diagonal);
	free(datar);
	free(error);
	free(datawp);
	free(covariance_matrix);

return(0);
}

void inv(double **a_matrix, int dimension){
	
	
	int i,j;
	
	
	
	gsl_matrix * C = gsl_matrix_calloc ( dimension, dimension );
	gsl_matrix * C_i = gsl_matrix_calloc ( dimension, dimension );
	
	for ( i=0;i<dimension;i++ )
		for ( j=0;j<dimension;j++ ){
			
			gsl_matrix_set ( C, i, j, a_matrix[i][j] );
		}
	
	
	int s;
	
	gsl_permutation *p = gsl_permutation_alloc ( dimension );
	
	gsl_linalg_LU_decomp ( C, p, &s );
	
	gsl_linalg_LU_invert ( C, p, C_i );
	
	for ( i=0;i<dimension;i++ )
  	{
    	for ( j=0;j<dimension;j++ ) 
		{
			a_matrix[i][j] = gsl_matrix_get( C_i, i, j ) ;
			
  		}
	}
	
	
	gsl_permutation_free( p );
	gsl_matrix_free( C );
	gsl_matrix_free( C_i );
	
	
}

void filter (double **matrix, int start_bin, int bins_final)
{

        double **matrix_filter;
        int i,j;

        matrix_filter = malloc(bins_final * sizeof(double *));

        for (i = 0; i < bins_final; i++)
                {
                         matrix_filter[i] = malloc(bins_final * sizeof(double *));
                 }




        for(i=0;i<bins_final;i++)
                for(j=0;j<bins_final;j++)
                {
                        matrix_filter[j][i]=matrix[i+start_bin][j+start_bin];
                }


        for(i=0;i<bins_final;i++)
                for(j=0;j<bins_final;j++)
                        matrix[j][i]=matrix_filter[j][i];

free(matrix_filter);


}

double chisquared_function (int dimen, double **covar, double *theory, double *data, double *error, double format)
{

double *intermediate_step,chisquared=0;
int i,j;

intermediate_step = (double * ) calloc(dimen,sizeof(double));

for(i=0;i<dimen;i++)
{
        intermediate_step[i] = (data[i]-theory[i])/(error[i]);

}

if(format==0)
{
        for(i=0;i<dimen;i++)
                {
                for(j=0;j<dimen;j++)
                        {
                                chisquared+=intermediate_step[i]*covar[i][j]*intermediate_step[j];
                        }
                }

}else{
        for(i=0;i<dimen;i++)
                        {
                                chisquared+=intermediate_step[i]*intermediate_step[i]; //cumulative chisquared using another matrix

                        }

}

free(intermediate_step);

return chisquared;


} 
