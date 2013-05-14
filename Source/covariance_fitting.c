/*

covariance_fitting datafile mock_errorfile stomp_format? invert? variance? diagonalerrors? start_bin end_bin <covariance_matvariancechisquared matrix 
 
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

int main(int argc, char *argv[])
{

   	int     i,j,k,l,input_format,use_variance,stomp_format,diagonal_errors;
        FILE    *fp1,*fp2;
        char    *data_file,*variance_file;
        void    inv (double**,int),filter(double**,int,int);
	double 	chisquared_function (int , double**, double*, double*, double*, double);
        double  **covariance_matrix;
        int     size_y=100,size_x=100,Nmax=100,N_file;
        double  *datatheta,*datawtheta;
        double  slope0,delta_slope,amplitude0, delta_amplitude,n_slope,n_amplitude;
        double  chisquared, chisquaredmin=1.E10, slope_best=-1000.0,amplitude_best=-1000.0;
        double  slope, amplitude,fake_chisquared,fake_chisquared_min=10000.0;
        double  *intermediate_step,trash;
        int     iarg,bins2=20,dof,start_bin,end_bin,bins_final;
        double  *theory,*variance,*error,*covar,*diagonal;



        data_file = argv[1];
        variance_file = argv[2];
        sscanf(argv[3],"%d",&stomp_format);
        sscanf(argv[4],"%d",&input_format);
        sscanf(argv[5],"%d",&use_variance);
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
				fprintf(stderr,"covariance > I expect YOU to have filtered data & variance files\n");
			}	
		if(bins2 < bins_final)
                        {
                                fprintf(stderr,"Bins total less than filtered bins! Bailing Out!\n");
                                return -1;
                        }
        sscanf(argv[10],"%lf",&slope0); //-1.4
        sscanf(argv[11],"%lf",&delta_slope); // 0.00075
        sscanf(argv[12],"%lf",&n_slope);//2000
                fprintf(stderr,"Your slope range is %lf to %lf\n",slope0, slope0 + delta_slope * n_slope);

        sscanf(argv[13],"%lf",&amplitude0);//0
        sscanf(argv[14],"%lf",&delta_amplitude);//0.01
        sscanf(argv[15],"%lf",&n_amplitude);//2000
                fprintf(stderr,"Your amplitude range is %lf to %lf\n",amplitude0,amplitude0 + delta_amplitude * n_amplitude);

        iarg=16;

        covar=(double * ) calloc(bins2*bins2,sizeof(double));
        diagonal=(double * ) calloc(bins2,sizeof(double));
        datatheta=(double *) calloc(bins2,sizeof(double)) ;
        error=(double *) calloc(bins2,sizeof(double)) ;
        datawtheta=(double *) calloc(bins2,sizeof(double)) ;
        intermediate_step=(double *) calloc(bins2,sizeof(double)) ;
        variance=(double *) calloc(bins2,sizeof(double)) ;
	theory=(double *) calloc(bins2,sizeof(double));

    	fp1=fopen(data_file,"r") ;
	assert(fp1 != NULL );

        l=0;

        while(fscanf(fp1,"%lf %lf %lf",&datatheta[l],&datawtheta[l],&error[l])!=EOF)
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

	fp2=fopen(variance_file,"r") ;
	assert(fp2 != NULL );
        l=0;
	
        while(fscanf(fp2,"%lf ",&variance[l])!=EOF)
		{
             		l++ ;
		}
	
        N_file = l ;

        if((N_file != bins_final)&&(use_variance==1))
                {
                        fprintf(stderr,"covariance > BIN NUMBER MISMATCH!! You DIDN'T filter the VARIANCE FILE!\n");
                        fprintf(stderr,"covariance > BAILING OUT!\n");
			return -1;                
                }
        fclose(fp2) ;	



	covariance_matrix = malloc(bins2 * sizeof(double *));
    	
	for (i = 0; i < bins2; i++)
     		{
			 covariance_matrix[i] = malloc(bins2 * sizeof(double *));
        	 }
	
 
	if (stomp_format==1) {
		
		l=0;
		
		while(fscanf(stdin,"%lf %lf %lf",&trash,&trash,&covar[l])!=EOF)
		{
		
			l++;
		}
	
		k=0;
		for(j=0;j<bins2;j++)
		{
			for(i=0;i<bins2;i++)
			{
				
				covariance_matrix[j][i]=covar[k];	
				k++;
			}
			diagonal[j]=covariance_matrix[j][j];
		}

		for(j=0;j<bins2;j++){
			for(i=0;i<bins2;i++){
				
				covariance_matrix[j][i] /= pow((diagonal[i] * diagonal[j]),0.5);
				
			}
			
		}

	}
	
	
	else{
		for(i=0;i<bins2;i++)
			for(j=0;j<bins2;j++)
			{
				fscanf(stdin,"%lf ",&covariance_matrix[i][j]);

			}

	}


// Filter Covar Matrix 
if(bins2 != bins_final)
        filter(covariance_matrix,start_bin,bins_final);
bins2=bins_final;

dof=bins2-2;

//End of filter
	
if(input_format == 1)
		inv(covariance_matrix,bins2);	







		
//Chisquared loop
	for(k=0;k<n_amplitude;k++)
	{	
		amplitude = amplitude0 + delta_amplitude*k; //loop over amplitudes
		
		
		for(l=0;l<n_slope;l++)
		{	
			slope = slope0 + delta_slope*l; //loop over slopes
			chisquared = 0;	//reset value for each combination of slope and amplitude
	
			
			for(i=0;i<bins2;i++) //loop over bins
			{
      			
				theory[i] = amplitude*pow(datatheta[i]/datatheta[(int)bins2/2],slope); //theory is a power law
			}	
	
			if(use_variance==0)
				chisquared=chisquared_function(bins2,covariance_matrix, theory, datawtheta, error, diagonal_errors);
			else
				chisquared=chisquared_function(bins2,covariance_matrix, theory, datawtheta, variance, diagonal_errors);

			fprintf(stdout,"%lf ",chisquared);	
			
			if(chisquared < chisquaredmin)
			{
				slope_best = slope;
				amplitude_best = amplitude;
				chisquaredmin = chisquared; 
			}
			
			
			
		}	
		
		fprintf(stdout,"\n");	
	}

if((slope_best == slope) | (amplitude_best == amplitude) | (slope_best == slope0) | (amplitude_best == amplitude0))
        {
                fprintf(stderr,"Bad Values for Slope or Amplitude\n");
        
        }



fprintf(stderr,"%s> Best Chisquared = %lf for Slope =%lf and Amplitude=%lf\n",data_file,chisquaredmin/dof,slope_best,amplitude_best);


free(covar);
free(diagonal);
free(datatheta);
free(error);
free(datawtheta);
free(variance);
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
