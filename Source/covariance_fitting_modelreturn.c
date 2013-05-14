/*

covariance_fitting datafile mock_errorfile stomp_format? invert? variance? diagonalerrors? start_bin end_bin <covariance_matrix >chisquared matrix 
 
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
#include "chisquared_function.h"
#include "filter.h"
#include "inv.h"
#define sigma (5.67E-5)
#define NMAX (500)
#define MAXLEN (1000)
#define sqr(x) ((x)*(x))

double main(int argc, char *argv[])
{

	int 	i,j,k,l,input_format,use_variance,stomp_format,diagonal_errors,count=0;
	FILE 	*fp1,*fp2,*fp3,*fp4;
	char 	*data_file,*variance_file,*randoms_file;
	char	inpath[MAXLEN],model_file[MAXLEN],prefix[MAXLEN];
	void 	inv (double**,int), filter(double**,int,int);
	double  chisquared_function (int , double**, double*, double*, double*, double);
	double 	**covariance_matrix;
	int 	size_y=100,size_x=100,Nmax=100,N_file;
	double 	*datatheta,*datawtheta,*theory_wtheta;	
	double 	gamma0,delta_gamma,fgal0, delta_fgal,n_gamma,n_fgal;
	double  chisquared, chisquaredmin=1.E10, gamma_best=-1000.0,fgal_best=-1000.0;
	double 	theory,gamma=0,fgal=0;
	double 	*intermediate_step,trash;
	int 	iarg,models_total,bins2=20,dof,start_bin,end_bin,bins_final;
	double 	*variance,*error,*covar,*diagonal,*RR;




	data_file = argv[1];
        variance_file = argv[2];
	randoms_file = argv[3];
	sscanf(argv[4],"%d",&stomp_format);
	sscanf(argv[5],"%d",&input_format);
        sscanf(argv[6],"%d",&use_variance);
        sscanf(argv[7],"%d",&diagonal_errors);
	sscanf(argv[8],"%d",&bins2);
	sscanf(argv[9],"%d",&start_bin);
	sscanf(argv[10],"%d",&end_bin);
		if(end_bin < 1)
			end_bin=bins2;
		bins_final=end_bin - start_bin;
		if(bins2 != bins_final)
			fprintf(stderr,"covariance > You have decided to filter out %d bins\n",bins2 - bins_final);
		if(bins2 < bins_final)
			{
				fprintf(stderr,"covariance > Bins total less than filtered bins! Bailing Out!\n");
				return -1;
			}

	snprintf(inpath, MAXLEN,"%s",argv[11]);
	snprintf(prefix, MAXLEN,"%s",argv[12]);

        covar=(double * ) calloc(bins2*bins2,sizeof(double));
        diagonal=(double * ) calloc(bins2,sizeof(double));
        datatheta=(double *) calloc(bins2,sizeof(double)) ;
        error=(double *) calloc(bins2,sizeof(double)) ;
        datawtheta=(double *) calloc(bins2,sizeof(double)) ;
        intermediate_step=(double *) calloc(bins2,sizeof(double)) ;
        theory_wtheta=(double *) calloc(bins2,sizeof(double)) ;
        RR=(double *) calloc(bins2,sizeof(double)) ;
	variance=(double *) calloc(bins2,sizeof(double));



	fp1=fopen(data_file,"r") ;
	assert( fp1 != NULL );	
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
	assert (fp2 != NULL );

        l=0;

        while(fscanf(fp2,"%lf ",&variance[l])!=EOF)
       		{
             		l++ ;
		}

        N_file = l ;

	if((N_file != bins_final) && (use_variance == 1))
                {
                        fprintf(stderr,"covariance > BIN NUMBER MISMATCH!! You DIDN'T filter the VARIANCE FILE!\n");
                        fprintf(stderr,"covariance > BAILING OUT!\n");
                        return -1;
                }
       
        fclose(fp2) ;
        fp3=fopen(randoms_file,"r") ;
	assert( fp3 != NULL ); 

        l=0;

        while(fscanf(fp3,"%lf ",&RR[l])!=EOF)
                {
                        l++ ;
                }

        N_file = l ;

        if(N_file != bins_final)
                {
                        fprintf(stderr,"covariance > BIN NUMBER MISMATCH!! You DIDN'T filter %d bins of the RANDOMS FILE!\n",N_file - bins_final);
                        fprintf(stderr,"covariance > BAILING OUT!\n");
                        return -1;
                }

        fclose(fp3) ;

	
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

	snprintf(model_file,MAXLEN,"%s/Wtheta_%s.average.f%5.4lf.g%5.4lf.wtheta",inpath,prefix,fgal,gamma);	
	fp4=fopen(model_file,"r");
	assert(fp4 != NULL);	
	
	while(fscanf(fp4,"%lf",&theory_wtheta[count])!=EOF)
	{
		theory_wtheta[count] = theory_wtheta[count]/RR[count] - 1.0;		
		count++;
        }
			
	fclose(fp4) ;			
	assert(count == bins2);			

	if(use_variance==0)
		chisquared=chisquared_function(bins2,covariance_matrix, theory_wtheta, datawtheta, error, diagonal_errors);
	else
		chisquared=chisquared_function(bins2,covariance_matrix, theory_wtheta, datawtheta, variance, diagonal_errors);

free(covar);
free(diagonal);
free(datatheta);
free(error);
free(datawtheta);
free(intermediate_step);
free(theory_wtheta);
free(RR);
free(variance);
free(covariance_matrix);
return(chisquared);
}

