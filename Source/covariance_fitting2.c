/*

covariance_fitting2 datafile randomsfile args  modelfiles >output


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

	int i,j,k,l,input_format,iarg,arg_number,diagonal_errors;
	FILE *fp1,*fp2,*fp3;
	char *data_file,*randoms_file,*model_file;
	void inv (double**,int);
	double  chisquared_function (int , double**, double*, double*, double*, double);
	double **covariance_matrix;
	int Nmax=100,N_file;
	double *datatheta,*datawtheta;
	double slope0=-1.4,delta_slope=0.00075,n_slope,n_amplitude,amplitude0=0, delta_amplitude=0.01, chisquared, chisquaredmin=1.E10, slope_best=-1000.0,amplitude_best=-1000;
	double slope, amplitude,junkf,junk;
	double *intermediate_step;
	int bins2=20,dof = 18,print_variance,junki;
	double *theory,wtheta1,*wtheta_average,*variance,*randoms,*Delta,*error;




	data_file = argv[1];
        randoms_file = argv[2];
	sscanf(argv[3],"%d",&bins2);
	sscanf(argv[4],"%d",&print_variance);
        sscanf(argv[5],"%d",&diagonal_errors);
	sscanf(argv[6],"%lf",&slope0); //-1.4
        sscanf(argv[7],"%lf",&delta_slope); // 0.00075
        sscanf(argv[8],"%lf",&n_slope);//2000
                fprintf(stderr,"Your slope range is %lf to %lf\n",slope0, slope0 + delta_slope * n_slope);

        sscanf(argv[9],"%lf",&amplitude0);//0
        sscanf(argv[10],"%lf",&delta_amplitude);//0.01
        sscanf(argv[11],"%lf",&n_amplitude);//2000
                fprintf(stderr,"Your amplitude range is %lf to %lf\n",amplitude0,amplitude0 + delta_amplitude * n_amplitude);
	
	arg_number=12;
	
        datatheta=(double *) calloc(bins2,sizeof(double)) ;
        datawtheta=(double *) calloc(bins2,sizeof(double)) ;
        intermediate_step=(double *) calloc(bins2,sizeof(double)) ;
        wtheta_average=(double *) calloc(bins2,sizeof(double)) ;
        variance=(double *) calloc(bins2,sizeof(double)) ;
        randoms=(double *) calloc(bins2,sizeof(double)) ;
        Delta=(double *) calloc(bins2,sizeof(double)) ;
        error=(double *) calloc(bins2,sizeof(double)) ;
	theory=(double *) calloc(bins2,sizeof(double));
	

	dof=bins2-2.;

	
//Reading in Data File and Storing Errors

	fp1=fopen(data_file,"r") ;
	assert(fp1 != NULL );	
	l=0;
	
    	while(l < bins2){
		fscanf(fp1,"%lf %lf %lf",&datatheta[l],&datawtheta[l],&error[l]);
        	fprintf(stderr,"%lf %lf %lf\n",datatheta[l],datawtheta[l],error[l]);
		l++ ;
       	}

	N_file = l ;
  
	if(N_file > Nmax)
		fprintf(stderr,"File too long, adjust parameters!\n");
	else
		fprintf(stderr,"Data File Read in Successfully.\n");
	fclose(fp1) ;
//Reading in Randoms File


     	fp2=fopen(randoms_file,"r") ;
	assert(fp2 != NULL );
        l=0;

        while(fscanf(fp2,"%lf ",&randoms[l])!=EOF){
             l++ ;
       }

        N_file = l ;

        if(N_file > bins2)
                fprintf(stderr,"File too long, adjust parameters!\n");
	else
		fprintf(stderr,"Randoms File Read in Successfully.\n");
        fclose(fp2) ;

//Allocating for covariance_matrix	
	covariance_matrix = malloc(bins2 * sizeof(double *));
    	
	for (i = 0; i < bins2; i++)
     		{
			 covariance_matrix[i] = malloc(bins2 * sizeof(double *));
        	 }



//Wtheta Average Loop
for(iarg=arg_number;iarg<argc;iarg++)
	{
		
		
	     	model_file = argv[iarg];
    		fprintf(stderr, "covariance> Reading Model file: %s\n", model_file);
		
		fp3=fopen(model_file,"r") ;
		assert(fp3 != NULL );
		l=0;
		while(fscanf(fp3,"%d %lf %lf %d ",&i,&junk,&wtheta1,&i)!=EOF)
       			{
				wtheta_average[l]+= (wtheta1/randoms[l] - 1.0);		
				l++;
				
			}
		fclose(fp3) ;	

	}
//average the average
for(i=0;i<bins2;i++)
	{
		
		wtheta_average[i]= wtheta_average[i]/(argc - arg_number);
		fprintf(stderr,"Wtheta Model = %lf, Wtheta Data = %lf\n",wtheta_average[i],datawtheta[i]);
	}


//Variance Loop
	for(iarg=arg_number;iarg<argc;iarg++)
	{
		
		
		model_file = argv[iarg];
		
		
		fp3=fopen(model_file,"r") ;
		l=0;
		while(fscanf(fp3,"%d %lf %lf %d ",&i,&junk,&wtheta1,&i)!=EOF)
		{
			variance[l]+=sqr((wtheta1/randoms[l] - 1.0)-wtheta_average[l]);
			l++;
		}
		fclose(fp3) ;
		
		
	}
	
	//Normalizing the Error from the Mocks 
	for(i=0;i<bins2;i++)
	{
		
		variance[i]= pow(variance[i]/(argc - arg_number - 1.0),0.5);
		fprintf(stderr,"%lf\n",variance[i]);	

	}

	if(print_variance==1){
		for(i=0;i<bins2;i++)
		{
			fprintf(stdout,"%lf\n",variance[i]);        
			
		}
		
		return 0;
	}


//Delta and Covariance Loop
for(iarg=arg_number;iarg<argc;iarg++)
        {


                model_file = argv[iarg];	

                fp3=fopen(model_file,"r") ;
                l=0;
                while(fscanf(fp3,"%d %lf %lf %d ",&i,&junk,&wtheta1,&i)!=EOF)
                        {
                               	Delta[l] = ((wtheta1/randoms[l] - 1.0) - wtheta_average[l])/variance[l]; //Fill in the Delta Array for EACH mock 
                        
				l++;
                        }
                fclose(fp3) ;

		for(i=0;i<bins2;i++)
			for(j=0;j<bins2;j++)
				{
					covariance_matrix[i][j]+= Delta[i]*Delta[j];//Fill in the cumulative covariance matrix
		
				}



	}


if(print_variance != 1){
	for(i=0;i<bins2;i++)
		{
		for(j=0;j<bins2;j++)
			{
				covariance_matrix[i][j]= covariance_matrix[i][j]/(argc - arg_number -1.0); //Normalize Covariance Matrix by number of mocks
				fprintf(stdout,"%lf ",covariance_matrix[i][j]);		
                	}
                        
			fprintf(stdout,"\n");
		}
}else{
	for(i=0;i<bins2;i++)
	{
        	for(j=0;j<bins2;j++)
                	{
                                covariance_matrix[i][j]= covariance_matrix[i][j]/(argc - arg_number -1.0); //Normalize Covariance Matrix by number of mocks
                       
                        }
 
                }

}

inv(covariance_matrix,bins2);

//Chisquared loop
        for(k=0;k<n_amplitude;k++)
        {
                amplitude = amplitude0 + delta_amplitude*k; //loop over amplitudes


                for(l=0;l<n_slope;l++)
                {
                        slope = slope0 + delta_slope*l; //loop over slopes
                        chisquared = 0; //reset value for each combination of slope and amplitude


                        for(i=0;i<bins2;i++) //loop over bins
                        {

                                theory[i] = amplitude*pow(datatheta[i]/datatheta[(int)bins2/2],slope); //theory is a power law
                        }

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
                fprintf(stderr,"Bad Values for Slope or Amplitude, slope = %lf amp=%lf\n",slope_best,amplitude_best);
                return (-1);
        }

fprintf(stderr,"%s>Best Chisquared = %lf for Slope =%lf and Amplitude=%lf\n",data_file,chisquaredmin/dof,slope_best,amplitude_best);


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

