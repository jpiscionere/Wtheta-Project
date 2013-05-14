/*

covariance_fitting3 datafile randomsfile outputfile gamma fgal modelfiles >output

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
#define NMAX (500)
#define sqr(x) ((x)*(x))

int main(int argc, char *argv[])
{

int i,j,k,l,input_format,iarg=6,arg_number=6;
FILE *fp1,*fp2,*fp3,*fp4;
char *data_file,*randoms_file,*model_file,*output_file;
void legs(double **,int ,double *,double *,int *);
void elgs (double** ,int ,int* );
void inv (double**,int);
double **covariance_matrix;
int size_y=100,size_x=100,Nmax=100,N_file;
double *datatheta,*datawtheta,chisquared;
double *intermediate_step,gamma,fgal;
int bins = 20,bins2=20,dof = 18;
double wtheta1,wtheta_average[bins],variance[bins],randoms[bins],Delta[bins],error[bins];

	sscanf(argv[4],"%lf",&gamma);
	sscanf(argv[5],"%lf",&fgal);

	
	datatheta=(double *) calloc(Nmax,sizeof(double)) ;
	datawtheta=(double *) calloc(Nmax,sizeof(double)) ;
	intermediate_step=(double *) calloc(Nmax,sizeof(double)) ;

	
//Reading in Data File and Storing Errors
	data_file = argv[1];
	fp1=fopen(data_file,"r") ;
	
	if(data_file == NULL)
		{
			fprintf(stderr,"Cannot Open %s\n",data_file) ;
			return (-1) ;
		}
	
	else
		{
			fprintf(stderr,"Opened %s\n",data_file) ;
       		}
	
	l=0;
	
	while(fscanf(fp1,"%lf %lf %lf",&datatheta[l],&datawtheta[l],&error[l])!=EOF)
       {
             l++ ;
       }

	N_file = l ;
  
	if(N_file > Nmax)
		fprintf(stderr,"File too long, adjust parameters!\n");


	fclose(fp1) ;
//Reading in Randoms File

     randoms_file = argv[2];
     fp2=fopen(randoms_file,"r") ;

        if(fp2 == NULL)
                {
                        fprintf(stderr,"Cannot Open %s\n",randoms_file) ;
                        return (-1) ;
                }

        else
		fprintf(stderr,"Opened %s\n",randoms_file) ;
      	

        l=0;

        while(fscanf(fp2,"%lf ",&randoms[l])!=EOF)
       {
             l++ ;
       }

        N_file = l ;

        if(N_file > Nmax)
                fprintf(stderr,"File too long, adjust parameters!\n");


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
		
		fp4=fopen(model_file,"r") ;
		l=0;
		while(fscanf(fp4,"%d %lf %lf ",&i,&j,&wtheta1)!=EOF)
       			{
				wtheta_average[l]+= (wtheta1/randoms[l] - 1.0);	
				l++;
				
			}
		fclose(fp4) ;	

	}
//average the average
for(i=0;i<bins;i++)
	{
		
		wtheta_average[i]= wtheta_average[i]/(argc - arg_number);

	}


//Variance Loop
for(iarg=arg_number;iarg<argc;iarg++)
        {


                model_file = argv[iarg];
             

                fp4=fopen(model_file,"r") ;
                l=0;
                while(fscanf(fp4,"%d %lf %lf ",&i,&j,&wtheta1)!=EOF)
                        {
				variance[l]+=sqr((wtheta1/randoms[l] - 1.0)-wtheta_average[l]);
                                l++;
                        }
                fclose(fp4) ;
                
               
           }

//Normalizing Variance
for(i=0;i<bins;i++)
        {

                variance[i]= pow(variance[i]/(argc - arg_number - 1.0),0.5);

        }

//Creating the bb matrix
for (i = 0; i < bins2; i++)
	{
        for (j = 0; j < bins2; j++)
                {
        
                	covariance_matrix[i][j] =0;
		}
        }



//Delta and Covariance Loop
for(iarg=arg_number;iarg<argc;iarg++)
        {


                model_file = argv[iarg];	

                fp4=fopen(model_file,"r") ;
                l=0;
                while(fscanf(fp4,"%d %lf %lf ",&i,&j,&wtheta1)!=EOF)
                        {
                               	Delta[l] = ((wtheta1/randoms[l] - 1.0) - wtheta_average[l])/variance[l]; //Fill in the Delta Array for EACH mock 
                        
				l++;
                        }
                fclose(fp4) ;

		for(i=0;i<bins2;i++)
			for(j=0;j<bins2;j++)
				{
					covariance_matrix[i][j]+= Delta[i]*Delta[j];//Fill in the cumulative covariance matrix
		
				}



	}


for(i=0;i<bins2;i++)
	{
	for(j=0;j<bins2;j++)
		{
			covariance_matrix[i][j]= covariance_matrix[i][j]/(argc - arg_number -1.0); //Normalize Covariance Matrix by number of mocks
			fprintf(stdout,"%lf ",covariance_matrix[i][j]);
			
			
                }
                        
		fprintf(stdout,"\n");

		
	}

inv(covariance_matrix,bins2);

		chisquared = 0;	//reset value for each combination of slope and amplitude

	
		for(i=0;i<bins2;i++) //loop over bins
    			{
      			
    			
				intermediate_step[i]= (datawtheta[i]-wtheta_average[i])/variance[i]; //Same as Deltas	

		
			}      			
			




			
		for(i=0;i<bins2;i++) //double loop to measure chisquared
   			{
      			for(j=0;j<bins2;j++)
        			{	
	       				chisquared+=intermediate_step[i]*covariance_matrix[i][j]*intermediate_step[j]; //cumulative chisquared for each combination of slope and amplitude using covariance matrix
       	
       				
				}
   			 }
		
			

    
//Write values to output file
	
	output_file = argv[3];
	fp3=fopen(output_file,"a") ;
	fprintf(fp3,"%lf %lf  %lf\n",gamma,fgal,chisquared/dof);
       	fclose(fp3) ;



return(0);
}

void legs(double **a_matrix,int n,double *b,double *x,int *indx)
{
	int i,j;
	void elgs();
	
	fprintf(stderr,"Legs has been called\n") ;
	
	elgs (a_matrix,n,indx);
	
	for(i = 0; i < n-1; ++i)
	{
		for(j = i+1; j < n; ++j)
		{
			b[indx[j]] = b[indx[j]]-a_matrix[indx[j]][i]*b[indx[i]];
		}
	}
	
	x[n-1] = b[indx[n-1]]/a_matrix[indx[n-1]][n-1];
	for (i = n-2; i>=0; i--)
	{
		x[i] = b[indx[i]];
		for (j = i+1; j < n; ++j)
		{
			x[i] = x[i]-a_matrix[indx[i]][j]*x[j];
		}
		x[i] = x[i]/a_matrix[indx[i]][i];
	}
}


void elgs(double **a_matrix,int n, int *indx)
{
  int i, j, k, itmp;
  double c1, pi, pi1, pj;
  double c_matrix[NMAX];

  if (n > NMAX)
  {
    printf("The matrix dimension is too large.\n");
    exit(1);
  }

/* Initialize the index */

  for (i = 0; i < n; ++i)
  {
    indx[i] = i;
  }

/* Find the rescaling factors, one from each row */

  for (i = 0; i < n; ++i)
  {
    c1 = 0;
    for (j = 0; j < n; ++j)
    {
      if (fabs(a_matrix[i][j]) > c1) c1 = fabs(a_matrix[i][j]);
    }
    c_matrix[i] = c1;
  }

/* Search the pivoting (largest) element from each column */

  for (j = 0; j < n-1; ++j)
  {
    pi1 = 0;
    for (i = j; i < n; ++i)
    {
      pi = fabs(a_matrix[indx[i]][j])/c_matrix[indx[i]];
      if (pi > pi1)
      {
        pi1 = pi;
        k = i;
      }
    }

/* Interchange the rows via indx[] to record pivoting order */

    itmp = indx[j];
    indx[j] = indx[k];
    indx[k] = itmp;
    for (i = j+1; i < n; ++i)
    {
      pj = a_matrix[indx[i]][j]/a_matrix[indx[j]][j];

/* Record pivoting ratios below the diagonal */

      a_matrix[indx[i]][j] = pj;

/* Modify other elements accordingly */

      for (k = j+1; k < n; ++k)
      {
        a_matrix[indx[i]][k] = a_matrix[indx[i]][k]-pj*a_matrix[indx[j]][k];
      }
    }
  }
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


