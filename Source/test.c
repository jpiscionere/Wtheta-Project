#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "/home/piscioja/NRC/include/nr.h"
#define sigma (5.67E-5)
#define NMAX (500)
#define sqr(x) ((x)*(x))

int main(int argc, char *argv[])
{

void legs(double **,int ,double *,double *,int *);
void elgs (double ** ,int ,int* );
void inv(double **, int);
double **amatrix,**bb,*x,*b;
int *indx;
int i,j,bins2=3;

indx=(int *) calloc(bins2,sizeof(int)) ;
b=(double *) calloc(bins2,sizeof(double)) ;
x=(double *) calloc(bins2,sizeof(double)) ;

        for(i=0;i<bins2;i++)
        {
                b[i] = 1.0;
        }

amatrix = malloc(bins2 * sizeof(double *));
bb = malloc(bins2 * sizeof(double *));

        for (i = 0; i < bins2; i++)
                {
                         amatrix[i] = malloc(bins2 * sizeof(double *));
			bb[i] =  malloc(bins2 * sizeof(double *)); 
               }



for (i = 0; i < bins2; i++)
        {
        for (j = 0; j < bins2; j++)
                {
                        amatrix[i][j] = 1.0;
               		bb[i][j] =1.0;
		 }	
        }


for(i=0;i<bins2;i++)
        {
		amatrix[i][i] = 2.0; 

        }




legs(amatrix,bins2,b,x,indx); 
//elgs(amatrix,bins2,indx);
//inv(amatrix,bins2);
for(i=0;i<bins2;i++)
        {
        for(j=0;j<bins2;j++)
                {

                        fprintf(stdout,"%lf ",amatrix[i][j]); //print it out to check


                }

                fprintf(stdout,"\n");


        }





return 0;

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
double tmp;


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
