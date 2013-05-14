/* 
Using GSL library to inverse matrix
LU decomposition is used, and inverse is calculated from LU decomposition.
Note this can be inaccurate (consult books on numerical linear algebra for details)
*/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define pi (3.1415926535898)


int main(int argc, char **argv){

  int i,j;
  
  int dimension;
  double tmp;

  if ( argc == 2 ){
    sscanf(argv[1],"%d",&dimension);
  }
  else{
    fprintf( stderr, "Usage:\n %s [dimension] < matrix_file > out_file\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  gsl_matrix * C = gsl_matrix_calloc ( dimension, dimension );
  gsl_matrix * C_i = gsl_matrix_calloc ( dimension, dimension );
  
  for ( i=0;i<dimension;i++ )
    for ( j=0;j<dimension;j++ ){
      //gsl_matrix_fscanf(stdin, C);
      fscanf( stdin,"%le",&tmp );
      gsl_matrix_set ( C, i, j, tmp );
    }

  for ( i=0;i<dimension;i++ ){
    for ( j=0;j<dimension;j++ ) fprintf( stdout, "%le\t", gsl_matrix_get( C, i, j ) );
    fprintf(stdout, "\n");
  }

  int s;

  gsl_permutation *p = gsl_permutation_alloc ( dimension );
  
  gsl_linalg_LU_decomp ( C, p, &s );

  gsl_linalg_LU_invert ( C, p, C_i );

  for ( i=0;i<dimension;i++ ){
    for ( j=0;j<dimension;j++ ) fprintf( stdout, "%le\t", gsl_matrix_get( C, i, j ) );
    fprintf(stdout, "\n");
  }

  for ( i=0;i<dimension;i++ ){
    for ( j=0;j<dimension;j++ ) fprintf( stdout, "%le\t", gsl_matrix_get( C_i, i, j ) );
    fprintf(stdout, "\n");
  }
  


  gsl_permutation_free( p );
  gsl_matrix_free( C );
  gsl_matrix_free( C_i );

  return 0;

}
