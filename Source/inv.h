#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>





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

