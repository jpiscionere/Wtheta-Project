#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#define sigma (5.67E-5)
#define NMAX (500)
#define MAXLEN (1000)
#define sqr(x) ((x)*(x))


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
