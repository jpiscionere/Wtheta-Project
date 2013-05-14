#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define sigma (5.67E-5)
#define NMAX (500)
#define MAXLEN (1000)
#define sqr(x) ((x)*(x))


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

