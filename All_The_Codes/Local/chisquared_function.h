#ifndef _CHISQUARED_FUNCTION_H
#define _CHISQUARED_FUNCTION_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "/home/piscioja/NRC/include/nr.h"
#define sigma (5.67E-5)
#define NMAX (500)
#define sqr(x) ((x)*(x))

double chisquared_function (int dimen, double **covar, double *theory, double *data, double *error, double format);
#endif

