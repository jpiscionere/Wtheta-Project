#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "chi2_fof.h"
#include "chi2_fof.c"

#define MAXLEN (100000)

int main(){


        double chi2_fof(int lum_sample,double box_choose,double logM1,double alpha,double gamma,double fgal);

        double chisquared;


        chisquared=chi2(19,1,13.2711595479,1.17625300352,1.82346906053,0.210540118623);
        fprintf(stderr,"chisquared = %lf\n",chisquared);


return 0;

}


