#/bin/bash
GSL_INCL=/usr/bin/include
GSL_LIB=/usr/lib64

gcc -O3 -Wall -I${GSL_INCL}  MatrixInverse.c -Wl,-R${GSL_LIB} -L${GSL_LIB} -lgsl -lgslcblas -lm  -o MatrixInverse
gcc -O3 -Wall -I${GSL_INCL} test.c -Wl,-R${GSL_LIB} -L${GSL_LIB} -lgsl -lgslcblas -lm  -o test
gcc -O3 -Wall -I${GSL_INCL} covariance_fitting2.c -Wl,-R${GSL_LIB} -L${GSL_LIB} -lgsl -lgslcblas -lm  -o covariance_fitting2
gcc -O3 -Wall -I${GSL_INCL} covariance_fitting3.c -Wl,-R${GSL_LIB} -L${GSL_LIB} -lgsl -lgslcblas -lm  -o covariance_fitting3
gcc -O3 -Wall -I${GSL_INCL} covariance_fitting4.c -Wl,-R${GSL_LIB} -L${GSL_LIB} -lgsl -lgslcblas -lm  -o covariance_fitting4
gcc -O3 -Wall -I${GSL_INCL} covariance_fitting.c -Wl,-R${GSL_LIB} -L${GSL_LIB} -lgsl -lgslcblas -lm  -o covariance_fitting
gcc -O3 -Wall -I${GSL_INCL} covariance_fitting_sigma.c -Wl,-R${GSL_LIB} -L${GSL_LIB} -lgsl -lgslcblas -lm  -o covariance_fitting_sigma
gcc -O3 -Wall -I${GSL_INCL} covariance_fitting_fgalgamma.c -Wl,-R${GSL_LIB} -L${GSL_LIB} -lgsl -lgslcblas -lm  -o covariance_fitting_fgalgamma
gcc -O3 -Wall -I${GSL_INCL} covariance_fitting_andreas.c -Wl,-R${GSL_LIB} -L${GSL_LIB} -lgsl -lgslcblas -lm -o covariance_fitting_andreas
