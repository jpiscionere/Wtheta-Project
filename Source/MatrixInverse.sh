#/bin/bash
#source ~/.bashrc


#gcc -O3 -Wall `gsl-config --cflags` MatrixInverse.c `gsl-config --libs` -lm  -o MatrixInverse
#gcc -O3 -Wall `gsl-config --cflags` test.c `gsl-config --libs` -lm  -o test
#gcc -O3 -Wall `gsl-config --cflags` covariance_fitting2.c `gsl-config --libs` -lm  -o covariance_fitting2
#gcc -O3 -Wall `gsl-config --cflags` covariance_fitting3.c `gsl-config --libs` -lm  -o covariance_fitting3
#gcc -O3 -Wall `gsl-config --cflags` covariance_fitting4.c `gsl-config --libs` -lm  -o covariance_fitting4
gcc -O3 -Wall `gsl-config --cflags` covariance_fitting.c `gsl-config --libs` -lm  -o covariance_fitting
#gcc -O3 -Wall `gsl-config --cflags` covariance_fitting_sigma.c `gsl-config --libs` -lm  -o covariance_fitting_sigma
#gcc -O3 -Wall `gsl-config --cflags` covariance_fitting_fgalgamma.c `gsl-config --libs` -lm  -o covariance_fitting_fgalgamma
#gcc -O3 -Wall `gsl-config --cflags` covariance_fitting_andreas.c `gsl-config --libs` -lm  -o covariance_fitting_andreas
gcc -O3 -Wall `gsl-config --cflags` chi2.c `gsl-config --libs` -lm  -o chi2


