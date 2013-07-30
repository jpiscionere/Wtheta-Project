double get_Mmin0(double rhogal, double siglogM, double logM0, double logM1, double alpha, char *bgc_file, char *mf_file);
double chisquared_function (int dimen, double **covar, double *theory, double *data, double *error, double format);
double chi2_fof(int lum_sample,int bender,double box_choose,double siglogM,double logM0,double logM1,double alpha,double gamma,double fgal);
PyMODINIT_FUNC init_chi2_so(void);
