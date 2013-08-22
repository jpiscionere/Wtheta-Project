#define NMAX (500)
#define MAXLEN (1000)
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "/home/piscioja/NRC/include/nr.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "chisquared_function.h"
#include "inv.h"
#define SQRT(x) (pow(x,0.5))

int main(int argc, char *argv[]){

double chi2(double, double, double, double, double, double, double, int);
double value;
double logMmin,siglogM,logM0,alpha,logM1,fgal,gamma;



logMmin=11.8013669053163;
siglogM=0.0973817154668221;
logM0=12.687852882335;
logM1=13.0400932413139;
alpha=1.05768884410071;
fgal = 0.205; 
gamma = 1.78;


value = chi2(logMmin,siglogM,logM0,alpha,logM1,fgal,gamma,18);

fprintf(stderr,"value=%lf\n",value);



return 0;

}




double chi2(double logMmin, double siglogM, double logM0, double alpha, double logM1, double fgal, double gamma, int lum_sample )
{

	//halobias parameters
	char halobias_path[MAXLEN];
	char halo_directory[MAXLEN],sample[MAXLEN], sample_prefix[MAXLEN];
	char fff_output[MAXLEN];
	double delta_vir, Mstar, redshift;
	double rhogal;
	char halobias_command[MAXLEN],HOD_inputs[MAXLEN];
	char dpp_char[MAXLEN],mf_file[MAXLEN];
	char bgc_file[MAXLEN],file_tag[MAXLEN],halos_file[MAXLEN];
	char sys_command_rm[MAXLEN], sys_command_mv[MAXLEN];
	int bender=0;
	int box=1;
	double box_choose=1.0;	
	//makemock parameter
	char makemock_path[MAXLEN];
	char map_file[MAXLEN],map_directory[MAXLEN],galaxy_file[MAXLEN];
	double max_redshift;
	char makemock_command[MAXLEN];

	//astro_stomp parameters
	char stomp_path[MAXLEN],output_tag[MAXLEN],ztag[MAXLEN];
	double maxtheta;
	char stomp_command[MAXLEN];
	int binning;

	//fitting parameters
	char	data_directory[MAXLEN],rand_directory[MAXLEN],model_file[MAXLEN],data_file[MAXLEN], randoms_file[MAXLEN], covariance_file[MAXLEN];
	int 	models_total,bins2=20,diagonal_errors=0,itrash;
	char	model_command[MAXLEN];
	FILE	*fp1, *fp2, *fp3, *fp4;

	double  *datatheta,*datawtheta,*theory_wtheta;
	double  *error,*covar,*diagonal,*RR,*intermediate_step,trash,chisquared;
	int     i,j,k,l,N_file,count;
	void    inv (double**,int);
	double  chisquared_function (int , double**, double*, double*, double*, double);
        double  **covariance_matrix;



	//Changeables

	delta_vir=200.0;
	Mstar= 2.29E12;
	maxtheta=0.1;
	binning=10;

	


	//halobias chars

	snprintf(halobias_path,MAXLEN,"/home/piscioja/halobias/halobias_fof_nfw");
	snprintf(HOD_inputs,MAXLEN,"%lf %lf %lf %lf %lf 1 %lf %lf",logMmin,siglogM,logM0,logM1,alpha,gamma,fgal); 

	//makemock chars
	snprintf(makemock_path,MAXLEN,"/home/piscioja/halobias/makemock");
	snprintf(map_directory,MAXLEN,"/home/piscioja/SDSSPix/Maps");
	snprintf(map_file,MAXLEN,"%s/lss_geometry.dr72.stripe_trim.pix",map_directory);


	//stomp chars
	snprintf(stomp_path,MAXLEN,"/home/piscioja/astrostomp/astro-stomp-read-only/examples/stomp_galaxy_autocorrelation");

	//covar chars
	snprintf(data_directory,MAXLEN,"/home/piscioja/Clustering/WthetaPaper/Data");	
        snprintf(data_file,MAXLEN,"%s/Wtheta/Wtheta_vollim_Mr%d_fib0.100rand.overlap.filter",data_directory,lum_sample); 
        snprintf(covariance_file,MAXLEN,"%s/Wtheta/CovarWtheta_vollim_Mr%d_fib0.100rand.overlap.filter",data_directory,lum_sample);
        snprintf(rand_directory,MAXLEN,"/hd0/Research/Clustering/Randoms");
	snprintf(randoms_file,MAXLEN,"%s/Wtheta_sdssmock_gamma_main%d.rand_10x_weighted.ns.rdcz.wtheta",rand_directory,lum_sample);


        if(lum_sample==18 || lum_sample==500){
                box= 4020 + (int)floor(box_choose);
                snprintf(sample_prefix,MAXLEN,"Consuelo");
                snprintf(ztag,MAXLEN,"0p054");
                Mstar=2.49E12;
                max_redshift=0.042;
                rhogal=0.03;
                snprintf(dpp_char,MAXLEN,"fdpp");
        }else if(lum_sample==19){
                box= 4020 + (int)floor(box_choose);
                snprintf(sample_prefix,MAXLEN,"Consuelo");
                snprintf(ztag,MAXLEN,"0p054");
                Mstar=2.29E12;
                max_redshift=0.067;
                rhogal=0.0157;
                snprintf(dpp_char,MAXLEN,"fdpp");

        }else if(lum_sample==20){
                box=3011 + (int)floor(box_choose);
                snprintf(sample_prefix,MAXLEN,"Esmeralda");
                snprintf(ztag,MAXLEN,"0p000");
                Mstar=2.49E12;
                max_redshift=0.106;
                rhogal=0.0063;
                snprintf(dpp_char,MAXLEN,"dpp");

        }else if(lum_sample==21){
                box=(int) 2020 + floor(box_choose);
                snprintf(sample_prefix,MAXLEN,"Carmen");
                snprintf(ztag,MAXLEN,"0p132");
                Mstar=2.49E12;
                max_redshift=0.165;
                rhogal=0.0012;
                snprintf(dpp_char,MAXLEN,"dpp");




        }else{
                fprintf(stderr,"Incorrect Lum Sample! Exiting......\n");
                return -1;

        }


        snprintf(mf_file,MAXLEN,"/home/piscioja/Clustering/WthetaPaper/Emcee_Outputs/Run_%s/%s_Main%d_mass_function.out",sample_prefix,sample_prefix,lum_sample);



        if(bender==1){
                snprintf(stomp_path,MAXLEN,"/home/piscioja/astrostomp/stomp_bender/examples/stomp_galaxy_autocorrelation");
                snprintf(halo_directory,MAXLEN,"/data0/LasDamas/%s/%d",sample_prefix,box);
                if(lum_sample != 21)
                        snprintf(bgc_file,MAXLEN,"%s/fof_b0p2/%s_%d_z%s_fof_b0p2.0000.bgc",halo_directory,sample_prefix,box,ztag);
                else   
                        snprintf(bgc_file,MAXLEN,"%s/%s_%d_z%s_fof_b0p2.0000.bgc",halo_directory,sample_prefix,box,ztag);

                snprintf(fff_output,MAXLEN,"/data2/jap/Emcee_test/fff/halobias_fof_nfw_%s_fff",file_tag);
                snprintf(galaxy_file,MAXLEN,"/data2/jap/Emcee_test/fff/halobis_fof_nfw_%s_galaxies",file_tag);



        }else{ 

                snprintf(stomp_path,MAXLEN,"/home/piscioja/Bin/stomp_galaxy_autocorrelation");
                snprintf(halo_directory,MAXLEN,"/ssd1/Research/halo_files/%s/Fof",sample_prefix);
                snprintf(bgc_file,MAXLEN,"/ssd1/Research/halo_files/%s/Fof/%s_%d_z%s_fof_b0p2.0000.bgc",sample_prefix,sample_prefix,box,ztag);
                snprintf(fff_output,MAXLEN,"%s/fff/halobias_fof_nfw_%s_fff",halo_directory,file_tag);
                snprintf(galaxy_file,MAXLEN,"%s/fff/halobis_fof_nfw_%s_galaxies",halo_directory,file_tag);


        }

        snprintf(halos_file,MAXLEN,"%s/%s_%d_z%s_fof_b0p2.%s.halos.ff",halo_directory,sample_prefix,box,ztag,dpp_char);










	

	//Model Wtheta Assembly

/*
	snprintf(halobias_command,MAXLEN,"%s 3 4 1 %s %d %lf %lf %lf trashfile.out -%i < %s > %s",halobias_path,HOD_inputs,-1*lum_sample,delta_vir,Mstar,redshift,box,sample,fff_output);
	fprintf(stderr,"Executing system command '%s' \n",halobias_command);
	system(halobias_command);

	snprintf(makemock_command,MAXLEN,"%s 1 1 0 0 0 0 0.02 %lf %s 0.6 0 < %s > %s",makemock_path,max_redshift,map_file,fff_output,galaxy_file);
	fprintf(stderr,"Executing system command '%s' \n",makemock_command);
	system(makemock_command);	

	snprintf(sys_command_rm,MAXLEN,"rm -r %s",fff_output);
	fprintf(stderr,"Executing system comman '%s' \n",sys_command_rm);
//	system(sys_command_rm);

	snprintf(stomp_command,MAXLEN,"%s --map_file=%s --galaxy_file=%s -output_tag=%s --theta_max=%lf --n_bins_per_decade=%d -single_index",stomp_path,map_file,galaxy_file,output_tag,maxtheta,binning);
	fprintf(stderr,"Executing system command '%s' \n",stomp_command);
	system(stomp_command);
	
	

	snprintf(sys_command_rm,MAXLEN,"rm -r %s",galaxy_file);	
	system(sys_command_rm);	

*/
	//Chisquared Fitting Routine Starts Here
	fprintf(stderr,"We have arrived at the chisquared routine\n");


	covar=(double * ) calloc(bins2*bins2,sizeof(double));
        diagonal=(double * ) calloc(bins2,sizeof(double));
        datatheta=(double *) calloc(bins2,sizeof(double)) ;
        error=(double *) calloc(bins2,sizeof(double)) ;
        datawtheta=(double *) calloc(bins2,sizeof(double)) ;
        intermediate_step=(double *) calloc(bins2,sizeof(double)) ;
        theory_wtheta=(double *) calloc(bins2,sizeof(double)) ;
        RR=(double *) calloc(bins2,sizeof(double)) ;

	
	fprintf(stderr,"Done allocating\n");
 	
	fp1=fopen(data_file,"r") ;
        assert( fp1 != NULL );
        l=0;
	fprintf(stderr,"Reading in data file %s\n",data_file);	
        while(fscanf(fp1,"%lf %lf %lf %lf %lf %lf",&datatheta[l],&datawtheta[l],&error[l],&trash,&trash,&trash)!=EOF)
       		{
        		fprintf(stderr,"%lf\n",datawtheta[l]);     
			l++ ;
       		}

        N_file = l ;
	bins2=N_file;	
	fclose(fp1);
	


	fp3=fopen(randoms_file,"r") ;
        assert( fp3 != NULL );
	fprintf(stderr,"Reading in %s\n",randoms_file);

        l=0;
        while(fscanf(fp3,"%d %lf %lf %lf ",&itrash,&trash,&RR[l],&trash)!=EOF)
                {
        
			l++ ;
                }
 	
	fclose(fp3);

	
	covariance_matrix = malloc(bins2 * sizeof(double *));
	
	for (i = 0; i < bins2; i++)
                {

			covariance_matrix[i] = malloc(bins2 * sizeof(double *));

                 
		}
	fprintf(stderr,"Reading in %s\n",covariance_file);

	fp2=fopen(covariance_file,"r");
	assert(fp2!= NULL);
	l=0;
        while(fscanf(fp2,"%lf %lf %lf",&trash,&trash,&covar[l])!=EOF)
		{
        		fprintf(stderr,"%lf\n",covar[l]);     
                	l++;
		}
	fprintf(stderr,"%d\n",l);

       k=0;

       for(j=0;j<bins2;j++)
       {
       		for(i=0;i<bins2;i++)
		{
			covariance_matrix[j][i]=covar[k];
	
			k++;
		}
	
		diagonal[j]=covariance_matrix[j][j];
		fprintf(stderr,"%lf\n",diagonal[j]); 
			
	}


	for(j=0;j<bins2;j++){
		for(i=0;i<bins2;i++){

			fprintf(stderr,"%lf\n",covariance_matrix[j][i]); 
			}

                }
	fprintf(stderr,"%d\n",bins2);	

//	inv(covariance_matrix,bins2);
	
	 for(i=0;i<bins2;i++){
                for(j=0;j<bins2;j++)
                        {
                                fprintf(stdout,"%lf ",covariance_matrix[i][j]/SQRT(diagonal[j]*diagonal[i]));

                        }
		fprintf(stdout,"\n");
	}
	
/*
	snprintf(model_file,MAXLEN,"Wtheta_%s.wtheta",output_tag);
        fp4=fopen(model_file,"r");
        assert(fp4 != NULL);
	count=0;
        fprintf(stderr,"Reading in model file %s\n",model_file);
	double try;
	while(fscanf(fp4,"%d %lf %lf %d",&itrash,&trash,&theory_wtheta[count],&itrash)!=EOF)
        {
               	
		theory_wtheta[count] = theory_wtheta[count]/RR[count] - 1.0;	 
		count++;
        }

	
        fclose(fp4) ;
        assert(count == bins2);


	chisquared=chisquared_function(bins2,covariance_matrix, theory_wtheta, datawtheta, error, diagonal_errors);
*/
	return 2;

}
