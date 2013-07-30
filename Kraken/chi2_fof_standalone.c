#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include "chisquared_function.h"
#include "get_Mmin0.h"

#include "get_Mmin0.c"
#include "chisquared_function.c"

#define MAXLEN (1000)



int main()
{
        double chi2_fof(int ,double, double,double ,double ,double ,double ,double );
        int i;
        double chi2;

        for(i=0;i<9;i++){
                chi2=chi2_fof(18, i,0.19, 9.81,12.7053039145 ,  1.14827081243 ,  0.255258314714 , 1.31705259849);
                fprintf(stderr,"%d %lf\n",i,chi2);
        }
return 0;

}

double chi2_fof(int lum_sample,double box_choose,double siglogM,double logM0,double logM1,double alpha,double gamma,double fgal)
{


	int box;
	//get_Mmin parameters
	char bgc_file[MAXLEN],mf_file[MAXLEN],ztag[MAXLEN];	
	double logMmin=11.0,rhogal=0.0;
	//halobias parameters
	char halobias_path[MAXLEN];
	char halo_directory[MAXLEN],halos_file[MAXLEN], sample_prefix[MAXLEN];
	char dpp_char[MAXLEN];
	char fff_output[MAXLEN];
	double delta_vir, Mstar;
	char halobias_command[MAXLEN],HOD_inputs[MAXLEN];
	char sys_command_rm[MAXLEN];//, sys_command_mv[MAXLEN];
	char file_tag[MAXLEN];
	//makemock parameter
	char makemock_path[MAXLEN];
	char map_file[MAXLEN],map_directory[MAXLEN],galaxy_file[MAXLEN];
	double max_redshift=0.0;
	char makemock_command[MAXLEN];

	//astro_stomp parameters
	char stomp_path[MAXLEN],output_tag[MAXLEN];
	double maxtheta;
	char stomp_command[MAXLEN];
	int binning;

	//fitting parameters
	char	data_directory[MAXLEN],rand_directory[MAXLEN],model_file[MAXLEN],data_file[MAXLEN], randoms_file[MAXLEN], covariance_file[MAXLEN];
	int 	bins2=20,diagonal_errors=0,itrash;
	FILE	*fp1, *fp2, *fp3, *fp4;

	double  *datatheta,*datawtheta,*theory_wtheta;
	double  *error,*RR,trash,chisquared;
	int     i,j,l,N_file,count;
	/* double  chisquared_function (int , double**, double*, double*, double*, double); */
	//double get_Mmin0(double, double , double, double, double, char*, char*); 
	double  **covariance_matrix=NULL;

	const double BAD_VALUE=DBL_MAX;
	int returnvalue;
	
	//Changeables
	delta_vir=377.0;
	Mstar= 2.29E12;
	maxtheta=0.1;
	binning=10;




	
	//halobias chars	
	snprintf(halobias_path,MAXLEN,"/lustre/scratch/jenpi/Source/halobias_fof_nfw_ff");
	
	snprintf(file_tag,MAXLEN,"%10.9f.%10.9f.%10.9f",siglogM,logM0,logM1);
//	snprintf(file_tag,MAXLEN,"tag");
	if(sizeof(file_tag) > MAXLEN)
	  return BAD_VALUE;
	//makemock chars
	snprintf(makemock_path,MAXLEN,"/lustre/scratch/jenpi/Source/makemock_double");
	snprintf(map_directory,MAXLEN,"/lustre/scratch/jenpi/Files");
//	snprintf(map_file,MAXLEN,"%s/lss_geometry.dr72.stripe_trim.pix",map_directory);
	snprintf(map_file,MAXLEN,"%s/lss_geometry.fullsphere_single_index.pix",map_directory);


	//covar chars
	snprintf(data_directory,MAXLEN,"/lustre/scratch/jenpi/Files");	
	snprintf(data_file,MAXLEN,"%s/Wtheta_vollim_Mr%d_fib0.100rand.overlap.filter",data_directory,lum_sample);
        snprintf(covariance_file,MAXLEN,"%s/CovarWtheta_vollim_Mr%d_fib0.100rand.overlap.filter.inv",data_directory,lum_sample);
	snprintf(randoms_file,MAXLEN,"%s/Wtheta_sdssmock_gamma_main%d.rand_100x_sphere.filter.wtheta",data_directory,lum_sample);


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
    	
		
       }



	snprintf(mf_file,MAXLEN,"%s/%s_Main%d_mass_function.out",data_directory,sample_prefix,lum_sample);
	snprintf(stomp_path,MAXLEN,"/lustre/scratch/jenpi/Source/astro-stomp-read-only/examples/.libs/stomp_galaxy_autocorrelation");
	snprintf(halo_directory,MAXLEN,"/lustre/scratch/jenpi/halo_files/%s/Fof",sample_prefix);
	snprintf(bgc_file,MAXLEN,"/lustre/scratch/jenpi/halo_files/%s/Fof/%s_%d_z%s_fof_b0p2.0000.bgc",sample_prefix,sample_prefix,box,ztag);
        snprintf(fff_output,MAXLEN,"%s/fff/halobias_fof_nfw_%s_fff",halo_directory,file_tag);	
	snprintf(galaxy_file,MAXLEN,"%s/fff/halobis_fof_nfw_%s_galaxies",halo_directory,file_tag); 
				



       	snprintf(halos_file,MAXLEN,"%s/%s_%d_z%s_fof_b0p2.%s.halos.ff",halo_directory,sample_prefix,box,ztag,dpp_char);		

		
	logMmin=get_Mmin0(rhogal,siglogM,logM0,logM1,alpha,bgc_file,mf_file);

		 
	snprintf(output_tag,MAXLEN,"stomp_sdssmock_main%d_%s.%s.f%5.4f.g%5.4f.fof.lss_geometry",lum_sample,sample_prefix,file_tag,fgal,gamma);
	
	snprintf(HOD_inputs,MAXLEN,"%lf %lf %lf %lf %lf 1 %lf %lf",logMmin,siglogM,logM0,logM1,alpha,gamma,fgal); 
	snprintf(halobias_command,MAXLEN,"%s 3 4 1 %s %lf %lf %lf -%i < %s > %s",halobias_path,HOD_inputs,-1.0*lum_sample,delta_vir,Mstar,box,halos_file,fff_output);
	



	
	fprintf(stderr,"Executing system command for halobias `%s' \n",halobias_command);
	returnvalue = system(halobias_command);
	//Check exit status
	if(returnvalue != 0)
	  return BAD_VALUE;
	



	snprintf(makemock_command,MAXLEN,"%s 0 0 0 0 0 0 0.02 %lf %s 0.6 0 < %s > %s",makemock_path,max_redshift,map_file,fff_output,galaxy_file);
//	fprintf(stderr,"Executing system command  for makemock `%s' \n",makemock_command);
	returnvalue = system(makemock_command);	
	if (returnvalue!=0)
	  return BAD_VALUE;

	snprintf(sys_command_rm,MAXLEN,"rm -r %s",fff_output);
//	fprintf(stderr,"Executing system command for rm '%s' \n",sys_command_rm);
	system(sys_command_rm);

	snprintf(stomp_command,MAXLEN,"%s --map_file=%s --galaxy_file=%s -output_tag=%s --theta_max=%lf --theta_min=0.002 --n_bins_per_decade=%d -single_index",stomp_path,map_file,galaxy_file,output_tag,maxtheta,binning);
	fprintf(stderr,"Executing system command '%s' \n",stomp_command);
	returnvalue = system(stomp_command);
	if(returnvalue != 0)
	  return BAD_VALUE;
	


	snprintf(sys_command_rm,MAXLEN,"rm -r %s",galaxy_file);	
	system(sys_command_rm);	


	//Chisquared Fitting Routine Starts Here
//	fprintf(stderr,"We have arrived at the chisquared routine\n");




        datatheta=(double *) calloc(bins2,sizeof(double)) ;
        error=(double *) calloc(bins2,sizeof(double)) ;
        datawtheta=(double *) calloc(bins2,sizeof(double)) ;
        theory_wtheta=(double *) calloc(bins2,sizeof(double)) ;
        RR=(double *) calloc(bins2,sizeof(double)) ;




	fp1=fopen(data_file,"r") ;
       	if(fp1==NULL)
	{
		fprintf(stderr,"Chi2 > data file is null..exiting %s\n",data_file);
		return -1;
	} 

        l=0;

//	fprintf(stderr,"Reading in data file %s\n",data_file);	
        while(fscanf(fp1,"%lf %lf %lf %lf %lf %lf",&datatheta[l],&datawtheta[l],&error[l],&trash,&trash,&trash)!=EOF)
       		{
//             		fprintf(stderr,"data %e\n",datatheta[l]);
		
			l++ ;
       		}

        N_file = l ;	
	bins2 = l ;
	fclose(fp1);
	
	
	
	fp3=fopen(randoms_file,"r") ;
      	if(fp3==NULL)
	{
		fprintf(stderr,"randoms_file is null...exiting\n");
		return -1;

	}
 
//	fprintf(stderr,"Reading in %s\n",randoms_file);

        l=0;
        while(fscanf(fp3,"%d %lf %lf %lf ",&itrash,&trash,&RR[l],&trash)!=EOF)
                {
//       			fprintf(stderr,"random %e\n",RR[l]); 
			l++ ;
                }
 	
	fclose(fp3);

	
	covariance_matrix = malloc(bins2 * sizeof(double *));
	for (i = 0; i < bins2; i++)
	  {
	    
	    covariance_matrix[i] = malloc(bins2 * sizeof(double *));
	  }
//	fprintf(stderr,"Reading in %s\n",covariance_file);

	fp2=fopen(covariance_file,"r");
	if(fp2==NULL)
	{
		fprintf(stderr,"covariance_file is null...exiting %s\n",covariance_file);
		return -1;
	}


	l=0;

	for(i=0;i<bins2;i++)
                for(j=0;j<bins2;j++)
                        {
                                fscanf(fp2,"%lf ",&covariance_matrix[i][j]);

                       }

	fclose(fp2);
	snprintf(model_file,MAXLEN,"/lustre/scratch/jenpi/Wtheta/Wtheta_%s.wtheta",output_tag);
        fp4=fopen(model_file,"r");

	if(fp4==NULL)
	{
		fprintf(stderr,"model file is null..exiting\n");
		return -1;
	}       
 
	count=0;
//        fprintf(stderr,"Reading in model file %s\n",model_file);
	int flag=0;
	while(fscanf(fp4,"%d %lf %lf %lf",&itrash,&trash,&theory_wtheta[count],&trash)!=EOF)
        {

		theory_wtheta[count] = theory_wtheta[count]/RR[count] - 1.0;	 	
		
		
		count++;
        }


	
        fclose(fp4) ;
      


	if(count!=bins2)
	{
		fprintf(stderr,"bin number mismatch...exiting\n");
		return -1;
	} 
	
	if(flag == 1){
	  for (i = 0; i < bins2; i++)
	    {
	      free(covariance_matrix[i]);
	    }
	  
	  free(covariance_matrix);
	  free(RR);
	  free(theory_wtheta);
	  free(error);
	  free(datawtheta);
	  return BAD_VALUE;
	}

	snprintf(sys_command_rm,MAXLEN,"rm -r %s",model_file);
        system(sys_command_rm);

	chisquared=chisquared_function(bins2,covariance_matrix, theory_wtheta, datawtheta, error, diagonal_errors);







	for (i = 0; i < bins2; i++)
	  {
	    free(covariance_matrix[i]);
	  }
	free(covariance_matrix);

	free(RR);
	free(theory_wtheta);
	free(error);
	free(datawtheta);

	fprintf(stderr,"\nEnd of chisquared loop Chi2= %lf\n",chisquared);	

	return chisquared;

}
