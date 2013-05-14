chisquaredmin=10000;
int theory2;
for(k=0;k<2000;k++)
{

        amplitude = amplitude0 + delta_amplitude*k; //loop over amplitudes

        for(l=0;l<2000;l++)
        {
                slope = slope0 + delta_slope*l; //loop over slopes
                chisquared = 0; //reset value for each combination of slope and amplitude


                for(i=0;i<bins2;i++) //loop over bins
                {
                        theory2 = amplitude*pow(datatheta[i]/datatheta[(int)bins2/2],slope); //theory is a power law
                        intermediate_step[i]= (datawtheta[i]-theory2)/variance[i]; //Same as Deltas      


                }


                for(i=0;i<bins2;i++) //double loop to measure chisquared
                {
                        for(j=0;j<bins2;j++)
                        {
                                 chisquared+=intermediate_step[i]*covariance_matrix[i][j]*intermediate_step[j]; //cumulative chisquared for each combination of slope and amplitude using covariance matrix                 

                        }
                }



                if(chisquared < chisquaredmin)
                {
                        slope_best = slope;
                        amplitude_best = amplitude;
                        chisquaredmin = chisquared;
                }




        }


}


if((slope_best == slope) | (amplitude_best == amplitude) | (slope_best == slope0) | (amplitude_best == amplitude0))
        {
                fprintf(stderr,"Bad Values for Slope or Amplitude, slope = %lf amp=%lf\n",slope_best,amplitude_best);
                return (-1);
        }

fprintf(stderr,"%s>Best Chisquared = %lf for Slope =%lf and Amplitude=%lf\n",data_file,chisquaredmin/dof,slope_best,amplitude_best);

