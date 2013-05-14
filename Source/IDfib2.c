
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI (3.141592)
#include <string.h>
#include "/home/piscioja/NRC/include/nr.h"
#define sqr(x) ((x)*(x))

int main(int argc, char *argv[])
{



/*---Declare-variables-------------------------------------*/

  int Nmax, Nfile_h, Nfile, Ntot, paircounter;
  void Printhelp(void) ;
  void vflip(int ,int *);
  void sort2ii(unsigned long , int *, int *);


 //General Counters
 int i,j,k,l,w,ipick;
 
 //Arrays
  double *hodcolumn1, *hodcolumn2, *hodcolumn3;
  double *CZ,*X,*Y,*Z;
  double dTheta;
  double *m_r,distance,M_r;

 
 //Random 
 float ran2(long *) ;
 long seed  ;
 seed = -42 ; 
 double n1, n2, n3;
 int match;
 
 //Inputs
  FILE *fp1;
  char *hod_file;
  int ifib=1,zcorrect=1;
 
 //output array counters
 int gcounteri,gcounterj,gcounterk,gcounterl,gcountermax ;
 
/*---Allocate-arrrays for RA DEC cZ format-------------------------------------*/
  if(argc < 4)
    {
      Printhelp() ;
      return -1 ;
    } 
  

  sscanf(argv[2],"%d",&ifib) ;
  sscanf(argv[3],"%ld",&seed) ;
  sscanf(argv[4],"%d",&zcorrect) ;


    fprintf(stderr,"Allocate Arrays\n") ;
/*---Allocate-arrrays for RA DEC cZ format-------------------------------------*/
  
  Nmax = 800000;
  hodcolumn1=(double *) calloc(Nmax,sizeof(double)) ; //RA
  hodcolumn2=(double *) calloc(Nmax,sizeof(double)) ; //DEC
  hodcolumn3=(double *) calloc(Nmax,sizeof(double)) ; //CZ
  
 
  
 
  CZ=(double *) calloc(Nmax,sizeof(double)) ;
  X=(double *) calloc(Nmax,sizeof(double)) ;
  Y=(double *) calloc(Nmax,sizeof(double)) ;
  Z=(double *) calloc(Nmax,sizeof(double)) ;
  
 fprintf(stderr,"Read in HOD File for RA DC cZ format\n") ;
   
  hod_file = argv[1];
  fp1=fopen(hod_file,"r") ;
  if (fp1 == NULL) 
   		{
   		fprintf(stderr,"Cannot Open %s\n",hod_file) ;
   		return (-1) ;
   		}
  else {
   fprintf(stderr,"Opened %s\n",hod_file) ;
       }
  l=0;
  while(fscanf(fp1,"%lf %lf %lf",&hodcolumn1[l],&hodcolumn2[l],&hodcolumn3[l])!=EOF) //RA DEC CZ M_R(evolved)
       {

		
		l++ ;
       }
       
  Nfile_h = l ;
  
  fclose(fp1) ;
  
  fprintf(stderr,"readinfile> HOD file has %d lines in it\n",Nfile_h) ;
  
  Nfile = Nfile_h;
  
  /*-----Change into Cartsian and Filter out all z < .02 ------------*/ 

  fprintf(stderr,"Changing Coordinates to Cartisian and Creating Master Galaxy Arrays\n") ;
  
  j=0;

  for(i=0;i<Nfile_h;i++)
  	{
  				CZ[i]=hodcolumn3[i];
  	  			X[i] = sin((90-hodcolumn2[i]) * PI/180)*cos(hodcolumn1[i] * PI/180) ;   
				Y[i] = sin((90-hodcolumn2[i]) * PI/180)*sin(hodcolumn1[i] * PI/180) ;
				Z[i] = cos((90-hodcolumn2[i]) * PI/180) ;
  
  	}
  	

 //\\---------Dynamically Allocate A 2-D Array-------\\//
  fprintf(stderr,"Allocating array space\n") ; 
  int **IDfib;     /* this is the array name */
  int size_x; /* this variable will be used for the first  dimension */
  int size_y; /* this variable will be used for the second dimension */

  /* suppose we want an array of int */
  size_x = 1000000;
  size_y = 20;

  /*  allocate storage for an array of pointers */
  IDfib = malloc(size_x * sizeof(int *));

  /* for each pointer, allocate storage for an array of ints */
  for (i = 0; i < size_x; i++) {
    IDfib[i] = malloc(size_y * sizeof(int));
  }

//Initalizing values of first column
for(i=0;i<Nfile;i++)
	{
		IDfib[i][0] = 0;
		
	}	
fprintf(stderr,"IDfib[1][3] = %d\n",IDfib[1][3]);

for(i=0;i<Nfile;i++)
        {
        for(j=1;j<size_y;j++)
                {

                        IDfib[i][j] = -1;
		 
		}
        }

IDfib[1][2]=-1;
IDfib[1][3]=-1;
IDfib[2][6]=-1;
IDfib[2][7]=-1;

paircounter = 0;

fprintf(stderr,"Starting Galaxy Loop...you better hope this works\n");

fprintf(stderr,"IDfib[1][3] = %d\n",IDfib[1][3]);

 //\\---------Galaxy Loop-------\\//

for(k=0;k<Nfile;k++)
	for(l=k+1;l<Nfile;l++)
		{
			dTheta = acos(X[k] * X[l] + Y[k] * Y[l] + Z[k] * Z[l]);


				if(dTheta <= (11.0*PI)/129600.0 )
		//		if(dTheta <= 0.000290888209 )
				
				{	
				




						IDfib[k][0]++;  //Counts the number of pairs for this galaxy particle
						
 						IDfib[l][0]++;
 						
 						gcounterk = IDfib[k][0];
 						
 						IDfib[k][gcounterk] = l ; //Indeces which galaxy it collides with
						
						gcounterl = IDfib[l][0] ;
						
						IDfib[l][gcounterl] = k ;
 						
						paircounter++ ;
						
						if((gcounterk >= 19) | (gcounterl >= 19))
						{
							fprintf(stderr,"Too many pairs! IDfib Array too Small!\n") ;
							return(-1);
							
						}
					
					}
	 	
	 	}

l=0;
while(IDfib[l][0]==0){
	l++;
	}
	
fprintf(stderr,"IDFib[%d][0] = %d\n IDfib[%d][1] = %d \n",l,IDfib[l][0],l,IDfib[l][1]);

match=0;
for(i=1;i<10;i++)
{
	if( IDfib[IDfib[l][1]][i] == l)
	{
		fprintf(stderr,"Match found...you got lucky...this time.\n");
		fprintf(stderr,"IDfib[%d][%d] = %d\n",IDfib[l][1],i,IDfib[IDfib[l][1]][i]);
		match = 1;
	 }
  }
  
if(match == 0) {
	fprintf(stderr,"Match NOT FOUND.  Fix this ish. \n");
	return(-1);
	}
	
fprintf(stderr,"Done with Galaxy loop!\n") ;
fprintf(stderr,"There are %d Pairs\n",paircounter) ;


 
 j = 0 ;
 i = 0 ;
 gcounteri = 0;
 gcounterj = 0;
 


fprintf(stderr,"Writing Output File\n");

Ntot = Nfile;


//The input for the grouping code needs a null value to be -1

for(i=0;i<Ntot;i++)
        {
        for(j=1;j<10;j++)
                {
                        if(IDfib[i][j] == 0)
                                {
                                IDfib[i][j] = -1 ;
                                }
                }
   }


gcountermax = 0;
int gcounter2=0,*numberofcollisions;
numberofcollisions=(int *) calloc(Ntot,sizeof(int)) ;

for (i=0; i<Ntot; i++)
	{	
		numberofcollisions[i]=IDfib[i][0];
		
		if((IDfib[i][0] > gcountermax))
			{ 	
				gcountermax = IDfib[i][0];
			}
		if(IDfib[i][0] > 1){
			gcounter2++;
			}
		
	}
	
if(ifib == 0 ){
for(i=0;i<Nfile_h;i++)
	{
	for(j=0;j<10;j++)
		{
			fprintf(stdout,"%d  ",IDfib[i][j]) ;
		}
	
	fprintf(stdout,"\n");
   }
}   
fprintf(stderr,"The max number of pairs is %d, the number > 1 = %d\n",gcountermax,gcounter2);

//This part searches through the 2D array and groups together all the connected galaxies - i.e., finds collision groups


int Ngal3,indx,jndx,indx2,Ni,done,Nfibmax,index_counter_counter=0;
int *fiblist,*idrank,*fibcollisions,*index_counter;

Ngal3 = Ntot;
Nfibmax = 50;


fprintf(stderr,"Building fiblist\n");
/*---IDfib[][]-is-ready.--Now-build-fiblist[]----------------*/

    fiblist=(int *) calloc(Nfibmax,sizeof(int));
    idrank=(int *) calloc(Nfibmax,sizeof(int)) ;   
    fibcollisions=(int *) calloc(Ngal3,sizeof(int)) ;
    index_counter=(int *) calloc(Ngal3,sizeof(int));
   fprintf(stderr,"Allocated Memory for fiblist\n"); 
   for(i=0;i<Ngal3;i++)
	{
	  fibcollisions[i] = -1 ;
	}

      for(i=0;i<Ngal3;i++)
        {
          if(IDfib[i][0]>0)
            {
              for(indx=0;indx<Nfibmax;indx++)
                fiblist[indx]=-1 ;

              fiblist[0] = i ;

              indx=jndx=indx2=0 ;
              while(fiblist[indx]!=-1)
                {
                  Ni = IDfib[fiblist[indx]][0] ;
                  for(jndx=1;jndx<Ni+1;jndx++)
                    {
                      done=0 ;
                      for(j=0;j<indx2+1;j++)
                        {
                          if(IDfib[fiblist[indx]][jndx]==fiblist[j]) done=1 ;
                        }
                      if(done==0)
                        {
                          indx2++ ;
                          fiblist[indx2] = IDfib[fiblist[indx]][jndx] ;
                        }
                    }
                  indx++ ;
                }

              for(j=0;j<indx;j++){
                idrank[j] = IDfib[fiblist[j]][0] ;
                }

              sort2ii(indx,(idrank-1),(fiblist-1)) ;
              vflip(indx,fiblist) ;
             	index_counter[index_counter_counter] = indx;
		index_counter_counter++; 

//The IDs of the galaxies in a collision group are stored in the array fiblist[i], and the array is sorted by the number of collisions each galaxy has.
//The next part decides which galaxies to throw out

/*---Eliminate-galaxies-starting-with-ones-with-most-links---*/

/*   @-@  */
              if(indx==2)
                {
                  if(ran2(&seed)<=0.5)
                    fibcollisions[fiblist[0]] = fiblist[1] ;
                  else
                    fibcollisions[fiblist[1]] = fiblist[0] ;
                }
              if(indx==3)
                {

/*   @-@
     |
     @    */
                  if(IDfib[fiblist[1]][0]==1)
                    {
                      fibcollisions[fiblist[0]] = IDfib[fiblist[0]][1] ;
                    }

/*   @-@
     |/
     @    */
                  else
                    {
                      if(ran2(&seed)<=0.333)
                        {
                          fibcollisions[fiblist[0]] = fiblist[2] ;
                          fibcollisions[fiblist[1]] = fiblist[2] ;
                        }
                      else if(ran2(&seed)<=0.666)
                        {
                          fibcollisions[fiblist[0]] = fiblist[1] ;
                          fibcollisions[fiblist[2]] = fiblist[1] ;
                        }
                      else
                        {
                          fibcollisions[fiblist[1]] = fiblist[0] ;
                          fibcollisions[fiblist[2]] = fiblist[0] ;
                        }
                    } 
                }
              if(indx==4)
                
                {
                  

/*   @-@
     |
     @-@  */
                if(IDfib[fiblist[0]][0]==2 && IDfib[fiblist[1]][0]==2 && IDfib[fiblist[2]][0]==1)
                    {
                      for(j=1;j<3;j++)
                        {
                          if(IDfib[fiblist[0]][j] != fiblist[1])
                            {
                              fibcollisions[fiblist[0]] = IDfib[fiblist[0]][j] ;
                              break ;
                            }
                        }
                      for(j=1;j<3;j++)
                        {
                          if(IDfib[fiblist[1]][j] != fiblist[0])
                            {
                              fibcollisions[fiblist[1]] = IDfib[fiblist[1]][j] ;
                              break ;
                            }
                        }
                    }

/*   @-@
     | |
     @-@  */
                 if(IDfib[fiblist[0]][0]==2 && IDfib[fiblist[1]][0]==2 && IDfib[fiblist[2]][0]==2)
                    {
                      for(k=1;k<4;k++)
                        {
                          done=0 ;
                          for(j=1;j<3;j++)
                            {
                              if(fiblist[k]==IDfib[fiblist[0]][j]) done=1 ;
                            }
                          if(done==0)
                            {
                              ipick = k ;
                              break ;
                            }
                        }
                      for(j=1;j<3;j++)
                        {
                          if(IDfib[fiblist[0]][j] != fiblist[ipick])
                            {
                              fibcollisions[fiblist[0]] = IDfib[fiblist[0]][j] ;
                              break ;
                            }
                        }
                      for(j=1;j<3;j++)
                        {
                          if(IDfib[fiblist[ipick]][j] != fiblist[0])
                            {
                              fibcollisions[fiblist[ipick]] = IDfib[fiblist[ipick]][j] ;
                              break ;
                            }
                        }
                    }

/*   @-@
     |\
     @ @  */
                 if(IDfib[fiblist[0]][0]==3 && IDfib[fiblist[1]][0]==1)
                    {
                      fibcollisions[fiblist[0]] = IDfib[fiblist[0]][1] ;
                    }

/*   @-@
     |\
     @-@  */
                  if(IDfib[fiblist[0]][0]==3 && IDfib[fiblist[1]][0]==2)
                    {
                      for(j=1;j<4;j++)
                        {
                          if(IDfib[fiblist[0]][j] != fiblist[1])
                            {
                              fibcollisions[fiblist[0]] = IDfib[fiblist[0]][j] ;
                              break ;
                            }
                        }
                      for(j=1;j<3;j++)
                        {
                          if(IDfib[fiblist[1]][j] != fiblist[0])
                            {
                              fibcollisions[fiblist[1]] = IDfib[fiblist[1]][j] ;
                              break ;
                            }
                        }
                    }

/*   @-@
     |\|
     @-@  */

                  if(IDfib[fiblist[0]][0]==3 && IDfib[fiblist[1]][0]==3 && IDfib[fiblist[2]][0]==2)
                    {
                      for(j=1;j<4;j++)
                        {
                          if(IDfib[fiblist[0]][j] != fiblist[1])
                            {
                              fibcollisions[fiblist[0]] = IDfib[fiblist[0]][j] ;
                              break ;
                            }
                        }
                      for(j=1;j<4;j++)
                        {
                          if(IDfib[fiblist[1]][j] != fiblist[0])
                            {
                              fibcollisions[fiblist[1]] = IDfib[fiblist[1]][j] ;
                              break ;
                            }
                        }
                    }

/*   @-@
     |X|
     @-@  */
                  if(IDfib[fiblist[0]][0]==3 && IDfib[fiblist[1]][0]==3 && IDfib[fiblist[2]][0]==3)
                    {
                      fibcollisions[fiblist[0]] = fiblist[3] ;
                      fibcollisions[fiblist[1]] = fiblist[3] ;
                      fibcollisions[fiblist[2]] = fiblist[3] ;
                    }
                        
               }

              if(indx >=5)
               {

		for(j=0;j<indx;j++)
               	 {
               	   fibcollisions[fiblist[j]] = fiblist[indx - 1] ;
               	 }
               	}
               		

              for(j=0;j<indx;j++)
                {
                  IDfib[fiblist[j]][0] = 0 ;
                }
            }
        }
    

//This array fibcollisions[i] now has an element for each galaxy that is -1 if the galaxy is kept and equal to the ID of the galaxy it collided with if it is thrown out.

/*---Print-out-results------------------------------------------------------*/
fprintf(stderr,"Writing Results File\n");
  j=k=w=0;
  int idmin;
  double thetamax;


for(w=0;w<10;w++)
	fprintf(stderr,"IDfib[1][%d] = %d\n",w,IDfib[1][w]);



 if(ifib==1)                                  /*---ifib=1-------------*/
    {
      for(i=0;i<Nfile_h;i++)
        {
          if(fibcollisions[i] == -1) 
            {

              fprintf(stdout,"%10.6f %10.6f  %8.2f\n",hodcolumn1[i],hodcolumn2[i],hodcolumn3[i]) ;
            }else{
		j++;
		idmin = fibcollisions[i];
		dTheta=acos(X[i] * X[idmin] + Y[i] * Y[idmin] + Z[i] * Z[idmin]);	
		if(zcorrect==1)
			fprintf(stdout,"%10.6f %10.6f  %8.2f\n",hodcolumn1[i],hodcolumn2[i],hodcolumn3[idmin]);

		}	
	
	}

        
        
      fprintf(stderr,"IDfib> Number of galaxies lost to collisions= %d, fraction= %f\n",j,((float)j/(float)Nfile_h));
    }
    

free(IDfib);
free(X);
free(Y);
free(Z);
free(hodcolumn1);
free(hodcolumn2);
free(hodcolumn3);


return(0);
}

void vflip(int n,int *x)
{
   int i ;
   int *tmp ;

   tmp=(int *) calloc(n,sizeof(int)) ;

   for(i=0;i<n;i++)
     *(tmp+i) = *(x+n-1-i) ;

   for(i=0;i<n;i++)
     *(x+i) = *(tmp+i) ;
   return ;
}

void Printhelp(void)
{
  fprintf(stderr,"===========================================================================\n") ;
  fprintf(stderr,"PROGRAM IDfib hodinput ifib seed zcorrect >output.dat\n") ;

}
