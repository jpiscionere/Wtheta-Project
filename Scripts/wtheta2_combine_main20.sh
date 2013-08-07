#! /bin/sh

wait=$1
file=$2
#------------------
datad=/data2/jap/GalaxyCatalogues/Overlap_Region_Test/Archive
figd=./
#------------------
labx=5
laby=9.5
dlaby=0.6
x1=-3.2
x2=1.2
y1=-2.1
y2=-2.5
y2=2.2
#--------------------------------------------------------------------
#--------------------------------------------------------------------
sm -s <<FIN
if ($wait) {xterm_s} else {laser_sp $figd/$file\n lweight 2}

CTYPE = CTYPE() concat 125 + 256*(125 + 256*125)
CTYPE = CTYPE(STRING) concat 'grey'

CTYPE = CTYPE() concat 50 + 256*(205 + 256*50)
CTYPE = CTYPE(STRING) concat 'goodgreen'

CTYPE = CTYPE() concat 128 + 256*(0 + 256*128)
CTYPE = CTYPE(STRING) concat 'newmagenta'


CTYPE = CTYPE() concat 199 + 256*(97 + 256*20)
CTYPE = CTYPE(STRING) concat 'orange' 

data Wtheta_sdssmock_gamma_main20.rand_10x_sphere.RR.wtheta
read {RR 3}

set dimen(wtheta)=0

macro plot 8 {

set wtheta =0

#-----------------------------------------
#	data \$2.wtheta
#	read {theta 2 wtheta 3 error 4}
 
	do i=3001,3010,1 {
		data Wtheta_stomp_sdssmock_main20_esmeralda.\$i.\$2.full_sphere.wtheta
		read {theta1 2 DD 3}
		set wtheta1 = DD/RR - 1

		set logtheta = lg(theta1)
		set logwtheta_nfw = lg(wtheta1)
		limits $x1 $x2 $y1 $y2
		ctype white
		if(\$1 == 2){
			connect logtheta logwtheta_nfw
			}
		set wtheta = wtheta + wtheta1
		}		
	ctype black
	set wtheta_nfw=wtheta/10
	set logtheta = lg(theta1)
	set logwtheta = lg(wtheta_nfw)
	print<wtheta_nfw>

#----------------------------------Plot-Points
	window 1 1 1 1
	limits $x1 $x2 $y1 $y2

	ctype \$3
        if(\$1==1) {
          set logtheta2=logtheta if(logtheta<\$(lg(55/3600)))
          set logwtheta2=logwtheta+0.523 if(logtheta<\$(lg(50/3600)))
          ptype 20 0 points logtheta2 logwtheta2 
 	  ptype 20 3 points logtheta logwtheta
          set A = (10**0.523)*error if(logtheta<\$(lg(55/3600)))        
          logerr logtheta2 logwtheta2 A

	} else {
        	lweight 2.25
		connect logtheta logwtheta
		lweight 1 
        }

		
#	logerr logtheta logwtheta error 
#	set theory = \$8*theta**\$7 if ( theta < 0.1)
	
#	set lgtheory = lg(theory)
#	ctype grey
#	lweight 4
#	connect logtheta lgtheory
	ctype black
	lweight 1.0
#----------------------------------Labels

#	window 2 2 \$5 \$6
	limits 0 10 0 10
	expand 1.2
        ctype \$3 
	if(\$1==0) {
	   relocate \$($labx+0.5) \$($laby-0.1-\$1*$dlaby) dot
        }
       
	if(\$1==1) {
	  relocate \$($labx+0.5) \$($laby-0.1- \$1*$dlaby) dot
	}

	if(\$1 > 1){

		relocate $labx \$($laby-\$1*$dlaby) draw \$($labx + 1) \$($laby-\$1*$dlaby)

	}

	ltype 0
	ctype default
	expand 1
	relocate \$($labx+1.3) \$($laby-\$1*$dlaby) label \$4 
	expand 1        
}	
        location 4000 30000 4000 30000
        limits  .456 4.85  390 400
        ticksize -1 10 1 5    
        box 4 4 1 4

	expand 1


	
        limits  .456 4.85  390 400
        ticksize -1 10 1 5    
#        box 4 4 1 4

expand 1.0
limits $x1 $x2 $y1 $y2
ticksize -1 10 -1 10 

box 1 2 4 0

limits 0 10 0 10

ctype default
expand 1.2
relocate -1 5 label \omega(\theta)
relocate 5 -0.75 label \theta
expand 1
relocate 4.5 10.5 label R (Kpc)
relocate 5 9 label SDSS Geometry
expand 1


#-------------------------------------------------
#----------------------------------Box
	


#-----------------------------------------
cd $datad

















plot 2 s_o red "Single Overlap" 1 l 1 1
plot 3 d_o blue "Double Overlap" 1 # 1 1
plot 4 n_c goodgreen "No Collisions" 1 l 1 1
plot 5 n_o newmagenta "No Overlap" 1 # 1 1 
#plot 6 s_o.quads_test black "Quads Test" 1 # 1 1 







#-----------------------------------------











if ($wait) {!sleep $wait} else {hardcopy}
FIN
#-----------------------------------------





