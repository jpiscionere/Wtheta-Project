library("MASS")
library("ggplot2")
library("gridExtra")
library("Hmisc")
library('magicaxis')
library("foreign")
library("gplots")
source("../Box_test/stat-ellipse.R")
library('gridExtra')
library(GGally)


#par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1),mgp = c(1.5, 1, 0) )
par(mfrow=c(2,2), mgp = c(1.5, 1, 0) )



data=read.table("/home/piscioja/Clustering/WthetaPaper/Data/Wtheta/Wtheta_vollim_Mr18_fib0.20rand.overlap.short")
 wtheta_data=data.frame(data)
plot(V2~V1,data=wtheta_data,log="xy",axes=FALSE,type="p",col="deeppink",ylim=c(1,100),xlab=expression(theta),ylab=expression(omega(theta)))

magaxis(log='xy')
box()
#par(fg = "deeppink") 
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="deeppink")
par(fg ="black")
names=c("SDSS Mr < -18","Slope =-0.634",expression(paste(  chi^2/Dof==0.87)))
colors=c("deeppink","dark gray",NA)
legend("topright",names,col=colors,lty=c(1,1,NA),pch=c(1,NA,NA),bty="n")
data=read.table("/hd0/Research/Clustering/Emcee_test/Power_Law/18_power_law")
for(i in c(1:100)) {theory=data$V1[i]*(wtheta_data$V1/wtheta_data$V1[length(wtheta_data$V1)/2 + 1])^data$V2[i];lines(wtheta_data$V1,theory,col=rgb(150,150,150,50,maxColorValue=255))}
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="deeppink")


data=read.table("/home/piscioja/Clustering/WthetaPaper/Data/Wtheta/Wtheta_vollim_Mr19_fib0.20rand.overlap.short")
wtheta_data=data.frame(data)
plot(V2~V1,data=wtheta_data,log="xy",axes=FALSE,type="p",col="deeppink",ylim=c(1,100),xlab=expression(theta),ylab=expression(omega(theta)))
magaxis(log='xy')
box()
#par(fg = "deeppink") 
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="deeppink")
par(fg = "black") 
names=c("SDSS Mr < -19  ","Slope = -0.7007",expression(paste(  chi^2/Dof==1.06)))
colors=c("deeppink","dark gray",NA)
legend("topright",names,col=colors,lty=c(1,1,NA),pch=c(1,NA,NA),bty="n")
data=read.table("/hd0/Research/Clustering/Emcee_test/Power_Law/19_power_law")
for(i in c(1:100)) {theory=data$V1[i]*(wtheta_data$V1/wtheta_data$V1[length(wtheta_data$V1)/2 + 1])^data$V2[i];lines(wtheta_data$V1,theory,col=rgb(150,150,150,50,maxColorValue=255))}
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="deeppink")



data=read.table("/home/piscioja/Clustering/WthetaPaper/Data/Wtheta/Wtheta_vollim_Mr20_fib0.20rand.overlap.short")
 wtheta_data=data.frame(data)
plot(V2~V1,data=wtheta_data,log="xy",axes=FALSE,type="p",col="deeppink",ylim=c(1,100),xlab=expression(theta),ylab=expression(omega(theta)))
magaxis(log='xy')
box()
#par(fg = "deeppink") 
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="deeppink")
par(fg = "black") 
names=c("SDSS Mr < -20  ","Slope = -0.757",expression(paste(  chi^2/Dof==0.71)))
colors=c("deeppink","dark gray",NA)
legend("topright",names,col=colors,lty=c(1,1,NA),pch=c(1,NA,NA),bty="n")
data=read.table("/hd0/Research/Clustering/Emcee_test/Power_Law/20_power_law")
for(i in c(1:100)) {theory=data$V1[i]*(wtheta_data$V1/wtheta_data$V1[length(wtheta_data$V1)/2 + 1])^data$V2[i];lines(wtheta_data$V1,theory,col=rgb(150,150,150,50,maxColorValue=255))}
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="deeppink")




data=read.table("/home/piscioja/Clustering/WthetaPaper/Data/Wtheta/Wtheta_vollim_Mr21_fib0.20rand.overlap.short")
 wtheta_data=data.frame(data)
plot(V2~V1,data=wtheta_data,log="xy",axes=FALSE,type="p",col="deeppink",ylim=c(1,100),xlab=expression(theta),ylab=expression(omega(theta)))
magaxis(log='xy')
box()
#par(fg = "deeppink") 
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="deeppink")
par(fg = "black") 
names=c("SDSS Mr < -21","Slope =-0.933", expression(paste(  chi^2/Dof==0.813)))
colors=c("deeppink","dark gray",NA)
legend("topright",names,col=colors,lty=c(1,1,NA),pch=c(1,NA,NA),bty="n")
data=read.table("/hd0/Research/Clustering/Emcee_test/Power_Law/21_power_law")
for(i in c(1:100)) {theory=data$V1[i]*(wtheta_data$V1/wtheta_data$V1[length(wtheta_data$V1)/2 + 1] )^data$V2[i];lines(wtheta_data$V1,theory,col=rgb(150,150,150,50,maxColorValue=255))}
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="deeppink")


data=read.table("/home/piscioja/Clustering/WthetaPaper/Data/Wtheta/Wtheta_vollim_Mr20_fib0.20rand.overlap.short")
 wtheta_data=data.frame(data)
plot(V2~V1,data=wtheta_data,log="xy",axes=FALSE,type="p",col="deeppink",ylim=c(1,100),xlab=expression(theta),ylab=expression(omega(theta)))
magaxis(log='xy')
box()

errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="deeppink")

names=c("SDSS Mr < -20  ",expression(paste(gamma == 1.8598)), expression(paste(f[gal]==0.15488)),expression(paste(chi^2/Dof==2.03)))
colors=c("deeppink","dark gray",NA,NA)
legend("topright",names,col=colors,lty=c(1,1,NA,NA),pch=c(1,NA,NA,NA),bty="n")
data_files<-list.files("/hd0/Research/Clustering/Emcee_test/Consuelo_Run/",pattern="0.036258163.12.002614185.13.301923187.f0.2463.g1.8598.",full.names=T)
randoms<-read.table("/hd0/Research/Clustering/Randoms/Wtheta_sdssmock_gamma_main20.rand_100x_sphere.wtheta")
for(i in 1:length(data_files)){a=read.table(data_files[i]);model=data.frame(a);wtheta=model$V3/randoms$V3 - 1;lines(model$V2,wtheta,type="l",col=colors()[170 + i*4],lwd=3)}
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="deeppink")

bp<-qplot(Gamma, data=wtheta, geom="density", fill=LumSample, alpha=I(.5),xlab=expression(paste(gamma)),ylab="Density")
bp + scale_fill_discrete(name="Luminosity Sample",breaks=c("Main18","Main19","Main20"),labels=c(expression(paste(M[r] < -18)),expression(paste(M[r] < -19)),expression(paste(M[r] < -20))))


data=read.table("/home/piscioja/Clustering/WthetaPaper/Data/Wtheta/Wtheta_vollim_Mr20_fib0.20rand.overlap.short")
 wtheta_data=data.frame(data)
plot(V2~V1,data=wtheta_data,log="xy",axes=FALSE,type="p",col="deeppink",ylim=c(1,100),xlab=expression(theta),ylab=expression(omega(theta)))
top=tan(wtheta_data$V1)*0.086*10^3
magaxis(log='xy')
box()

errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="deeppink")
randoms<-read.table("~/Clustering/WthetaPaper/NFW/Wtheta_sdssmock_gamma_main20.rand_10x_weighted.overlap.rdcz.wtheta")
data_files<-list.files("~/Clustering/WthetaPaper/NFW",pattern="overlap",full.names=T)
wtheta_mean=wtheta*0
for(i in 1:length(data_files)){a=read.table(data_files[i]);model=data.frame(a);wtheta=model$V3/randoms$V3 - 1;wtheta_mean=wtheta_mean + wtheta;lines(wtheta_data$V1,wtheta,col=rgb(150,150,150,50,maxColorValue=255))}
wtheta_mean = wtheta_mean/length(data_files)
lines(wtheta_data$V1,wtheta_mean,col="springgreen3",lwd="3")
y=c(0.00001: 10000)	
x =y/y *0.014
lines(x,y,col="cyan",lty=4)


