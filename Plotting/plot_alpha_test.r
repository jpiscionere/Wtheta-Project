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

data=read.table("/home/piscioja/Clustering/WthetaPaper/Data/Wtheta/Wtheta_vollim_Mr20_fib0.20rand.overlap.short")
 wtheta_data=data.frame(data)
plot(V2~V1,data=wtheta_data,log="xy",axes=FALSE,type="p",col="red",ylim=c(1,1500),xlab=expression(theta),ylab=expression(omega(theta)),main="Esmeralda Varying LogM1")
magaxis(log='xy')
box()
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="red")
data_files<-list.files(".",pattern="average.fof.logM1",full.names=T)
for(i in 1:length(data_files)){a=read.table(data_files[i]);model=data.frame(a);lines(V2~V1,data=model,type="l",col=colors()[170 + i],lwd=3)}
logM1=seq(12.75,13.75,by=0.1)
for(i in 1:length(logM1)){a=read.table(data_files[i]);model=data.frame(a);text(model$V1[i+4],model$V2[i+4],logM1[i],col=colors()[170 + i*4],pos=3,cex=0.75,srt=-45)}
 data<-read.table("Wtheta_main20_esmeralda.average.fof.logM113.05.wtheta")
model=data.frame(data)
lines(V2~V1,data=model,type="l",col="blue",lwd=3)
names=c("SDSS Mr < -20","CKM FOF Wp Best Fit","Varying logM1")
colors=c("red","blue","dark gray")
legend("topright",names,col=colors,lty=c(1,1,1),pch=c(1,NA,NA))

library("Hmisc")
library('magicaxis')
library("foreign")

data=read.table("/home/piscioja/Clustering/WthetaPaper/Data/Wtheta/Wtheta_vollim_Mr20_fib0.20rand.overlap.short")
 wtheta_data=data.frame(data)
plot(V2~V1,data=wtheta_data,log="xy",axes=FALSE,type="p",col="red",ylim=c(1,50),xlab=expression(theta),ylab=expression(omega(theta)),main="Esmeralda Varying LogM1")
magaxis(log='xy')
box()
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="red")
data_files<-list.files(".",pattern="average.so",full.names=T)
for(i in 1:length(data_files)){a=read.table(data_files[i]);model=data.frame(a);lines(V2~V1,data=model,type="l",col=colors()[170 + i*4],lwd=3)}
logM1=seq(12.75,13.75,by=0.1)
for(i in 1:length(logM1)){a=read.table(data_files[i]);model=data.frame(a);text(model$V1[i],model$V2[i],logM1[i],col=colors()[170 + i*4],pos=4,cex=0.75,srt=0)}
 data<-read.table("Wtheta_main20_esmeralda.average.so.logM113.05.wtheta")
model=data.frame(data)
lines(V2~V1,data=model,type="l",col="blue",lwd=3)
names=c("SDSS Mr < -20","CKM FOF Wp Best Fit","Varying logM1")
colors=c("red","blue","dark gray")
legend("topright",names,col=colors,lty=c(1,1,1),pch=c(1,NA,NA))

data=read.table("fof_pnm_file")
data_hist<-hist(data$V2)
breaks = data_hist$breaks
breaks=append(breaks,15.6)
plot(data_hist$mids,log_counts,col="white",log="xy",main="Friends of Friends vs. Spherical Overdensity Mass Functions",xlab="LogMass",ylab="dLogN/dlogM")

data_files<-list.files(".",pattern="fof_pnm",full.names=T)
for(i in 1:length(data_files)){a=read.table(data_files[i]);data_hist<-hist(a$V2,plot=FALSE,breaks=breaks);log_counts=log10(data_hist$counts);lines(log_counts~data_hist$mids,type="l",col="red",lwd=1,lty=3)}
data_files<-list.files(".",pattern="so_pnm",full.names=T)
for(i in 1:length(data_files)){a=read.table(data_files[i]);data_hist<-hist(a$V2,plot=FALSE,breaks=breaks);log_counts=log10(data_hist$counts);lines(log_counts~data_hist$mids,type="l",col="blue",lwd=1,lty=3)}
names=c("Fof Halos","SO Halos")
colors=c("red","blue")
legend("topright",names,col=colors,lty=c(3,3),pch=c(NA,NA))


data=read.table("xi_so_nfw_3011.out") 
xi=data.frame(r=data$V1,xi=data$V4,err=data$V5)
plot(xi~r,data=xi,log="xy",axes=FALSE,type="p",col="white",xlab="r (Mpc/h)",ylab=expression(xi),main="Esmeralda Fof vs. So",ylim=c(10^-1.5,10^3))
magaxis(log="xy")
box()
#errbar(xi$r,xi$xi,xi$xi+xi$err,xi$xi-xi$err,log="xy",add=TRUE,col="red")
data_files<-list.files(".",pattern="xi_so",full.names=T)
average=0*c(1:length(data$V4))
for(i in 1:length(data_files)){a=read.table(data_files[i]);model=data.frame(r=a$V1,xi=a$V4,err=a$V5);average=average+model$xi;lines(xi~r,data=model,type="l",col=colors()[170 + i*4],lwd=3)}
average=average/length(data_files)
lines(xi$r,average,col="red")
average=0
data_files<-list.files(".",pattern="xi_fof",full.names=T)
for(i in 1:length(data_files)){a=read.table(data_files[i]);model=data.frame(r=a$V1,xi=a$V4,err=a$V5);average=average+model$xi;lines(xi~r,data=model,type="l",col=colors()[170 + i*4],lwd=3)}
average=average/length(data_files)
lines(xi$r,average,col="blue")
names=c("So Halos","Fof Halos","xi=(r/5.0)^-1.8")
colors=c("red","blue","black")
legend("topright",names,col=colors,lty=c(1,1,3))
xi_model=(xi$r/5)^-1.8
lines(xi$r,xi_model,lty=3)


#plotting walkers
data=read.table("../Box_test/tmp2",head=T)
 walker=rep(1:500,length(data$Prob)/500)
RedChi2=2*abs(data$Prob)/14
fit<-data.frame(walker=walker,RedChi2,time=data$Time)
p<-ggplot(fit,aes(x=time,y=RedChi2,group=walker)) + ylim(1.5,10)
 p + geom_line(aes(colour=walker)) + scale_colour_gradient(low="red") + xlab("Step") + ylab(Chi^2~"/DOF")

data=read.table("/home/piscioja/Clustering/WthetaPaper/Data/Wtheta/Wtheta_vollim_Mr20_fib0.20rand.overlap.short")
 wtheta_data=data.frame(data)
plot(V2~V1,data=wtheta_data,log="xy",axes=FALSE,type="p",col="red",ylim=c(1,100),xlab=expression(theta),ylab=expression(omega(theta)),main="Esmeralda Best Fit HOD")
magaxis(log='xy')
box()
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="red")
data_files<-list.files(".",pattern="best_fit2.wtheta",full.names=T)
randoms<-read.table("/hd0/Research/Clustering/Randoms/Wtheta_sdssmock_gamma_main20.rand_100x_sphere.wtheta")
for(i in 1:length(data_files)){a=read.table(data_files[i]);model=data.frame(a);wtheta=model$V3/randoms$V3 - 1;lines(model$V2,wtheta,type="l",col=colors()[170 + i*4],lwd=3)}
#data_files<-list.files(".",pattern=".nfw.fof.lss_geometry2.wtheta",full.names=T)
#for(i in 1:length(data_files)){a=read.table(data_files[i]);model=data.frame(a);lines(V3~V2,data=model,type="l",col=colors()[170 + i*4],lwd=3)}


data=read.table("/home/piscioja/Clustering/WthetaPaper/Data/Wtheta/Wtheta_vollim_Mr19_fib0.100rand.overlap.filter")
 wtheta_data=data.frame(data)
plot(V2~V1,data=wtheta_data,log="xy",axes=FALSE,type="p",col="red",xlab=expression(theta),ylab=expression(omega(theta)),main="Consuelo Best Fit Fof HOD")
magaxis(log='xy')
box()
errbar(wtheta_data$V1,wtheta_data$V2,wtheta_data$V2+wtheta_data$V3,wtheta_data$V2-wtheta_data$V3,log="xy",add=TRUE,col="red")
data_files<-list.files("/net/bender/data2/jap/Emcee_test/Consuelo_Run",pattern="main19_consuelo",full.names=T)
randoms<-read.table("/home/piscioja/Clustering/WthetaPaper/Data/Randoms/Wtheta_sdssmock_gamma_main19.rand_100x_sphere.filter.wtheta")
for(i in 1:length(data_files)){a=read.table(data_files[i]);model=data.frame(a);wtheta=model$V3/randoms$V3 - 1;lines(model$V2,wtheta,type="l",col=colors()[170 + i*4],lwd=3)}
#for(i in 1:length(data_files)){a=read.table(data_files[i]);model=data.frame(a);wtheta=model$V3;lines(model$V2,wtheta,type="l",col=colors()[170 + i*4],lwd=3)}

#Convergence Test
Iteration=floor(data$Time/500)
stdev_measure=c(1:length(unique(Iteration)))
for(i in 1:length(stdev_measure)){stdev_measure[i]=sd(data$Fgal[which(Iteration==i)])}
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, colour="red", geom=geom, width=0.2, ...)
}


fgal<-ggplot(fit,aes(Iteration,Fgal)) + stat_smooth()+ stat_sum_df("mean_cl_normal",geom="smooth") 
gamma<-ggplot(fit,aes(Iteration,Gamma)) + stat_smooth(level=0.999)+ stat_sum_df("mean_cl_normal",geom="smooth",color="blue")
chi2<-ggplot(fit,aes(Iteration,RedChi2)) + stat_sum_df("mean_cl_normal",geom="smooth")
siglogm<-ggplot(fit,aes(Iteration,SigLogM)) + stat_smooth(level=0.999) + stat_sum_df("mean_cl_normal",geom="smooth")
 grid.arrange(fgal,gamma,chi2,siglogm)


siglogm<-ggplot(fit,aes(Iteration,SigLogM)) + stat_sum_df("mean_cl_normal",geom="smooth")
logM0<-ggplot(fit,aes(Iteration,LogM0)) + stat_sum_df("mean_cl_normal",geom="smooth")
logM1<-ggplot(fit,aes(Iteration,LogM1)) + stat_sum_df("mean_cl_normal",geom="smooth")
alpha<-ggplot(fit,aes(Iteration,Alpha)) + stat_sum_df("mean_cl_normal",geom="smooth")
gamma<-ggplot(fit,aes(Iteration,Gamma)) + stat_sum_df("mean_cl_normal",geom="smooth")
fgal<-ggplot(fit,aes(Iteration,Fgal)) + stat_sum_df("mean_cl_normal",geom="smooth")

grid.arrange(fgal,gamma,chi2,siglogm,logM0,logM1,alpha)



myFunc = function(x) { 
	result=c(mean(x) - sd(x),mean(x)+sd(x))
	names(result)=c("ymin","ymax") 
	result
	}


#possible ellispes
 qplot(data=data,x=Fgal,y=Gamma) + stat_ellipse(geom="polygon",alpha=1,level=0.68,color="red",aes(fill="red")) + stat_ellipse(geom="polygon",alpha=0.51,level=0.95,color="pink",aes(fill="pink")) 
ci2d(data$Fgal,data$Gamma)
ggpairs(wtheta,columns=c(3:8), lower=list(continuous="density"), diag=list(continuous="density"),axisLabels="show")
