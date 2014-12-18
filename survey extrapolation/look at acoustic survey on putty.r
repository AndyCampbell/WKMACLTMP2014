
rm(list=ls())


setwd("/media/m/old stuff/chile project/simulator/V12 herring again extended survey/survey extrapolation")

load("NL_data_5nm.RData")
load("NO_data_5nm.RData")
load("UK_data_5nm.RData")

NL_data_5nm$country<-"NL"
NO_data_5nm$country<-"NO"
UK_data_5nm$country<-"UK"


#NO_data_5nm[,1]<-NO_data_5nm[,1]/1000


dat<-rbind(NL_data_5nm,NO_data_5nm,UK_data_5nm)
names (dat) <- c("Long_M","Lat_M","Nasc","Year","country")
dat$Nasc[is.na(dat$Nasc)]<-0
dat<-aggregate(dat$Nasc,by=list(lat=dat$Lat_M,long=dat$Long_M,year=dat$Year,country=dat$country),FUN=sum)
names(dat)[5]<-"abund"
 library(lattice)
xyplot(lat~long|year,groups=country,data=dat)


range(dat$lat)    #53.653 61.940

range(dat$long)   # -4.307167  7.977372




Ymin<-0
Ymax<-(62-53)*60

Xmin<-0
Xmax<-round((8+4.5)*60*cos(53*2*pi/360))

dat$x<-dat$long
dat$y<-dat$long

dat$y<-round((dat$lat-53)*60,0)
dat$x<-round((dat$long+4.5)*60*cos(dat$lat*2*pi/360),0)

par(mfrow=c(1,2))
plot(dat$long,dat$lat)
plot(dat$x,dat$y)


# select one year
yr<-2005

dat<-dat[dat$year==yr,]


symbols(dat$x,dat$y,circles=dat$abund)
points(dat$x,dat$y,col="red",cex=0.2)


ppp<-aggregate(dat$abund,by=list(X=dat$x,Y=dat$y),FUN=sum)
car1<-paste("X",ppp$X,"Y",ppp$Y,sep="")


symbols(ppp$X,ppp$Y,circles=ppp$x)
points(ppp$x,ppp$y,col="red",cex=0.2)


xy<- expand.grid(Xmin:Xmax,Ymin:Ymax)
names(xy)<-c("X","Y")
xy$Z<-rep(0,length(xy$X))
car2<-paste("X",xy$X,"Y",xy$Y,sep="")

xy$Z[is.element(car2,car1)]<-ppp$x

 
symbols(xy$X,xy$Y,circles=xy$Z)
points(xy$x,xy$y,col="red",cex=0.2)



image(0:451,0:540,-matrix(log(xy$Z),452,541))

xy2<-matrix(xy$Z,452,541)
dimnames(xy2)<-list(X=0:451,Y=0:540)


########################################################################################################################################################################
#############                                     try to do some interpolation 
########################################################################################################################################################################
# 
#
# method : for each point on the grid look what is on the survey track for a scare of X nm of side  
# space between transects is about 30 nm 

#     sample(all the points , 1, replace = FALSE, prob = inverse of the distance)
# distance = ()



pnul<-0
cuales<-sample(c(1:length(xy$Z)),round(length(xy$Z)*(1-pnul),0),replace=F)

xy$Z<-0

for (i in cuales)  
{
dists<-((xy$X[i]-ppp$X)^2+(xy$Y[i]-ppp$Y)^2)^0.5
dists[dists==0]<-0.01
xy$Z[i]<-sample(ppp$x,1,replace=F, prob=1/dists) *(1-pnul)

}


xy2<-matrix(xy$Z,452,541)
dimnames(xy2)<-list(X=0:451,Y=0:540)

save(xy2,file=paste("interpolHerring",yr,".RData",sep="")
