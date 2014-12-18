#-------------------------------------------------------------------------------
# Mackerel management plan evaluation
#
# Author: Thomas Brunel    (based on Niels Hintzen code for North Sea herring management plan evaluation)
#         IMARES, The Netherlands
#
# Performs an MSE of NEA Mackerel under different TAC scenario's  :
#   
# This script simulates biological parameteres (stocks and catch weigth, maturity, F and Mprop) using a AUTOCORperm model where predictions
# are forced to conect with the last years of the observations
#
# Date: May/June-2014
#
# Build for R2.13.2, 32bits
# RUN with R2.13.2    (32-bit)
# packages :
# FLCore 2.4
# FLAssess 2.4
# FLSAM 0.99-991
#-------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   define location, load objects, define the dimentions
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#  !!!  requires R3.0.2 to run the  fArma library
rm(list=ls())
library(FLCore)


inPath        <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/Data/"
codePath      <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/R code/"
assess.name <-  "NEA-Mac-WGWIGE-2014-V2"



load("W:/IMARES/Data/ICES-WG/WKMACLTMP/Results/stocks.RData")
load("W:/IMARES/Data/ICES-WG/WKMACLTMP/Results/fishery.RData")

an  <-function(x) {as.numeric(x)} 
ac  <-function(x) {as.character(x)} 

 histMinYr<-an(dimnames(stock.wt(stocks))$year[1])
 histMaxYr<-2013
 histPeriod<-  histMinYr:histMaxYr
 futureMaxYr<-2100
 projPeriod<-(histMaxYr+1):futureMaxYr
 nits<-1000
 ages<-ac(0:12)
 



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   simulation of stock weights 
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
 sw<-read.table(file.path(inPath,"/",assess.name,"data","sw.dat"),skip=5)
 sw$Year<-1980:2014            #  in the stockassessment.org, stock weights are needed for the current year, but thats' just a copy of th elast data year
 sw<-sw[sw$Year!=2014,]
 
 
cv<-apply(sw[,-14],2, FUN = function(x) { sd(log(x))} )
cv[1]<-1
rho <- apply(sw[,-14],2, FUN = function(x) { an(unlist(acf(log(x),plot=F)[1])[1])})
rho <- mean(rho,na.rm=T)

mu<-apply(sw[,-14],2, FUN = function(x) { mean(rev(x)[1:3])} )
mu[1]<-0
 
 # prepare the output object
 wt<-array(data=NA,dim=c(length(ages),length(histMinYr:futureMaxYr),1,1,1,nits),dimnames=list(age=ages,year=histMinYr:futureMaxYr,unit="unique",season="all",area="uniqu",iters=1:nits))
 wt[,ac(histPeriod),,,,] <-  t(sw[,-14])
 
 
# muu<-wt[,ac(projPeriod),,,,]
# muu[]<-mu
 wt[,ac(projPeriod),,,,][] <-  mu   
 
# sweep( wt[,ac(projPeriod),,,,],c(1:6),c(mu),"+")
 
 x<-wt[1,ac(projPeriod),,,,]
 x[]<- rnorm(prod(dim(x)),0,1)
 for (i in 2:dim(x)[1]) x[i,] <- rho * x[i-1,] + sqrt(1-rho^2) * x[i,]
 x2<-wt[,ac(projPeriod),,,,]
 for (a in 1:13) x2[a,,]<-x 

eps<-sweep(x2,c(1,3),cv,"*") 
 
wt[,ac(projPeriod),,,,] <-  wt[,ac(projPeriod),,,,] * exp(eps)

swt<-stock.wt(stocks)
swt<-window(swt,start=1980,end=futureMaxYr)
swt<-propagate(swt,1000)
wt<-FLQuant(wt,dim=dim(swt),dimnames=dimnames(swt))

apply(wt,1:5,function(x) {sum(is.na(x))})

plot(wt)
plot(iter(wt,1))



save(wt,file=paste(inPath,"AUTOCORpermstock.wt.RData"))



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   simulation of catch weights 
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#look at the difference between catch and stock weights
 sw<-read.table(file.path(inPath,"/",assess.name,"data","sw.dat"),skip=5)
 sw$Year<-1980:2014
 sw<-sw[sw$Year!=2014,]
 
 swt<-expand.grid(Year=1980:2013, Age=0:12)
 swt$swt<-NA
 for (ag in 0:12)  swt$swt[swt$Age==ag]<-sw[,ag+1]

 cw<-read.table(file.path(inPath,"/",assess.name,"data","cw.dat"),skip=5)
 cw$Year<-1980:2013
 
 w<-expand.grid(Year=1980:2013, Age=0:12)
 w$cwt<-NA
 for (ag in 0:12)  w$cwt[w$Age==ag]<-cw[,ag+1]

 w<-merge(swt,w)

 w$dev<-w$cwt-w$swt

library(lattice)
xyplot(dev~Year|Age,w)
histogram(~dev|Age,data=w)

devmu<-aggregate(dev~Age ,w,mean)
devsd<-aggregate(dev~Age ,w,sd)

# generate random deviations based on the historic ones, and add them to the stock weights to get the catch weights corresponding to the stock 
# weights. 
dev<-wt

  for(a in 0:12)  dev[ac(a),]<- rnorm(nits*dim(dev)[2],devmu$dev[devmu$Age==a],devsd$dev[devmu$Age==a])
cw<-wt+dev 
cw[,ac(histPeriod)]   <-  stocks@catch.wt[,ac(histPeriod)]

apply(cw,1:5,function(x) {sum(is.na(x))})
plot(cw)

save(cw,file=paste(inPath,"AUTOCORpermcatch.wt.RData"))
 
 
 
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   simulation of maturity  option 1 : model as an AUTOCORperm process independent of the weights
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
 mat<-read.table(file.path(inPath,"/",assess.name,"data","mo.dat"),skip=5)
 mat$Year<-1980:2014
 mat<-mat[mat$Year!=2014,]
 
 w<-expand.grid(Year=1980:2013, Age=0:12)
 w$mat<-NA
 for (ag in 0:12)  w$mat[w$Age==ag]<-mat[,ag+1]
 
 # prepare the output object
 mato<-array(data=NA,dim=c(length(ages),length(histMinYr:futureMaxYr),1,1,1,nits),dimnames=list(age=ages,year=histMinYr:futureMaxYr,unit="unique",season="all",area="uniqu",iters=1:nits))
 mato[,ac(histPeriod),,,,] <-  t(mat[,-14])
 
 
 
cv<-apply(mat[,-14],2, FUN = function(x) { sd(log(x))} )
cv[1]<-1
rho <- apply(mat[,2:4],2, FUN = function(x) { an(unlist(acf(log(x),plot=F)[1])[1])})
rho <- mean(rho,na.rm=T)

mu<-apply(mat[,-14],2, FUN = function(x) { mean(rev(x)[1:3])} )

 
 
# muu<-wt[,ac(projPeriod),,,,]
# muu[]<-mu
mato[,ac(projPeriod),,,,][] <-  mu   
 
# sweep( wt[,ac(projPeriod),,,,],c(1:6),c(mu),"+")
 
 x<-mato[1,ac(projPeriod),,,,]
 x[]<- rnorm(prod(dim(x)),0,1)
 for (i in 2:dim(x)[1]) x[i,] <- rho * x[i-1,] + sqrt(1-rho^2) * x[i,]
 x2<-mato[,ac(projPeriod),,,,]
for (a in 1:13) x2[a,,]<-x 

eps<-sweep(x2,c(1,3),cv,"*") 
 
 
mato[,ac(projPeriod),,,,] <-  mato[,ac(projPeriod),,,,] * exp(eps)
mato[1,,,,,]<-0
mato[5:13,,,,,]<-1
 
mat<-mat(stocks)
mat<-window(mat,start=1980,end=futureMaxYr)
mat<-propagate(mat,1000)
mat<-FLQuant(mato,dim=dim(mat),dimnames=dimnames(mat))

apply(mat,1:5,function(x) {sum(is.na(x))})
plot(mat)

save(mat,file=paste(inPath,"AUTOCORpermmat.RData"))



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   simulation of m prop
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
 dt<-read.table(file.path(inPath,"/",assess.name,"data","pm.dat"),skip=5)
 dt$Year<-1980:2014
 dt<-dt[dt$Year!=2014,]
 
 w<-expand.grid(Year=1980:2013, Age=0:12)
 w$pm<-NA
 for (ag in 0:13)  w$pm[w$Age==ag]<-dt[,ag+1]
 
 # prepare the output object
 mp<-array(data=NA,dim=c(length(ages),length(histMinYr:futureMaxYr),1,1,1,nits),dimnames=list(age=ages,year=histMinYr:futureMaxYr,unit="unique",season="all",area="uniqu",iters=1:nits))
 mp[,ac(histPeriod),,,,] <-  t(dt[,-14])
 
 # fill in a 0 for age 0 
  
cv<-apply(dt[,-14],2, FUN = function(x) { sd(log(x))} )

rho <- apply(dt[,-14],2, FUN = function(x) { an(unlist(acf(log(x),plot=F)[1])[1])})
rho <- mean(rho,na.rm=T)

mu<-apply(dt[,-14],2, FUN = function(x) { mean(rev(x)[1:3])} )

 
 # prepare the output object
 wt<-array(data=NA,dim=c(length(ages),length(histMinYr:futureMaxYr),1,1,1,nits),dimnames=list(age=ages,year=histMinYr:futureMaxYr,unit="unique",season="all",area="uniqu",iters=1:nits))
 wt[,ac(histPeriod),,,,] <-  t(sw[,-14])
 
 
# muu<-mp[,ac(projPeriod),,,,]
# muu[]<-mu
 mp[,ac(projPeriod),,,,][] <-  mu   
 
# sweep( mp[,ac(projPeriod),,,,],c(1:6),c(mu),"+")
 
 x<-mp[1,ac(projPeriod),,,,]
 x[]<- rnorm(prod(dim(x)),0,1)
 for (i in 2:dim(x)[1]) x[i,] <- rho * x[i-1,] + sqrt(1-rho^2) * x[i,]
 x2<-mp[,ac(projPeriod),,,,]
 for (a in 1:13) x2[a,,]<-x 

eps<-sweep(x2,c(1,3),cv,"*") 
 
mp[,ac(projPeriod),,,,] <-  mp[,ac(projPeriod),,,,] * exp(eps)


 
mprop<-m.spwn(stocks)
mprop<-window(mprop,start=1980,end=futureMaxYr)
mprop<-propagate(mprop,1000)
mprop<-FLQuant(mp,dim=dim(mprop),dimnames=dimnames(mprop))

apply(mprop,1:5,function(x) {sum(is.na(x))})
plot(mprop)
save(mprop,file=paste(inPath,"AUTOCORpermmprop.RData"))


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   simulation of f prop using arma model  / computing it from the simulated m prop would require knowing (or simulating) how can are distributed in the year per age class.
#                                          / there is therefore no connection between changes in f and m prop here   
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


 dt<-read.table(file.path(inPath,"/",assess.name,"data","pf.dat"),skip=5)
 dt$Year<-1980:2014
 dt<-dt[dt$Year!=2014,]
 
 w<-expand.grid(Year=1980:2013, Age=0:12)
 w$pf<-NA
 for (ag in 0:12)  w$pf[w$Age==ag]<-dt[,ag+1]
 
 # prepare the output object
 fp<-array(data=NA,dim=c(length(ages),length(histMinYr:futureMaxYr),1,1,1,nits),dimnames=list(age=ages,year=histMinYr:futureMaxYr,unit="unique",season="all",area="uniqu",iters=1:nits))
 fp[,ac(histPeriod),,,,] <-  t(dt[,-14])
 

# muu<-fp[,ac(projPeriod),,,,]
# muu[]<-mu
   
cv<-apply(dt[,-14],2, FUN = function(x) { sd(log(x))} )
cv[1]<-1
rho <- apply(dt[,-14],2, FUN = function(x) { an(unlist(acf(log(x),plot=F)[1])[1])})
rho <- mean(rho[-1],na.rm=T)

mu<-apply(dt[,-14],2, FUN = function(x) { mean(rev(x)[1:3])} )


 fp[,ac(projPeriod),,,,][] <-  mu   
 
# sweep( fp[,ac(projPeriod),,,,],c(1:6),c(mu),"+")
 
 x<-fp[1,ac(projPeriod),,,,]
 x[]<- rnorm(prod(dim(x)),0,1)
 for (i in 2:dim(x)[1]) x[i,] <- rho * x[i-1,] + sqrt(1-rho^2) * x[i,]
 x2<-fp[,ac(projPeriod),,,,]
 for (a in 1:13) x2[a,,]<-x 

eps<-sweep(x2,c(1,3),cv,"*") 
 
fp[,ac(projPeriod),,,,] <-  fp[,ac(projPeriod),,,,] * exp(eps)


 
   
fprop<-harvest.spwn(stocks)
fprop<-window(fprop,start=1980,end=futureMaxYr)
fprop<-propagate(fprop,1000)
fprop<-FLQuant(fp,dim=dim(fprop),dimnames=dimnames(fprop))

apply(fprop,1:5,function(x) {sum(is.na(x))})

plot(fprop)

save(fprop,file=paste(inPath,"AUTOCORpermfprop.RData"))



