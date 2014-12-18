#-------------------------------------------------------------------------------
# WKHERMPII
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
# Performs an MSE of North Sea Herring under different TAC scenario's
#
# Date: 03-Oct-2011
#
# Build for R2.8.1
#-------------------------------------------------------------------------------

rm(list=ls())
library(FLCore)
library(FLAssess)
library(FLSAM)
library(MASS)
library(msm)

wine <- F
 path          <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/"
inPath        <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/Data/"
codePath      <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/R code/"
outPath       <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/Results/"
plotPath      <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/Results/Plots input/"
if(substr(R.Version()$os,1,3)== "lin"){
  path        <- sub("W:/","/media/n/",path)
  inPath      <- sub("W:/","/media/n/",inPath)
  codePath    <- sub("W:/","/media/n/",codePath)
  outPath     <- sub("W:/","/media/n/",outPath)
}


home<-F
if(home)
{
path          <- "D://MSE/"
inPath        <- "D://MSE/Data/"
codePath      <- "D://MSE/R code/"
outPath       <- "D://MSE/Results/"
}



perm<-T
if (perm)cat("!!! scenario with permanent changes")



if (perm)  outPathp <- paste(outPath,"perm",sep="")
#- Load objects
load(file=paste(outPath,"Mac.RData",            sep=""))
load(file=paste(outPath,"Macctrl.RData",        sep=""))
load(file=paste(outPathp,"biol.RData",          sep=""))
load(file=paste(outPathp,"fishery.RData",       sep=""))
load(file=paste(outPath,"propN.RData",          sep=""))
load(file=paste(outPath,"propWt.RData",         sep=""))
load(file=paste(outPath,"ctch.RData",          sep=""))
load(file=paste(outPath,"surveys.RData",       sep=""))
load(file=paste(outPathp,"stocks.RData",        sep=""))
load(file=paste(outPath,"settings.RData",       sep=""))
load(file=paste(outPath,"resNFinal.RData",      sep=""))
load(file=paste(outPath,"resFFinal.RData",      sep=""))
load(file=paste(outPath,"SRmod.RData",          sep=""))
load(file=paste(outPath,"resFFinal.RData",      sep=""))
for(i in 1:length(settings)) assign(x=names(settings)[i],value=settings[[i]])


source(paste(codePath,"functions.r",            sep=""))
source(paste(codePath,"04_forecastScenarios.r", sep=""))

RecType<-settings$RecType
projPeriod<-ac(2014:2042)

scen          <- c("LTMP2.2mt")              # 
opt           <- 8                      # for multiple scenario combinations a counter
TACvarlim     <- T                      # whether or not to apply the 20% limit on TAC change
Fvarlim       <- T                      # whether or not to apply the 10% limit on Fbar change
BBscen        <- "noBB"        # banking borrowing options :
                                                     # "Banking"          : always bank 
                                                     # "Borrowing"        : always borrow
                                                     # "AlternateBank"    : bank first and then alternate
                                                     # "AlternateBorrow"  : borrow first and then alternate
                                                     # "MinVar"           : use BB to minise TAC variability
                                                     
LastRecOp     <- "geom"                  # option to replace the last estimated recruitment : "SAM", "geom", "RCT3"
                                                     # "SAM"  = don't overwrite the SAM estimmate
                                                     # "SAM"  = don't overwrite the SAM estimmate
                                                     # "geom" = replace by geomean 1990/(TaY-1)
                                                     # "RCT3" = replace by RCT3 output
                                                     
                                                     
                                                     
mpOptions<-list(opt=opt,TACvarlim=TACvarlim,Fvarlim=Fvarlim,BBscen=BBscen)

sc.name<-paste("HCR sensitivity recovery from Blim/",scen,opt,"_TACvarlim",TACvarlim,"_Fvarlim",Fvarlim,"_",BBscen,"_LastRec",LastRecOp,sep="")

outPath<-paste(outPath,sc.name,"/",sep="")

test<-F
if(test)
{
outPath<-paste(path,"Results/test/",sep="")
scen<-"LTMP2.2mt15"
}


load(file=paste(outPath,"/",scen,opt,"_Finalf.RData",           sep=""))
load(file=paste(outPath,"/",scen,opt,"_Finalbiol.RData",           sep=""))
load(file=paste(outPath,"/",scen,opt,"_Finalstocks.RData",           sep=""))
load(file=paste(outPath,"/",scen,opt,"_Finalpercievedstocks.RData",           sep=""))
load(file=paste(outPath,"/",scen,opt,"_FinalTAC.RData",           sep=""))
load(file=paste(outPath,"/",scen,opt,"_FinalfSTF.RData",           sep=""))
load(file=paste(outPath,"/",scen,opt,"_Finalfishery.RData",           sep=""))
load(file=paste(outPath,"/",scen,opt,"_FinalmpPoints.RData",           sep=""))
load(file=paste(outPath,"/",scen,opt,"_Finaldepleted.RData",           sep=""))
load(file=paste(outPath,"/",scen,opt,"_Finalsettings.RData",           sep=""))
load(file=paste(outPath,"/",scen,opt,"_FinalSSB.RData",           sep=""))
source(paste(codePath,"functions.r",            sep=""))
source(paste(codePath,"04_forecastScenarios.r", sep=""))

RecType<-settings$RecType
nits<- dim(stocks@stock.n)[6]


projPeriod<-settings$projPeriod
projPeriod<-projPeriod[-length(projPeriod)]


##-------------------------------------------------------------------------------
## Figures on stochastic variables
##-------------------------------------------------------------------------------

ca    <- 1.1
cl    <- 1.3
fonts <- 2

  #-------------------------------------------------------------------------------
  # 0): Figure of biol future weights and recruitment
  #-------------------------------------------------------------------------------
plotPath<-paste(outPath,"Plots input/",sep="")
dir.create(plotPath)




#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Figures on results
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
plotPath2       <- paste(outPath,"Plot results/",sep="")
dir.create(plotPath2)





#-------------------------------------------------------------------------------
# stock trajectories figures
#-------------------------------------------------------------------------------


############### 

rSSBp <- apply(ssbb(biol,unitSums(f),stocks)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rSSBs <- apply(ssb(stockstore)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rFp   <- apply(apply(unitSums(f)[ac(4:8),],2:6,mean,na.rm=T)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rFs   <- apply(fbar(stockstore)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rFstf <- apply(apply(unitSums(fSTF)[ac(4:8),],2:6,mean,na.rm=T)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rLandf<- apply(computeLandings(fishery)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rLands<- apply(computeLandings(stockstore)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rTAC  <- apply(TAC@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)

  #- Plot settings
yrs   <- 2001:max(an(projPeriod))
cl    <- 1.2
ca    <- 1.1
fonts <- 2

 
 
#-------------------------------------------------------------------------------
# 3): Plot results of SSB stock and SSB pop
#-------------------------------------------------------------------------------

# define the target for recovery
Recov<-mpPoints$Bpa     

Ssb<-ssbb(biol,f,stocks)[,projPeriod]

par(mfrow=c(3,3),mar=c(3,3,3,3),oma=c(3.1,3.1,1,1))
its<-sample(1:nits,9,replace=F)
for (i in its)
{
xrange  <- range(an(projPeriod))
yrange  <- c(0,range(pretty(c(rSSBp[,,ac(yrs),,,],rSSBs[,,ac(yrs),,,])),na.rm=T)[2])

plot(c(iter(Ssb,i))~an(projPeriod),type="l",xlab="Years",ylab="",
     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1,col=i)
mtext(side=2,at=yrange[2]/2,text="SSB",outer=F,line=4,font=fonts,cex=ca)
abline(h=mpPoints$Blim,lwd=2,col="red",lty="dotted")
abline(h=Recov,lwd=2,col="green",lty="dotted")
legend(col=c("green","red"),lty=c("dotted","dotted"),x="topright",c("Recovery","Blim"),bg="white")
grid(); box()
}
savePlot(file=paste(plotPath2,"1.SSBtrajectories_iters.png",sep=""),type="png"); dev.off()



# recovery time
dp<-Ssb<Recov


RT<-rep(NA,nits)
for (its in 1:nits) 
{
idx <-  which(iter(depleted,its)==1)
yrsdpl<-  an( dimnames(depleted)$year[idx])[1]
idx <-  which(iter(dp,its)==0)
yrsndpl<-  an( dimnames(dp)$year[idx])
yrsndpl<-yrsndpl[yrsndpl>yrsdpl][1]
RT[its]<-yrsndpl-yrsdpl
}
plot(table(RT),xlab="recovery time",ylab="number of iterations",main="recovery from Blim")








 
# recovery time
RT<-rep(NA,nits)
for (its in 1:nits) 
{
idx <-  which(iter(depleted,its)==1)
yrsdpl<-  an( dimnames(depleted)$year[idx])[1]
idx <-  which(iter(depleted,its)==0)
yrsndpl<-  an( dimnames(depleted)$year[idx])
yrsndpl<-yrsndpl[yrsndpl>yrsdpl][1]
RT[its]<-yrsndpl-yrsdpl
}
plot(table(RT),xlab="recovery time",ylab="number of iterations",main="recovery from Blim")
savePlot(file=paste(plotPath2,"2.recovery time.png",sep=""),type="png"); dev.off()


par(mar=c(5.1,5.1,4.1,2.1))
xrange  <- range(yrs)
yrange  <- c(0,range(pretty(c(rSSBp[,,ac(yrs),,,],rSSBs[,,ac(yrs),,,])),na.rm=T)[2])
  #---------
  #- Biology
  #---------  
plot(rSSBp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="Years",ylab="",
     type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
mtext(text="SSB (thousand tonnes)",side=2,at=0.5,outer=T,cex=cl,line=-1)
grid(); box()
#-Reference level Blim & Bpa
abline(h=1.84e6,col="blue",lwd=2,lty=2);
abline(h=2.36e6,col="darkgreen",lwd=2,lty=2);
mtext(text="Blim",side=4,at=1.84e6,las=1,cex=0.8,col="blue",line=0.5,font=fonts)
mtext(text="Bpa",side=4,at=2.36e6,las=1,cex=0.8,col="darkgreen",line=0.5,font=fonts)
lines(rSSBp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rSSBp["5%",,ac(yrs),,,]~yrs,lty=3)
lines(rSSBp["95%",,ac(yrs),,,]~yrs,lty=3)

  #---------
  #- Stock
  #---------
lines(rSSBs["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
lines(rSSBs["5%",,ac(yrs),,,]~yrs,lty=3,col="red")
lines(rSSBs["95%",,ac(yrs),,,]~yrs,lty=3,col="red")

  #---------
  #- Legend
  #---------

legend("bottomright",legend=c("SSB assessed stock","SSB true population"),
       col=c("red","black"),lwd=3,lty=1,box.lty=0)
savePlot(file=paste(plotPath2,"3.SSB.png",sep=""),type="png");dev.off()








#-------------------------------------------------------------------------------
# 4): Plot results of F stock and true F
#-------------------------------------------------------------------------------
par(mar=c(5.1,5.1,4.1,2.1))
xrange  <- range(yrs)
yrange  <- c(0,range(pretty(c(rFp[,,ac(yrs),,,],rFs[,,ac(yrs),,,],rFstf[,,ac(yrs),,,])),na.rm=T)[2])

  #---------
  #- Biology
  #---------
plot(rFp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="Years",ylab="",
     type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
mtext(text="Fishing mortality (ages 4-8)",side=2,at=0.5,outer=T,cex=cl,line=-1)
grid(); box()
#-Reference level Fpa
abline(h=0.26,col="darkgreen",lwd=2,lty=2);
mtext(text="Fpa",side=4,at=0.26,las=1,cex=0.8,col="darkgreen",line=0.5,font=fonts)
lines(rFp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rFp["5%",,ac(yrs),,,]~yrs,lty=3)
lines(rFp["95%",,ac(yrs),,,]~yrs,lty=3)

  #---------
  #- Stock
  #---------
lines(rFs["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
lines(rFs["5%",,ac(yrs),,,]~yrs,lty=3,col="red")
lines(rFs["95%",,ac(yrs),,,]~yrs,lty=3,col="red")

  #---------
  #- STF
  #---------
lines(rFstf["50%",,ac(2014:rev(projPeriod)[2]),,,]~ac(2014:rev(projPeriod)[2]),lty=1,lwd=2,col="blue")
lines(rFstf["5%",,ac(2014:rev(projPeriod)[2]),,,]~ac(2014:rev(projPeriod)[2]),lty=3,col="blue")
lines(rFstf["95%",,ac(2014:rev(projPeriod)[2]),,,]~ac(2014:rev(projPeriod)[2]),lty=3,col="blue")

  #---------
  #- Legend
  #---------

legend("bottomright",legend=c("F assessed stock","F true population","F short term forecast"),
       col=c("red","black","blue"),lwd=3,lty=1,box.lty=0)
savePlot(file=paste(plotPath2,"4.F.png",sep=""),type="png");dev.off()










#-------------------------------------------------------------------------------
# 5): Plot results of landings by fleet & TAC on top
#-------------------------------------------------------------------------------


  #---------
  #- Fishery
  #---------
  par(mar=c(5.1,5.1,4.1,2.1))
iFsh<-1
  xrange  <- range(yrs)
  yrange  <- c(0,range(pretty(c(rLandf[,,ac(yrs),iFsh,,],rTAC[,,ac(yrs),iFsh,,])),na.rm=T)[2])

    #- Landings
  plot(rLandf["50%",,ac(yrs),iFsh,,]~yrs,xlim=xrange,ylim=yrange,xlab="Year",ylab="Landings (thousand tonnes)",
       type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
  axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
  grid(); box()
  lines(rLandf["50%",,ac(yrs),iFsh,,]~yrs,lty=1,lwd=2)
  lines(rLandf["5%",,ac(yrs),iFsh,,]~yrs,lty=3)
  lines(rLandf["95%",,ac(yrs),iFsh,,]~yrs,lty=3)
    #- TAC
  lines(rTAC["50%",,ac(yrs),iFsh,,]~yrs,lty=1,lwd=2,col="red")
  lines(rTAC["5%",,ac(yrs),iFsh,,]~yrs,lty=3,col="red")
  lines(rTAC["95%",,ac(yrs),iFsh,,]~yrs,lty=3,col="red")



legend("bottomright",legend=c("Landings","Advice TAC"),
       col=c("black","red"),lwd=3,lty=1,box.lty=0)
savePlot(file=paste(plotPath2,"5.CatchTAC.png",sep=""),type="png");dev.off()





#-------------------------------------------------------------------------------
# 6): Plot results total landings fishery and stock
#-------------------------------------------------------------------------------
par(mar=c(5.1,5.1,4.1,2.1))
xrange  <- range(yrs)
yrange  <- c(0,range(pretty(c(rLandf[,,ac(yrs),,,],rLands[,,ac(yrs),,,])),na.rm=T)[2])

  #---------
  #- Biology
  #---------
plot((rLandf["50%",,ac(yrs),,,])~yrs,xlim=xrange,ylim=yrange,xlab="Years",ylab="",
     type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
mtext(text="Landings (tonnes)",side=2,at=0.5,outer=T,cex=cl,line=-1)
grid(); box()
lines((rLandf["50%",,ac(yrs),,,])~yrs,lty=1,lwd=2)
lines((rLandf["5%",,ac(yrs),,,])~yrs,lty=3)
lines((rLandf["95%",,ac(yrs),,,])~yrs,lty=3)

  #---------
  #- Stock
  #---------
lines(rLands["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
lines(rLands["5%",,ac(yrs),,,]~yrs,lty=3,col="red")
lines(rLands["95%",,ac(yrs),,,]~yrs,lty=3,col="red")

  #---------
  #- Legend
  #---------

legend("bottomright",legend=c("Catch assessed stock","Catch true population"),
       col=c("red","black"),lwd=3,lty=1,box.lty=0)
savePlot(file=paste(plotPath2,"6.Catch.png",sep=""),type="png");dev.off()






#-------------------------------------------------------------------------------
# 7): Plot results of TAC
#-------------------------------------------------------------------------------

par(mar=c(2.1,2.1,2.1,2.1),oma=c(3,3,2,0))

  #---------
  #- Fishery
  #---------
iFsh<-1
  xrange  <- range(yrs)
  yrange  <- c(0,range(pretty(rTAC[,,ac(yrs),iFsh,,]),na.rm=T)[2])

  plot(rTAC["50%",,ac(yrs),iFsh,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
       type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
  axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
  mtext(text="TAC (thousand tonnes)",side=2,at=0.5,outer=T,cex=cl,line=1.5)
  mtext(text="Years",side=1,at=0.5,outer=T,cex=cl,line=1.5)
  grid(); box()
  lines(rTAC["50%",,ac(yrs),iFsh,,]~yrs,lty=1,lwd=2)
  lines(rTAC["5%",,ac(yrs),iFsh,,]~yrs,lty=3)
  lines(rTAC["95%",,ac(yrs),iFsh,,]~yrs,lty=3)


savePlot(file=paste(plotPath2,"7.TAC.png",sep=""),type="png");dev.off()



#-------------------------------------------------------------------------------
# 8): Plot results of TAC
#-------------------------------------------------------------------------------
par(mar=c(2.1,2.1,2.1,2.1),oma=c(3,3,2,0))

B   <-  ssbb(biol,unitSums(f),stocks)
ys  <-  dimnames(B)$year[which((iter(B,1))>0)]
B   <-  B[,ys]
risk<-  apply(B<mpPoints$Blim,c(1:5),sum)/nits

xrange  <- range(yrs)
yrange  <- c(0,1)

  plot(c(risk[1,ac(yrs)]@.Data)~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
       type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
  axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
  mtext(text="Risk B<Blim",side=2,at=0.5,outer=T,cex=cl,line=1.5)
  mtext(text="Years",side=1,at=0.5,outer=T,cex=cl,line=1.5)
  grid(); box()
  abline(h=0.05,lty="dotted",col="red")

ShortT    <-  ac(2014:2018)
MidT      <-  ac(2019:2028)
LongT     <-  ac(2029:2052)

  lines(an(ShortT),rep(max(c(risk[1,ShortT]@.Data)),length(ShortT)),col="blue")
  points(range(an(ShortT)),rep(max(c(risk[1,ShortT]@.Data)),2),col="blue",cex=1.2,pch="S")
  lines(an(MidT),rep(max(c(risk[1,MidT]@.Data)),length(MidT)),col="blue")
  points(range(an(MidT)),rep(max(c(risk[1,MidT]@.Data)),2),col="blue",cex=1.2,pch="M")
  lines(an(LongT),rep(max(c(risk[1,LongT]@.Data)),length(LongT)),col="blue")
  points(range(an(LongT)),rep(max(c(risk[1,LongT]@.Data)),2),col="blue",cex=1.2,pch="L")

savePlot(file=paste(plotPath2,"8.Risk.png",sep=""),type="png");dev.off()


#-------------------------------------------------------------------------------
# 8): Plot age and weights indicators 
#-------------------------------------------------------------------------------

# mean weight in the spawning stock and in the catches
smw<-quantSums((biol@n*biol@wt*biol@fec)) / quantSums((biol@n*biol@fec))
ymw<-quantSums((fishery@landings.n*fishery@landings.wt))/quantSums((fishery@landings.n))
rsmw  <- apply(smw@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rymw  <- apply(ymw@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)


par(mar=c(2.1,2.1,2.1,2.1),oma=c(3,3,2,0))
iFsh<-1
  xrange  <- range(yrs)
  yrange  <- c(range(pretty(rsmw[,,ac(yrs),iFsh,,]),na.rm=T)[1],range(pretty(rymw[,,ac(yrs),iFsh,,]),na.rm=T)[2])

  plot(rsmw["50%",,ac(yrs),iFsh,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
       type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
  axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
  mtext(text="mean weight (in kg)",side=2,at=0.5,outer=T,cex=cl,line=1.5)
  mtext(text="Years",side=1,at=0.5,outer=T,cex=cl,line=1.5)
  grid(); box()
  lines(rsmw["50%",,ac(yrs),iFsh,,]~yrs,lty=1,lwd=2)
  lines(rsmw["5%",,ac(yrs),iFsh,,]~yrs,lty=3)
  lines(rsmw["95%",,ac(yrs),iFsh,,]~yrs,lty=3)
  lines(rymw["50%",,ac(yrs),iFsh,,]~yrs,lty=1,lwd=2,col="red")
  lines(rymw["5%",,ac(yrs),iFsh,,]~yrs,lty=3,col="red")
  lines(rymw["95%",,ac(yrs),iFsh,,]~yrs,lty=3,col="red")

  legend(x="topright", col=c("black","red"),lty= c(1,1), c("in the spawning stock","in the catches"),bg="white")

savePlot(file=paste(plotPath2,"9.meanW.png",sep=""),type="png");dev.off()

# mean age in the spawning stock and in the catches
sma<-quantSums(sweep((biol@n*biol@fec),c(1:6) ,c(0:12),"*"))/quantSums((biol@n*biol@fec)) 
yma<-quantSums(sweep((fishery@landings.n),c(1:6) ,c(0:12),"*"))/quantSums((fishery@landings.n))
rsma  <- apply(sma@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
ryma  <- apply(yma@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)


par(mar=c(2.1,2.1,2.1,2.1),oma=c(3,3,2,0))
iFsh<-1
  xrange  <- range(yrs)
  yrange  <- c(range(pretty(rsma[,,ac(yrs),iFsh,,]),na.rm=T)[1],range(pretty(ryma[,,ac(yrs),iFsh,,]),na.rm=T)[2])

  plot(rsma["50%",,ac(yrs),iFsh,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
       type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
  axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
  mtext(text="mean Age",side=2,at=0.5,outer=T,cex=cl,line=1.5)
  mtext(text="Years",side=1,at=0.5,outer=T,cex=cl,line=1.5)
  grid(); box()
  lines(rsma["50%",,ac(yrs),iFsh,,]~yrs,lty=1,lwd=2)
  lines(rsma["5%",,ac(yrs),iFsh,,]~yrs,lty=3)
  lines(rsma["95%",,ac(yrs),iFsh,,]~yrs,lty=3)
  lines(ryma["50%",,ac(yrs),iFsh,,]~yrs,lty=1,lwd=2,col="red")
  lines(ryma["5%",,ac(yrs),iFsh,,]~yrs,lty=3,col="red")
  lines(ryma["95%",,ac(yrs),iFsh,,]~yrs,lty=3,col="red")

  legend(x="topright", col=c("black","red"),lty= c(1,1), c("in the spawning stock","in the catches"),bg="white")

savePlot(file=paste(plotPath2,"10.meanAge.png",sep=""),type="png");dev.off()




#-------------------------------------------------------------------------------
# 9): Plot trajectories of ssb(biol), rec(biol), f(biol), TAC(A)
#-------------------------------------------------------------------------------

par(mfrow=c(2,2),mar=c(3,3,3,3),oma=c(3.1,3.1,1,1))
 its<-sample(1:nits,10,replace=F)
#- Plot ssb based on biol
xrange  <- range(yrs)
yrange  <- range(0,pretty(range(c(iter(ssbb(biol[,ac(yrs)],unitSums(f[,ac(yrs)]),stocks[,ac(yrs)]),1:10))/1000),1))
plot(c(iter(ssbb(biol[,ac(yrs)],unitSums(f[,ac(yrs)]),stocks[,ac(yrs)]),1))/1000~yrs,type="l",xlab="Years",ylab="",
     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1)
mtext(side=2,at=yrange[2]/2,text="Spawning Stock Biomass (kt)",outer=F,line=4,font=fonts,cex=ca)
grid(); box()
for(i in its) lines(c(iter(ssbb(biol[,ac(yrs)],unitSums(f[,ac(yrs)]),stocks[,ac(yrs)]),i))/1000~yrs,col=i,lwd=2)

#- Plot rec based on biol
xrange  <- range(yrs)
yrange  <- c(0,rev(pretty(range(c(iter(biol@n[1,ac(yrs)],1:10))/1e6),1))[1])
plot(c(iter(biol@n[1,ac(yrs)],1))/1e6~yrs,type="l",xlab="Years",ylab="",
     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1)
mtext(side=2,at=yrange[2]/2,text="Recruitment (millions)",outer=F,line=4,font=fonts,cex=ca)
grid();box()
for(i in its) lines(c(iter(biol@n[1,ac(yrs)],i))/1e6~yrs,col=i,lwd=2)

#- f-ages 2-6
xrange  <- range(yrs)
yrange  <- c(0,pretty(range(c(iter(quantMeans(unitSums(f[ac(4:8),ac(yrs)])),1:10))),1)[3])
plot(c(iter(quantMeans(unitSums(f[ac(4:8),ac(yrs)])),1))~yrs,type="l",xlab="Years",ylab="",
     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1)
mtext(side=2,at=yrange[2]/2,text="Fishing mortality",outer=F,line=4,font=fonts,cex=ca)
grid();box()
for(i in its) lines(c(iter(quantMeans(unitSums(f[ac(4:8),ac(yrs)])),i))~yrs,col=i,lwd=2)

#- TAC of fleet A
xrange  <- range(yrs)
yrange  <- c(0,pretty(range(c(iter(TAC[,ac(yrs),"A"]/1000,1:10))),1)[2])
plot(c(iter(TAC[,ac(yrs),"A"]/1000,1))~yrs,type="l",xlab="Years",ylab="",
     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1)
mtext(side=2,at=yrange[2]/2,text="TAC (kt)",outer=F,line=4,font=fonts,cex=ca)
grid();box()
for(i in its) lines(c(iter(TAC[,ac(yrs),"A"]/1000,i))~yrs,col=i,lwd=2)

savePlot(file=paste(plotPath2,"11.iterations.png",sep=""),type="png");dev.off()




#-------------------------------------------------------------------------------
# 9): Report figures
#-------------------------------------------------------------------------------

yrs   <- 2005:an(rev(projPeriod)[2])
cl    <- 1.1
ca    <- 1
fonts <- 1
yrangeSSB <- c(0,6e6)
yrangeLan <- c(0,2e6)
#yrangeLan2<- c(0,2.5e4)
yrangeF   <- c(0,0.5)


par(mfrow=c(3,1),oma=c(6,6,2,3),mar=c(1,0,0,0))

  #---------
  #- Landings
  #---------

xrange  <- range(yrs)
yrange  <- yrangeLan

plot(rLandf["50%",,ac(yrs),1,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
mtext(text=expression(paste("Landings (",10^3," tonnes)",sep="")),side=2,at=(yrange[2]-yrange[1])/2+yrange[1],outer=F,cex=cl,line=4,font=fonts)
grid(); box()
lines(rLandf["50%",,ac(yrs),1,,]~yrs,lty=1,lwd=2)
lines(rLandf["5%",,ac(yrs),1,,]~yrs,lty=3,lwd=1)
lines(rLandf["95%",,ac(yrs),1,,]~yrs,lty=3,lwd=1)

  #------------------
  # True F
  #------------------

xrange  <- range(yrs)
yrange  <- yrangeF
plot(rFp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
mtext(text=expression(paste(F[2-6]," (",year^-1,")",sep="")),side=2,at=(yrange[2]-yrange[1])/2+yrange[1],outer=F,cex=cl,line=4,font=fonts)
grid(); box()
lines(rFp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rFp["5%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
lines(rFp["95%",,ac(yrs),,,]~yrs,lty=3,lwd=1)

  #-Reference level Fpa
abline(h=0.26,col="darkgreen",lwd=1,lty=2);
mtext(text="Fpa",side=4,at=0.26,las=1,cex=0.65,col="darkgreen",line=0.5,font=fonts)
text(x=xrange[1],y=yrange[2],pos=1,labels="(B)",font=fonts,cex=cl)

  #------------------
  # True SSB
  #------------------

xrange  <- range(yrs)
yrange  <- yrangeSSB
plot(rSSBp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
axis(1,las=1,at=pretty(xrange),labels=pretty(xrange),cex=ca,font=fonts)
mtext(text=expression(paste("SSB (",10^3," tonnes)",sep="")),side=2,at=(yrange[2]-yrange[1])/2+yrange[1],outer=F,cex=cl,line=4,font=fonts)
grid(); box()
lines(rSSBp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rSSBp["5%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
lines(rSSBp["95%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
  #-Reference level Blim & Bpa
abline(h=1.84e6,col="blue",lwd=1,lty=2);
abline(h=2.36e6,col="darkgreen",lwd=1,lty=2);
mtext(text="Blim",side=4,at=1.84e6,las=1,cex=0.65,col="blue",line=0.5,font=fonts)
mtext(text="Bpa",side=4,at=2.36e6,las=1,cex=0.65,col="darkgreen",line=0.5,font=fonts)
text(x=xrange[1],y=yrange[2],pos=1,labels="(C)",font=fonts,cex=cl)



  #- Labels x-axis
mtext(text=expression(Years),side=1,at=(xrange[2]-xrange[1])/2+xrange[1],outer=F,cex=cl,line=4,font=fonts)
savePlot(paste(plotPath2,"12.Truth.png",sep=""),type="png") ; dev.off()
