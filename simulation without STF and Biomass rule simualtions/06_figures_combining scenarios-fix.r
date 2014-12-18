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

#library(FLSAM)
library(MASS)
wine <- F

library(FLCore)
library(FLFleet)
#library(PBSadmb)
library(lattice)
library(MASS)

ac<-function(x) {return(as.character(x))}
an<-function(x) {return(as.numeric(x))}

path          <- getwd()      #"W:/IMARES/Data/ICES-WG/WKMACLTMP/"
inPath        <- "Data/"      #"W:/IMARES/Data/ICES-WG/WKMACLTMP/Data/"
codePath      <- "R code/"    #"W:/IMARES/Data/ICES-WG/WKMACLTMP/R code/"
outPath       <- "Results/"   #"W:/IMARES/Data/ICES-WG/WKMACLTMP/Results/"
plotPath      <- "Plots"

                              # if(substr(R.Version()$os,1,3)== "lin"){
                              #   path        <- sub("W:/","/media/n/",path)
                              #   inPath      <- sub("W:/","/media/n/",inPath)
                              #   codePath    <- sub("W:/","/media/n/",codePath)
                              #   outPath     <- sub("W:/","/media/n/",outPath)
                              # }
                              
                              
                              # home<-F
                              # if(home)
                              # {
                              # path          <- "D://MSE/"
                              # inPath        <- "D://MSE/Data/"
                              # codePath      <- "D://MSE/R code/"
                              # outPath       <- "D://MSE/Results/"
                              # }



perm<-T
if (perm)cat("!!! scenario with permanent changes")


outPathp<-outPath
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
#load(file=paste(outPath,"resNFinal_1000.RData",      sep=""))
#load(file=paste(outPath,"resFFinal_1000.RData",      sep=""))
load(file=paste(outPath,"SRmod.RData",          sep=""))

for(i in 1:length(settings)) assign(x=names(settings)[i],value=settings[[i]])

source(paste(codePath,'functions_modified.r',sep=""))
source(paste(codePath,"04_forecastScenarios.r", sep=""))

RecType<-settings$RecType
projPeriod<-ac(2014:2042)

outPathbase       <- "R code/GRID-RESULTS-1000"


windows(w=12,h=9)
mylayout=layout(matrix(c(
                  0,25,26,27,28,
                  29,1,2,3,4,
                  30,5,6,7,8,
                  31,9,10,11,12,
                  32,13,14,15,16,
                  33,17,18,19,20,
                  34,21,22,23,24
                    ),
                  ncol=7),
                  width=c(1,3.2,rep(3,5)),
                  height=c(0.5,3,3,3,3.3))
layout.show(mylayout)

fb<-rep(NA,6)

opts<-c(4,6,7,11,16,21)
if (perm==F)    opts<-c(4,6,10,11,16,21)             ### selection of the management scenarios to plot



for (nb in 1:6)
{

#scen          <- c("LTMP2.2mt")              # 
  scen          <- c("LTMP3.0mt") 
opt           <- opts[nb]                      # for multiple scenario combinations a counter
TACvarlim     <- T                             # whether or not to apply the 20% limit on TAC change
Fvarlim       <- F                             # whether or not to apply the 10% limit on Fbar change
BBscen        <- "noBB"                        # banking borrowing options :
LastRecOp     <- "geom"                  # option to replace the last estimated recruitment : "SAM", "geom", "RCT3"
mpOptions<-list(opt=opt,TACvarlim=TACvarlim,Fvarlim=Fvarlim,BBscen=BBscen)
sc.name<-paste("perm",scen,opt,"_TACvarlim",TACvarlim,"_Fvarlim",Fvarlim,"_",BBscen,"_LastRec",LastRecOp,"_fix.prop",sep="")
#outPath<-paste(outPathbase,"HCR base case/",sc.name,"/",sep="")
outPath<-paste(outPathbase,"/",sc.name,sep="")
source(paste(codePath,"07_scenarioDescription_fix.r", sep=""))
mpPoints      <- get(scen)[[which(names(get(scen))==paste("opt",opt,sep=""))]]
if (perm==F)  outPath<-paste(outPathbase,"HCR sensitivity perm/",sc.name,"/",sep="")

 load(file=paste(outPath,"/",scen,opt,"_Finalf.RData",           sep=""))
 load(file=paste(outPath,"/",scen,opt,"_Finalbiol.RData",           sep=""))
 load(file=paste(outPath,"/",scen,opt,"_Finalfishery.RData",           sep=""))
 load(file=paste(outPath,"/",scen,opt,"_FinalmpPoints.RData",           sep=""))
 load(file=paste(outPath,"/",scen,opt,"_Finalsettings.RData",           sep=""))
source(paste(codePath,'functions_modified.r',sep=""))
source(paste(codePath,"04_forecastScenarios.r", sep=""))
fb[nb]<-mpPoints$alpha
RecType<-settings$RecType
nits<- dim(stocks@stock.n)[6]


projPeriod<-settings$projPeriod
projPeriod<-projPeriod[-length(projPeriod)]


ca    <- 1.1
cl    <- 1.3
fonts <- 2

  #-------------------------------------------------------------------------------
  # 0): Figure of biol future weights and recruitment
  #-------------------------------------------------------------------------------
plotPathsc<-paste( "Plots/",scen,opt,"_fix.prop_",sep="")
dir.create(plotPathsc)



############### 

rSSBp <- apply(ssbb(biol,unitSums(f),stocks)@.Data,1:5,quantile,probs=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95),na.rm=T)
rFp   <- apply(apply(unitSums(f)[ac(4:8),],2:6,mean,na.rm=T)@.Data,1:5,quantile,probs=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95),na.rm=T)
rLandf<- apply(computeLandings(fishery)@.Data,1:5,quantile,probs=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95),na.rm=T)
rB2 <-biol@n[-c(1,2),]*biol@wt[-c(1,2),]


  #- Plot settings
yrs   <- 2014:max(an(projPeriod))
cl    <- 1.2
ca    <- 1.2
fonts <- 2

 
 

#-------------------------------------------------------------------------------
# 3): Plot results of SSB stock and SSB pop
#-------------------------------------------------------------------------------
par(mar=c(1,ifelse(nb==1,2,1),0.1,0.1))
xrange  <- range(yrs)
yrange  <- c(0,8e6)
#if (perm ==F)    yrange  <- c(0,9e6)
  #---------
  #- Biology
  #---------  
plot(rSSBp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
     type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")

if(nb==1) axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1e6,cex.lab=cl,cex.axis=ca,font=fonts)
grid(lty="solid"); box()
#-Reference level Blim & Bpa
abline(h=1.84e6,col="green",lwd=2,lty=2);
#abline(h=2.36e6,col="darkgreen",lwd=2,lty=2);
#mtext(text="Blim",side=4,at=1.84e6,las=1,cex=0.8,col="blue",line=0.5,font=fonts)
#mtext(text="Bpa",side=4,at=2.36e6,las=1,cex=0.8,col="darkgreen",line=0.5,font=fonts)

polygon(c(ac(yrs),rev(ac(yrs))),c(rSSBp["5%",,ac(yrs),,,],rev(rSSBp["95%",,ac(yrs),,,])),col=rgb(0,0,1,alpha=0.08),lty=0)
polygon(c(ac(yrs),rev(ac(yrs))),c(rSSBp["10%",,ac(yrs),,,],rev(rSSBp["90%",,ac(yrs),,,])),col=rgb(0,0,1,alpha=0.15),lty=0)
polygon(c(ac(yrs),rev(ac(yrs))),c(rSSBp["25%",,ac(yrs),,,],rev(rSSBp["75%",,ac(yrs),,,])),col=rgb(0,0,1,alpha=0.22),lty=0)
ssbit<-ssbb(biol[,ac(yrs),,,,1],unitSums(f[,ac(yrs),,,,1]),stocks[,ac(yrs),,,,1])@.Data
lines(c(ssbit)~yrs,lty=1,lwd=2,col="grey55")
lines(rSSBp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=3,col="blue")
lines(rSSBp["5%",,ac(yrs),,,]~yrs,lty=1,lwd=1,col="blue")



  #---------
  #- Fishery
  #---------
par(mar=c(1,ifelse(nb==1,2,1),0.1,0.1))
iFsh<-1
  xrange  <- range(yrs)
  yrange  <- c(0,1.5e6)
#if (perm ==F)    yrange  <- c(0,2e6)
    #- Landings
  plot(rLandf["50%",,ac(yrs),iFsh,,]~yrs,xlim=xrange,ylim=yrange,xlab="Year",ylab="",
       type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
if(nb==1) axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1e6,cex.lab=cl,cex.axis=ca,font=fonts)
grid(lty="solid"); box()
polygon(c(ac(yrs),rev(ac(yrs))),c(rLandf["5%",,ac(yrs),,,],rev(rLandf["95%",,ac(yrs),,,])),col=rgb(1,0,0,alpha=0.08),lty=0)
polygon(c(ac(yrs),rev(ac(yrs))),c(rLandf["10%",,ac(yrs),,,],rev(rLandf["90%",,ac(yrs),,,])),col=rgb(1,0,0,alpha=0.15),lty=0)
polygon(c(ac(yrs),rev(ac(yrs))),c(rLandf["25%",,ac(yrs),,,],rev(rLandf["75%",,ac(yrs),,,])),col=rgb(1,0,0,alpha=0.22),lty=0)
fbarit<-        computeLandings(fishery)[,ac(yrs),,,,1]
lines(c(fbarit)~yrs,lty=1,lwd=2,col="grey55")
lines(rLandf["50%",,ac(yrs),,,]~yrs,lty=1,lwd=3,col="red")
lines(rLandf["5%",,ac(yrs),,,]~yrs,lty=1,lwd=1,col="red")


 
 
 
 

#-------------------------------------------------------------------------------
# 4): Plot results of F stock and true F
#-------------------------------------------------------------------------------
par(mar=c(1,ifelse(nb==1,2,1),0.1,0.1))
xrange  <- range(yrs)
yrange  <- c(0,0.6)

  #---------
  #- Biology
  #---------
plot(rFp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="Years",ylab="",
     type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")

if(nb==1) axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex.lab=cl,cex.axis=ca,font=fonts)
grid(lty="solid"); box()
#-Reference level Fpa
polygon(c(ac(yrs),rev(ac(yrs))),c(rFp["5%",,ac(yrs),,,],rev(rFp["95%",,ac(yrs),,,])),col=rgb(0,1,0,alpha=0.2),lty=0)
polygon(c(ac(yrs),rev(ac(yrs))),c(rFp["10%",,ac(yrs),,,],rev(rFp["90%",,ac(yrs),,,])),col=rgb(0,1,0,alpha=0.3),lty=0)
polygon(c(ac(yrs),rev(ac(yrs))),c(rFp["25%",,ac(yrs),,,],rev(rFp["75%",,ac(yrs),,,])),col=rgb(0,1,0,alpha=0.4),lty=0)
fbarit<-        quantMeans(f[ac(4:8),ac(yrs),,,,1])@.Data
lines(c(fbarit)~yrs,lty=1,lwd=2,col="grey55")
lines(rFp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=3,col="chartreuse4")
lines(rFp["5%",,ac(yrs),,,]~yrs,lty=1,lwd=1,col="chartreuse4")






#-------------------------------------------------------------------------------
# 8): Plot results of TAC
#-------------------------------------------------------------------------------
 par(mar=c(2,ifelse(nb==1,2,1),0.1,0.1))

B   <-  ssbb(biol,unitSums(f),stocks)
ys  <-  dimnames(B)$year[which((iter(B,1))>0)]
B   <-  B[,ys]
risk<-  apply(B<mpPoints$Blim,c(1:5),sum)/nits

xrange  <- range(yrs)
yrange  <- c(0,0.15)

  plot(c(risk[1,ac(yrs)]@.Data)~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
       type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxt="n")
if(nb==1)  axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex.lab=cl,cex.axis=ca,font=fonts)
#  mtext(text="Risk B<Blim",side=2,at=0.5,outer=T,cex=cl,line=1.5)
#  mtext(text="Years",side=1,at=0.5,outer=T,cex=cl,line=1.5)
grid(lty="solid"); box()
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

}




par(mar=c(0,1,1,0.5))
plot(0,0,xaxt="n",yaxt="n",cex=0,bty="n",bg="grey60")
text(-0.5,0,"SSB (mt)",srt = 90,cex=1.2)

par(mar=c(0,1,1,0.5))
plot(0,0,xaxt="n",yaxt="n",cex=0,bty="n",bg="grey60")
text(0,0,"Catches (mt)",srt = 90,cex=1.2)

par(mar=c(0,1,1,0.5))
plot(0,0,xaxt="n",yaxt="n",cex=0,bty="n",bg="grey60")
text(0,0,"Fbar4-8",srt = 90,cex=1.2)

par(mar=c(0,1,1,0.5))
plot(0,0,xaxt="n",yaxt="n",cex=0,bty="n",bg="grey60")
text(0,0,"Risk p(SSB<Blim)",srt = 90,cex=1.2)



for (nb in 1:6)
{
par(mar=c(0,0.5,1,0.1))
plot(0,0,xaxt="n",yaxt="n",cex=0,bty="n",bg="grey60")
text(0,0,paste("Alpha = ",fb[nb]),cex=1.2)
}



savePlot(file=paste(plotPathsc,"trajec.comp.Btrig",mpPoints$Btrigger,"_lowtermweights.png",sep=""),type="png")
 dev.off()

