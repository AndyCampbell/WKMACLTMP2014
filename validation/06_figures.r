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
#load(file=paste(outPath,"resNFinal_1000.RData",      sep=""))
#load(file=paste(outPath,"resFFinal_1000.RData",      sep=""))
load(file=paste(outPath,"SRmod.RData",          sep=""))
nits<-1000
for(i in 1:length(settings)) assign(x=names(settings)[i],value=settings[[i]])


source(paste(codePath,"functions.r",            sep=""))
source(paste(codePath,"04_forecastScenarios.r", sep=""))

RecType<-settings$RecType
projPeriod<-ac(2014:2042)

scen          <- c("surrogateFbar")              # 
opt           <- ""                      # for multiple scenario combinations a counter
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
mpPoints<-list()                                                     
mpPoints$Fequ<-"FS"                                                     
prlength<-34                                                     
mpOptions<-list(opt=opt,TACvarlim=TACvarlim,Fvarlim=Fvarlim,BBscen=BBscen)

sc.name<-paste(scen,opt,"_TACvarlim",TACvarlim,"_Fvarlim",Fvarlim,"_",BBscen,"_LastRec",LastRecOp,sep="")

outPath<-paste(outPath,"validation/",sep="")





load(file=paste(outPath,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalf.RData",        sep=""))
load(file=paste(outPath,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalbiol.RData",        sep=""))
load(file=paste(outPath,"surrogateFbarFS1000its34yrs_Finalpercievedstocks.RData",        sep=""))
load(file=paste(outPath,"surrogateFbarFS1000its34yrs_Finalfishery.RData",        sep=""))


source(paste(codePath,"functions.r",            sep=""))
source(paste(codePath,"04_forecastScenarios.r", sep=""))

RecType<-settings$RecType
nits<- dim(stocks@stock.n)[6]


projPeriod<-ac(2014:2046)


##-------------------------------------------------------------------------------
## Figures on stochastic variables
##-------------------------------------------------------------------------------

ca    <- 1.1
cl    <- 1.3
fonts <- 2

  #-------------------------------------------------------------------------------
  # 0): Figure of biol future weights and recruitment
  #-------------------------------------------------------------------------------
plotPathsc<-paste( "W:/IMARES/Data/ICES-WG/WKMACLTMP/Results/validation/",sep="")
dir.create(plotPathsc)
plotPath<-paste(plotPathsc,"/Plots input/",sep="")
dir.create(plotPath)



    #- Recruitment lognormal default fit
hist(biol@n[1,projPeriod],breaks=100,xlab="Recruitment",ylab="Simulated frequency",main="",cex.lab=cl,cex.axis=ca,font=fonts)
abline(v=exp(mean(log(biol@n[1,recrPeriod]))),col="red",lty=3,lwd=3)
legend("topright",legend=c("Simulated recruitment","Geometric mean recruitment"),lty=c(0,2),pch=c(22,0),
       pt.bg="white",pt.cex=c(1,0),box.lty=0,col=c("black","red"),lwd=c(0,3))
legend("bottomright",lty=c(0,0),c(paste("Btrigger=",mpPoints$Btrigger),paste("Ftarget=",mpPoints$Ftarget)))
savePlot(file=paste(plotPath,"input_recruitment.png",sep=""),type="png") ; dev.off()

    #- Recruitment cummulative distribution
windows(w=12,h=6)
xl<-c(0,2e7)
par(mfrow=c(1,2))
plot(x=sort(c(rec(Mac[,ac(recrPeriod)]))),y=seq(0,1,length.out=length(c(rec(Mac[,ac(recrPeriod)])))+2)[-c(1,length(c(rec(Mac[,ac(recrPeriod)])))+2)],
     ylim=c(0,1),col=2,pch=15,xlab="Recruitment",ylab="Cumulative Probability",
     main="Observed/ simulated values",las=1,xlim=xl,yaxs="i")
points(x=sort(c(biol@n[1,ac(projPeriod)])),y=seq(0,1,length.out=length(sort(c(biol@n[1,ac(projPeriod)])))+2)[c(-1,-length(sort(c(biol@n[1,ac(projPeriod)]))))],pch=19,cex=0.3)
legend("bottomright",legend=c("Simulated recruitment","Observed recruitment"),lty=c(0,0),pch=c(16,15),
       pt.cex=c(0.5,1),box.lty=0,col=c("black","red"),lwd=c(0,0))

ys<-sort(c(rec(Mac[,ac(recrPeriod)])))
rs<-c(biol@n[1,ac(projPeriod)])
qs<-c(1:c(length(recrPeriod)))/c(length(recrPeriod)+1)
xs<-quantile(rs,probs=qs)

plot(xs,ys,cex.lab=cl,col="red",pch=15,xlim=xl,ylim=xl,main="QQ plot", xlab="Simulated", ylab="Observed")
abline(0,1)
savePlot(file=paste(plotPath,"input_recruitmentDistri.png",sep=""),type="png") ; dev.off()

    #- SR plot with simulated data
Recs<-c(biol@n[1,ac(projPeriod)]@.Data)
biol@n[1,"2047"]<-1
SSBs<-c(ssbb(biol[,ac(projPeriod)],f[,ac(projPeriod)],stockstore[,ac(projPeriod)])@.Data)/1e6
RecO<-c(rec(Mac[,ac(recrPeriod)]))
SSBO<-c(ssb(Mac[,ac(recrPeriod)]))/1e6


png(file=paste(plotPath,"input_simulated recruitments.png",sep=""),width=8,height=6,units="in",res=300)
plot(SSBs,Recs,xlim=c(0,8),xlab="SSB (10^6t)",ylab="Rec (thousands)",cex=0.4,pch=19,main=paste("Simulated values"))
seqMin <- 0.5
seqStep <- 0.5
Xmax <-8
up=seq(seqMin,Xmax,seqStep)
lw=seq(seqMin,Xmax,seqStep)
md=seq(seqMin,Xmax,seqStep)
ssb=seq(seqMin,Xmax,seqStep)
loopNum <- length(ssb)    
for (j in 1:loopNum){     
  up[j]=quantile(Recs[((SSBs>up[j])*(SSBs<up[j+1]))>0],probs=.95,na.rm=TRUE)
  lw[j]=quantile(Recs[((SSBs>lw[j])&(SSBs<lw[j+1]))>0],probs=.05,na.rm=TRUE)
  md[j]=quantile(Recs[((SSBs>md[j])&(SSBs<md[j+1]))>0],probs=.5,na.rm=TRUE)
  }
lines(ssb,md[1:loopNum],col="green",lwd=3)
lines(ssb,up[1:loopNum],col=4,lwd=3)
lines(ssb,lw[1:loopNum],col=4,lwd=3)
points(SSBO,RecO,col="red",cex=1,pch=19)

dev.off()





Recs<-biol@n[1,]
SSBs<-ssbb(biol,f,stockstore)/1e6


png(file=paste(plotPath,"input_simulated recruitments_iterations.png",sep=""),width=12,height=8,units="in",res=300)
par(mfrow=c(3,4))
for (its in sample(1:1000,12,replace=F))
{
plot(SSBs[,ac(projPeriod),,,,its],Recs[,ac(projPeriod),,,,its],xlim=c(0,6)
                ,ylim=c(0,2e7) ,pch=19,
                xlab="SSB (in mt)",ylab="Rct (thousands)",
                main=paste("iteration",its))
points(SSBs[,ac(recrPeriod),,,,its],Recs[,ac(recrPeriod),,,,its],col="red",pch=19 )
bs<-seq(0,6e6,length.out=1000)
mod<-B2Rsimple(bs,SRmod[its,])
lines(bs/1e6,mod)

}
dev.off()
#
#
#RecO<-c(rec(Mac[,ac(recrPeriod)]))
#SSBO<-c(ssb(Mac[,ac(recrPeriod)]))/1e6
#Rec<-c(rec(Mac))
#SSB<-c(ssb(Mac))/1e6
#
#plot(SSB,Rec,xlim=c(0,6),ylim=c(0,1.2e7),xlab="SSB (10^6t)",ylab="Rec (thousands)",cex=0.4,pch=19,main="")
#points(SSBO,RecO,col="red",cex=1,pch=19)
#





    #- Numbers-at-age
   
par(mfrow=c(4,6),mar=c(0.3,0.1,0.3,0.1),oma=c(6,6,2,2),las=1)
  for(iYr in 1990:2013){
    if(iYr %in% seq(1990,2007,6)) boxplot(data ~ as.factor(age),data=as.data.frame(biol@n[,ac(iYr)]),xaxt="n",xlab="",ylab="",ylim=c(0,rev(pretty(max(as.data.frame(biol@n[,ac(1990:2013)])$data,na.rm=T)))[1]))
    if(iYr %in% 2009:2013)        boxplot(data ~ as.factor(age),data=as.data.frame(biol@n[,ac(iYr)]),yaxt="n",xlab="",ylab="",ylim=c(0,rev(pretty(max(as.data.frame(biol@n[,ac(1990:2013)])$data,na.rm=T)))[1]))
    if(iYr == 2008)               boxplot(data ~ as.factor(age),data=as.data.frame(biol@n[,ac(iYr)]),xlab="",ylab="",ylim=c(0,rev(pretty(max(as.data.frame(biol@n[,ac(1990:2013)])$data,na.rm=T)))[1]))
    if(!iYr %in% c(seq(1990,2007,6),2008:2013)) boxplot(data ~ as.factor(age),data=as.data.frame(biol@n[,ac(iYr)]),xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,rev(pretty(max(as.data.frame(biol@n[,ac(1990:2013)])$data,na.rm=T)))[1]))
    text(x=9,y=rev(pretty(max(as.data.frame(biol@n[,ac(1990:2013)])$data,na.rm=T)))[1],pos=1,labels=iYr,font=2,cex=1)
}
savePlot(file=paste(plotPath,"input_natage.png",sep=""),type="png"); dev.off()


  #- Weight at age
boxplot(as.data.frame(biol@wt[,projPeriod])$data~as.factor(as.data.frame(biol@wt[,projPeriod])$age),
        boxwex=0.5,col="grey90",medlty=0,medlwd=1,medpch=19,medcol="black",boxlty=1,outlty=0,outcex=0.5,outpch=1,staplelty=1,
        xlab="Age",ylab="Weight in the stock(kg)",cex.axis=ca,cex.lab=cl,font=fonts)
points(as.data.frame(iter(Mac@stock.wt[,ac(2001:2012)],1))$data~as.factor(as.data.frame(iter(Mac@stock.wt[,ac(2001:2012)],1))$age),pch=19,cex=0.5,col="red")
legend("bottomright",legend=c("Simulated weights","Observed weights"),col=c("black","red"),
       pch=c(22,19),
       pt.bg="grey90",pt.cex=c(1,1),box.lty=0)
savePlot(file=paste(plotPath,"input_wtatage.png",sep=""),type="png"); dev.off()

xyplot(data ~ year | as.factor(age),data=subset(subset(as.data.frame(biol@wt[ac(1:6)]),year %in% 2001:futureMaxYr),iter %in% 1:5),type="l",group=iter,
       prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
       scales=list(alternating=1,y=list(relation="free",rot=0)))
savePlot(file=paste(plotPath,"input_wtscenarios.png",sep=""),type="png"); dev.off()

  #- Maturity at age
boxplot(as.data.frame(biol@fec[,projPeriod])$data~as.factor(as.data.frame(biol@fec[,projPeriod])$age),
        boxwex=0.5,col="grey90",medlty=0,medlwd=1,medpch=19,medcol="black",boxlty=1,outlty=0,outcex=0.5,outpch=1,staplelty=1,
        xlab="Age",ylab="Proportion mature",cex.axis=ca,cex.lab=cl,font=fonts)
points(as.data.frame(iter(Mac@mat[,ac(2001:2010)],1))$data~as.factor(as.data.frame(iter(Mac@mat[,ac(2001:2010)],1))$age),pch=19,cex=0.5,col="red")
legend("bottomright",legend=c("Simulated maturity","Observed maturity"),col=c("black","red"),
       pch=c(22,19),
       pt.bg="grey90",pt.cex=c(1,1),box.lty=0)
savePlot(file=paste(plotPath,"input_matatage.png",sep=""),type="png"); dev.off()

xyplot(data ~ year | as.factor(age),data=subset(subset(as.data.frame(biol@fec),year %in% 2001:2022),iter %in% 1:5),type="l",group=iter,
       prepanel=function(...) {list(ylim=range(pretty(list(...)$y)))},
       scales=list(alternating=1,y=list(relation="free",rot=0)))
savePlot(file=paste(plotPath,"input_matscenarios.png",sep=""),type="png"); dev.off()



  #-------------------------------------------------------------------------------
  # 1): Figure of fisheries future landing weights and selectivity
  #-------------------------------------------------------------------------------

  #- Landings weight

for(iFsh in dimnames(fishery@landings.wt)$unit){
  boxplot(as.data.frame(fishery@landings.wt[,projPeriod,iFsh])$data~as.factor(as.data.frame(fishery@landings.wt[,projPeriod,iFsh])$age),
          boxwex=0.5,col="grey90",medlty=0,medlwd=1,medpch=19,medcol="black",boxlty=1,outlty=0,outcex=0.5,outpch=1,staplelty=1,
          xlab="Age",ylab="Weight (kg)",cex.lab=cl,cex.axis=ca,font=fonts,main=paste("Fleet",iFsh))
  points(as.data.frame(iter(fishery@landings.wt[,ac(2001:2010),iFsh],1))$data~as.factor(as.data.frame(iter(fishery@landings.wt[,ac(2001:2010),iFsh],1))$age),pch=19,cex=0.5,col="red")
  legend("bottomright",legend=c("Simulated weights","Observed weights"),col=c("black","red"),
         pch=c(22,19),
         pt.bg="grey90",pt.cex=c(1,1),box.lty=0)
}
savePlot(file=paste(plotPath,"input_fish_wtatage.png",sep=""),type="png"); dev.off()

    #- Selectivity

  boxplot(as.data.frame(landings.sel(fishery)[,projPeriod])$data~as.factor(as.data.frame(landings.sel(fishery)[,projPeriod])$age),
          boxwex=0.5,col="grey90",medlty=0,medlwd=1,medpch=19,medcol="black",boxlty=1,outlty=0,outcex=0.5,outpch=1,staplelty=1,
          xlab="Age",ylab="Selectivity",cex.lab=cl,cex.axis=ca,font=fonts,main=paste("Fleet",iFsh))
  lines(as.data.frame(apply(landings.sel(fishery)[,projPeriod],1,median))$data~as.factor(as.data.frame(apply(landings.sel(fishery)[,projPeriod],1,median))$age),
        lwd=2,lty=2)

  lines(as.data.frame(iter(landings.sel(fishery)[,ac(histMaxYr)],1))$data~as.factor(as.data.frame(iter(landings.sel(fishery)[,ac(histMaxYr)],1))$age),
        col="red",lty=1,lwd=2)

savePlot(file=paste(plotPath,"input_fish_selatage.png",sep=""),type="png"); dev.off()

xyplot(data ~ age | as.factor(year),data=subset(subset(as.data.frame(landings.sel(fishery)[,,1]),year %in% 2001:futureMaxYr),iter %in% 1:5),type="l",group=iter,
       prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
       scales=list(alternating=1,y=list(relation="free",rot=0)))
savePlot(file=paste(plotPath,"input_selAscenarios.png",sep=""),type="png"); dev.off()


sel<-sweep(landings.sel(fishery)[ac(0:7),,,,,2:5],c(2:6),landings.sel(fishery)[ac(7),,,,,2:5],"/")
xyplot(data ~ year | iter,group= as.factor(age),data=sel[,ac(2000:2032)],type="l",auto.key=list(columns = 2))
savePlot(file=paste(plotPath,"input_selAscenarios2.png",sep=""),type="png"); dev.off()
  #-------------------------------------------------------------------------------
  # 2): Figure of survey residual pattern
  #-------------------------------------------------------------------------------


xyplot(data ~ year,subset(subset(as.data.frame(Rindex[1,]),year %in% 1998:2012),iter %in% 1:5),type="l",group=iter,
       prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
       scales=list(alternating=1,y=list(relation="free",rot=0)))
savePlot(file=paste(plotPath,"input_Rindex_errors.png",sep=""),type="png"); dev.off()



  #-------------------------------------------------------------------------------
  # 3): Figures of the stock
  #-------------------------------------------------------------------------------
units(stocks)$harvest<-"f"
units(stockstore)$harvest<-"f"
plot(trim(stocks,year=an(1990:histMaxYr)))
savePlot(file=paste(plotPath,"input_stocks.png",sep=""),type="png"); dev.off()


#
##- Assessment stock.n error
#xyplot(data ~ year|as.factor(age),data=as.data.frame(devN),pch=19,col="black",cex=0.4,xlab="Years",ylab="Error (log ratio)",main="Assessment numbers-at-age error")
#savePlot(file=paste(plotPath,"input_stock_errornatage.png",sep=""),type="png"); dev.off()
#
##- Assessment harvest error
#xyplot(data ~ year|as.factor(age),data=as.data.frame(devF),pch=19,col="black",cex=0.4,xlab="Years",ylab="Error (log ratio)",main="Assessment f-at-age error")
#savePlot(file=paste(plotPath,"input_stock_errorfatage.png",sep=""),type="png"); dev.off()
#
##- Catch number residuals
#xyplot(data ~ year|as.factor(age),data=as.data.frame(catchResids),pch=19,col="black",cex=0.4,xlab="Years",ylab="Error multiplier",main="Catch numbers-at-age error")
#savePlot(file=paste(plotPath,"input_stock_errorcatchatage.png",sep=""),type="png"); dev.off()
#




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
plotPath2<-paste(plotPathsc,"/Plots results/",sep="")
dir.create(plotPath2)

#-------------------------------------------------------------------------------
# assessment errors figures
#-------------------------------------------------------------------------------

# error on real scale
devSSB<-(ssb(stockstore)[,projPeriod]-ssbb(biol,f,stockstore)[,projPeriod])
devSSB<-devSSB[,ac(projPeriod)]

histogram(1~data, data=(devSSB),col="grey",main="assess errors SSB(TaY)")
savePlot(file=paste(plotPath2,"1.SSBinTaY_error.png",sep=""),type="png"); dev.off()

par(mfrow=c(3,3),mar=c(3,3,3,3),oma=c(3.1,3.1,1,1))
for (i in 1:9)
{xrange  <- range(an(projPeriod))
yrange  <- range(c(-1.5e6,1.5e6))
plot(c(iter(devSSB,i))~an(projPeriod),type="l",xlab="Years",ylab="",
     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1,col=i)
mtext(side=2,at=yrange[2]/2,text="error on SSB (TaY)",outer=F,line=4,font=fonts,cex=ca)
abline(h=0,lwd=2)
grid(); box()
}
savePlot(file=paste(plotPath2,"1.2.SSBinTaY_error_iters.png",sep=""),type="png"); dev.off()



# error on log scale
devSSB<-log(ssb(stockstore)[,projPeriod]/ssbb(biol,f,stockstore)[,projPeriod])
devSSB<-devSSB[,1:(dim(devSSB[,projPeriod])[2]-1)]


cvSSBTaY<-apply(devSSB,c(1,3:6),sd)
histogram(1~data, data=(devSSB),col="grey",main="assess errors logSSB(TaY)")
savePlot(file=paste(plotPath2,"2.logSSBinTaY_error_CVis",round(mean(cvSSBTaY),2),".png",sep=""),type="png"); dev.off()

 nits<-1000
 
# autocor error on log scale
rhoSSBdev<-rep(NA,nits)
for (i in 1:nits) rhoSSBdev[i]<- an(unlist(acf(iter(devSSB,i),plot=F)[1])[1])
histogram(1~rhoSSBdev,col="grey",main="RHO errors logSSB(TaY)")
savePlot(file=paste(plotPath2,"3.RHO_logSSBinTaYerrors_mean.is",round(mean(rhoSSBdev),2),".png",sep=""),type="png"); dev.off()


# error on the advisory year
ys<-dimnames(SSB)$year[which(!is.na(iter(SSB,1)))]
devSSB<-log(SSB[,ys]/ssbb(biol,f,stocks)[,ys])
devSSB<-devSSB[,1:(dim(devSSB)[2]-1)]
cvSSBAdY<-apply(devSSB,c(1,3:6),sd)
histogram(1~data, data=(devSSB),col="grey",main="assess errors logSSB(AdY)")
savePlot(file=paste(plotPath2,"4.logSSBinAdY_error_CVis",round(mean(cvSSBAdY),2),".png",sep=""),type="png"); dev.off()

# autocor error in advice year on log scale
rhoSSBdev<-rep(NA,nits)
for (i in 1:nits) rhoSSBdev[i]<- an(unlist(acf(iter(devSSB,i),plot=F)[1])[1])
histogram(1~rhoSSBdev,col="grey",main="RHO errors logSSB(AdY)")
savePlot(file=paste(plotPath2,"5.RHO_logSSBinAdYerrors_mean.is",round(mean(rhoSSBdev),2),".png",sep=""),type="png"); dev.off()






devFbar<-log(quantMeans(harvest(stockstore)[ac(4:8),])/quantMeans(f[ac(4:8),]))
devFbar<-devFbar[,1:(dim(devFbar)[2]-1)]
devFbar<-devFbar[,ac(projPeriod)]
cvFbarTaY<-apply(devFbar,c(1,3:6),sd)
rhoFdev<-rep(NA,nits)
for (i in 1:nits) rhoFdev[i]<- an(unlist(acf(iter(devFbar,i),plot=F)[1])[1])

#
#xyplot(data ~ year ,data=subset(subset(as.data.frame(devFbar),year %in% 2012:2032),iter %in% 1:5),type="l",group=iter,
#       prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
#       scales=list(alternating=1,y=list(relation="free",rot=0)))
##savePlot(file=paste(plotPath,"input_Fbar_errorscenario.png",sep=""),type="png"); dev.off()



histogram(1~data, data=(devFbar),col="grey",main="assess errors logFbar(TaY)")
savePlot(file=paste(plotPath2,"6.logFbarinTaY_error_CVis",round(mean(cvFbarTaY),2),".png",sep=""),type="png"); dev.off()

histogram(1~rhoFdev,col="grey",main="RHO errors logFbar(TaY)")
savePlot(file=paste(plotPath2,"7.RHO_logFbarinTaYerrors_mean.is",round(mean(rhoFdev),2),".png",sep=""),type="png"); dev.off()









#-------------------------------------------------------------------------------
# stock trajectories figures
#-------------------------------------------------------------------------------



############### 

rSSBp <- apply(ssbb(biol,unitSums(f),stockstore)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rSSBs <- apply(ssb(stockstore)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rFp   <- apply(apply(unitSums(f)[ac(4:8),],2:6,mean,na.rm=T)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rFs   <- apply(fbar(stockstore)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)

rLandf<- apply(computeLandings(fishery)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rLands<- apply(computeLandings(stockstore)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)


  #- Plot settings
yrs   <- 1990:max(an(projPeriod))
cl    <- 1.2
ca    <- 1.1
fonts <- 2

 
 
#-------------------------------------------------------------------------------
# 3): Plot results of SSB stock and SSB pop
#-------------------------------------------------------------------------------
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
savePlot(file=paste(plotPath2,"8.SSB.png",sep=""),type="png");dev.off()








#-------------------------------------------------------------------------------
# 4): Plot results of F stock and true F
#-------------------------------------------------------------------------------
rFp   <- apply(apply(unitSums(f)[ac(4:8),],2:6,mean,na.rm=T)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rFs   <- apply(fbar(stockstore)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)



xrange  <- range(yrs)
yrange  <- c(0,range(pretty(c(rFp[,,ac(yrs),,,],rFs[,,ac(yrs),,,])),na.rm=T)[2])

  #---------
  #- Biology
  #---------
plot(rFp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="Years",ylab="",
     type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
mtext(text="Fishing mortality (ages 4-8)",side=2,at=0.5,outer=T,cex=cl,line=-1)
grid(); box()
#-Reference level Fpa
#abline(h=0.26,col="darkgreen",lwd=2,lty=2);
#mtext(text="Fpa",side=4,at=0.26,las=1,cex=0.8,col="darkgreen",line=0.5,font=fonts)
#lines(rFp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rFp["5%",,ac(yrs),,,]~yrs,lty=3)
lines(rFp["95%",,ac(yrs),,,]~yrs,lty=3)

  #---------
  #- Stock
  #---------
rFs["50%",,ac(2015:2046),,,]<-rFs["50%",,ac(2014:2045),,,]
rFs["5%",,ac(2015:2046),,,]<-rFs["5%",,ac(2014:2045),,,]
rFs["95%",,ac(2015:2046),,,]<-rFs["95%",,ac(2014:2045),,,]

lines(rFs["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
lines(rFs["5%",,ac(yrs),,,]~yrs,lty=3,col="red")
lines(rFs["95%",,ac(yrs),,,]~yrs,lty=3,col="red")

  #---------
  #- Legend
  #---------

legend("bottomright",legend=c("F assessed stock","F true population"),
       col=c("red","black"),lwd=3,lty=1,box.lty=0)
savePlot(file=paste(plotPath2,"9.F.png",sep=""),type="png");dev.off()










#-------------------------------------------------------------------------------
# 5): Plot results of landings by fleet & TAC on top
#-------------------------------------------------------------------------------


  #---------
  #- Fishery
  #---------
  par(mar=c(5.1,5.1,4.1,2.1))
iFsh<-1
  xrange  <- range(yrs)
  yrange  <- c(0,range(pretty(c(rLandf[,,ac(yrs),iFsh,,])),na.rm=T)[2])

    #- Landings
  plot(rLandf["50%",,ac(yrs),iFsh,,]~yrs,xlim=xrange,ylim=yrange,xlab="Year",ylab="Landings (thousand tonnes)",
       type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
  axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
  grid(); box()
  lines(rLandf["50%",,ac(yrs),iFsh,,]~yrs,lty=1,lwd=2)
  lines(rLandf["5%",,ac(yrs),iFsh,,]~yrs,lty=3)
  lines(rLandf["95%",,ac(yrs),iFsh,,]~yrs,lty=3)
    #- TAC
 

legend("bottomright",legend=c("Landings","Advice TAC"),
       col=c("black","red"),lwd=3,lty=1,box.lty=0)
savePlot(file=paste(plotPath2,"10.CatchTAC.png",sep=""),type="png");dev.off()





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
savePlot(file=paste(plotPath2,"11.Catch.png",sep=""),type="png");dev.off()









#-------------------------------------------------------------------------------
# 9): Report figures
#-------------------------------------------------------------------------------

yrs   <- 1990:an(rev(projPeriod)[2])
cl    <- 1.1
ca    <- 1
fonts <- 1
yrangeSSB <- c(0,7e6)
yrangeLan <- c(0,2e6)
#yrangeLan2<- c(0,2.5e4)
yrangeF   <- c(0,0.7)


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
text(x=xrange[1],y=yrange[2],pos=1,labels="(A)",font=fonts,cex=cl)
abline(v=2013.5,col="green",lwd=2)
itr<-sample(1:1000,2)
ct<-computeLandings(fishery)[,ac(yrs)]
lines(yrs,c(iter(ct,itr[1])),col="red")
lines(yrs,c(iter(ct,itr[2])),col="blue")
lines(rLandf["50%",,ac(yrs),1,,]~yrs,lty=1,lwd=2)
  #------------------
  # True F
  #------------------

xrange  <- range(yrs)
yrange  <- yrangeF
plot(rFs["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
mtext(text=expression(paste(F[2-6]," (",year^-1,")",sep="")),side=2,at=(yrange[2]-yrange[1])/2+yrange[1],outer=F,cex=cl,line=4,font=fonts)
grid(); box()
lines(rFs["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rFs["5%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
lines(rFs["95%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
fb<-quantMeans(f[ac(4:8),ac(yrs)])
lines(yrs,c(iter(fb,itr[1])),col="red")
lines(yrs,c(iter(fb,itr[2])),col="blue")
lines(rFs["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
  #-Reference level Fpa
#abline(h=0.26,col="darkgreen",lwd=1,lty=2);
#mtext(text="Fpa",side=4,at=0.26,las=1,cex=0.65,col="darkgreen",line=0.5,font=fonts)
text(x=xrange[1],y=yrange[2],pos=1,labels="(B)",font=fonts,cex=cl)
abline(v=2013.5,col="green",lwd=2)
  #------------------
  # True SSB
  #------------------

xrange  <- range(yrs)
yrange  <- yrangeSSB
plot(rSSBs["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
axis(1,las=1,at=pretty(xrange),labels=pretty(xrange),cex=ca,font=fonts)
mtext(text=expression(paste("SSB (",10^3," tonnes)",sep="")),side=2,at=(yrange[2]-yrange[1])/2+yrange[1],outer=F,cex=cl,line=4,font=fonts)
grid(); box()
lines(rSSBs["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rSSBs["5%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
lines(rSSBs["95%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
sb<-ssbb(biol,unitSums(f),stockstore)[,ac(yrs)]
lines(yrs,c(iter(sb,itr[1])),col="red")
lines(yrs,c(iter(sb,itr[2])),col="blue")
lines(rSSBs["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)

  #-Reference level Blim & Bpa
abline(h=1.84e6,col="blue",lwd=1,lty=2);
abline(h=2.36e6,col="darkgreen",lwd=1,lty=2);
mtext(text="Blim",side=4,at=1.84e6,las=1,cex=0.65,col="blue",line=0.5,font=fonts)
mtext(text="Bpa",side=4,at=2.36e6,las=1,cex=0.65,col="darkgreen",line=0.5,font=fonts)
text(x=xrange[1],y=yrange[2],pos=1,labels="(C)",font=fonts,cex=cl)
abline(v=2013.5,col="green",lwd=2)


  #- Labels x-axis
mtext(text=expression(Years),side=1,at=(xrange[2]-xrange[1])/2+xrange[1],outer=F,cex=cl,line=4,font=fonts)
savePlot(paste(plotPath2,"17.Truth.png",sep=""),type="png") ; dev.off()

