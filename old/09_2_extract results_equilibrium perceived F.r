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
#library(FLAssess)
#library(FLSAM)
library(MASS)
#library(msm)

wine <- F

path          <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/"
inPath        <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/Data/"
codePath      <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/R code/"
outPath       <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/Results/"
if(substr(R.Version()$os,1,3)== "lin"){
  path        <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",path)
  inPath      <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",inPath)
  codePath    <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",codePath)
  outPath     <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",outPath)
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


source(paste(codePath,"functions.r",            sep=""))
load(file=paste(outPathp,"stocks.RData",         sep=""))
 outPath2      <- paste(outPath,"equilibrium/perceived/",sep="")
 scen          <- "equilPerc"              # 
 opt           <- ""                      # for multiple scenario combinations a counter


SSBsT<-c()
SSBsP<-c()
RecsT<-c()
RecsP<-c()
Yields<-c()
YieldsP<-c()
trueF<-c()
perceivedF<-c()
RiskBpa<-c()
RiskBlim<-c()
              
Risk1ShortT  <-c() 
Risk1MidT    <-c() 
Risk1LongT   <-c() 
              
Risk2ShortT  <-c() 
Risk2MidT    <-c() 
Risk2LongT   <-c() 
              
Risk3ShortT  <-c() 
Risk3MidT    <-c() 
Risk3LongT   <-c() 
              







Bpa   <-  2.36e6
Blim  <-  1.84e6


# choice of the set up
nits        <-  200      # number of iterations
prlength    <-  70       # length of the projection period
nYrMean     <-  40       # nb of years over which to take the mean
prlength2    <- 70       # length of the projection period  used for MSY calculation


############### plot the run for one of the F values to check that we are really at equ
checkequ<-F
if(checkequ==T)
{
fequ<-0.25
 mpPoints      <- list(Fequ=fequ)
load(file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalbiol.RData",        sep=""))
load(file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalfishery.RData",     sep=""))
load(file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalpercievedstocks.RData",      sep=""))
load(file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalf.RData",           sep=""))
projPeriod <- ac(2014:(2013+prlength2) )
#

png(file=paste(outPath2,"check equilibrium.png",sep=""),width=8,height=8,units="in",res=200)
plot(ssbb(biol,f,stockstore)[,ac(2000:rev(projPeriod)[1])],ylab="SSB")
dev.off()
}






############### compile the results for all F values
fs<- seq(0.0,0.5,0.01)

for (fequ in fs)
{          

cat(fequ,"\n")    
 mpPoints      <- list(Fequ=fequ)
  
load(file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalbiol.RData",        sep=""))
load(file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalfishery.RData",     sep=""))
load(file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalpercievedstocks.RData",      sep=""))
load(file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalf.RData",           sep=""))

nits<-dim(f)[6]
projPeriod <- ac(2014:(2013+prlength2) )
if (length(projPeriod)<20) cat("! these results are projecting only 20 years head !","\n")    
stocks<-stockstore[,,,,,nits]

stocks<-trim(stocks,year=an(projPeriod))
stockstore<-trim(stockstore,year=an(projPeriod))
biol<-trim(biol,year=an(projPeriod))
fishery<-trim(fishery,year=an(projPeriod))
f<-trim(f,year=an(projPeriod))


tp<-ssbb(biol,f,stocks)
tp<-c(yearMeans(tp[,(dim(tp)[2]-nYrMean):dim(tp)[2]])@.Data)
SSBsT<-c(SSBsT,tp)


tp<-ssbb(biol,f,stocks)
#tp<-tp<Blim
#tp<-yearSums(tp)
#tp<-tp>=1
#RiskBlim<-c(RiskBlim,sum(tp)/nits)
#cat("Blim",length(RiskBlim),"\n")
#
ShortT    <-  ac(2014:2018)
MidT      <-  ac(2019:2028)
LongT     <-  ac(2045:2085)


risk<-  apply(tp<Blim,c(1:5),sum)/nits

Risk1ShortT             <- c(Risk1ShortT,mean(risk[,ShortT])*100)       # percentage of iteration that reach Blim
Risk1MidT               <- c(Risk1MidT,mean(risk[,MidT])  *100)
Risk1LongT              <- c(Risk1LongT,mean(risk[,LongT]) *100)

Risk2ShortT             <- c(Risk2ShortT,length(unique(which(tp[,ShortT]<Blim,arr.ind=T)[,"dim6"]))/dims(tp)$iter  *100  )     # percentage of iteration that reach Blim
Risk2MidT               <- c(Risk2MidT,length(unique(which(tp[,MidT]  <Blim,arr.ind=T)[,"dim6"]))/dims(tp)$iter  *100  )
Risk2LongT              <- c(Risk2LongT,length(unique(which(tp[,LongT] <Blim,arr.ind=T)[,"dim6"]))/dims(tp)$iter  *100  )

Risk3ShortT             <- c(Risk3ShortT,max(risk[,ShortT])*100)       # percentage of iteration that reach Blim
Risk3MidT               <- c(Risk3MidT,max(risk[,MidT])  *100)
Risk3LongT              <- c(Risk3LongT,max(risk[,LongT]) *100)



tp<-ssb(stockstore)
tp<-c(yearMeans(tp[,(dim(tp)[2]-nYrMean):dim(tp)[2]])@.Data)
SSBsP<-c(SSBsP,tp)

tp<-n(biol)[ac(1),]
tp<-c(yearMeans(tp[,(dim(tp)[2]-nYrMean):dim(tp)[2]])@.Data)
RecsT<-c(RecsT,tp)

tp<-stock.n(stockstore)[ac(1),]
tp<-c(yearMeans(tp[,(dim(tp)[2]-nYrMean):dim(tp)[2]])@.Data)
RecsP<-c(RecsP,tp)


tp<-quantSums(fishery@landings.n*fishery@landings.wt)
tp<-c(yearMeans(tp[,(dim(tp)[2]-nYrMean):dim(tp)[2]])@.Data)
Yields<-c(Yields,tp) 

tp<-quantSums (computeLandings(stockstore))
tp<-c(yearMeans(tp[,(dim(tp)[2]-nYrMean):dim(tp)[2]])@.Data)
YieldsP<-c(YieldsP,tp) 

tp<-quantMeans(f[ac(4:8),])
tp<-c(yearMeans(tp[,(dim(tp)[2]-nYrMean):dim(tp)[2]])@.Data)
trueF<-c(trueF,tp)

tp<-quantMeans(stockstore@harvest[ac(4:8),])
tp<-c(yearMeans(tp[,(dim(tp)[2]-nYrMean):dim(tp)[2]])@.Data)
perceivedF<-c(perceivedF,tp)


save.image(file=paste(outPath2,scen,opt,nits,"its",prlength2,"yrs_resMSY.RData",        sep=""))
}




##################################################################################################################################################################
##################
##################
##################
##################
##################################################################################################################################################################








load(file=paste(outPath2,scen,opt,nits,"its",prlength2,"yrs_resMSY.RData",        sep=""))

outPath       <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/Results/"
 outPath2      <- paste(outPath,"equilibrium/perceived/",sep="")


nits<-dim(f)[6]

iter<-rep(1:nits,length(fs))

res<-data.frame(trueF=trueF,perceivedF=perceivedF,SSBsT=SSBsT,SSBsP=SSBsP,Yields=Yields,RecT=RecsT,RecP=RecsP,iter=iter)
res$perceivedF<-round(res$perceivedF,2)

res$perceivedF[res$trueF==0]<-0
res$SSBsP[res$trueF==0]<-res$SSBsT[res$trueF==0]

#
## plot as done by Einar
#png(file=paste(outPath2,"plot equilibrium",prlength2,"_",nYrMean,".png",sep=""),width=8,height=8,units="in",res=200)
#par(mfrow=c(2,2),mar=c(4,4,4,4))
## plot for the recs
#q95<-aggregate(RecT~perceivedF,data=res,function(x) {quantile(x,0.95)})
#q75<-aggregate(RecT~perceivedF,data=res,function(x) {quantile(x,0.75)})
#q50<-aggregate(RecT~perceivedF,data=res,function(x) {quantile(x,0.50)})
#q25<-aggregate(RecT~perceivedF,data=res,function(x) {quantile(x,0.25)})
#q05<-aggregate(RecT~perceivedF,data=res,function(x) {quantile(x,0.05)})
#
#plot(q50$perceivedF,q50$RecT,type="l",col="red", lwd=2, xlab="perceived Fbar",ylab="Rec",main="Recruitment",ylim=range(c(q05,q95)))
#polygon(c(q25$perceivedF,rev(q75$perceivedF)),c(q25$Rec,rev(q75$RecT)),col=rgb(1,0,0,0.6),border=F)
#polygon(c(q05$perceivedF,rev(q95$perceivedF)),c(q05$Rec,rev(q95$RecT)),col=rgb(1,0,0,0.3),border=F)
#lines(q50$perceivedF,q50$RecT,lwd=2)
#
#q50<-aggregate(RecP~perceivedF,data=res[res$perceivedF<=0.6,],function(x) {quantile(x,0.50)})
#lines(q50$perceivedF,q50$RecP,lwd=2,lty="dotted")
#legend(x="topright",c("true","perceived"),lty=c("solid","dotted"),bty="n")
#
#
#
#
## plot for the Fbars
#q95<-aggregate(trueF~perceivedF,data=res,function(x) {quantile(x,0.95)})
#q75<-aggregate(trueF~perceivedF,data=res,function(x) {quantile(x,0.75)})
#q50<-aggregate(trueF~perceivedF,data=res,function(x) {quantile(x,0.50)})
#q25<-aggregate(trueF~perceivedF,data=res,function(x) {quantile(x,0.25)})
#q05<-aggregate(trueF~perceivedF,data=res,function(x) {quantile(x,0.05)})
#
#plot(q50$perceivedF,q50$trueF,type="l",col="red", lwd=2, xlab="perceived Fbar",ylab="true Fbar",main="Fishing Mortality",ylim=range(c(q05,q95)))
#polygon(c(q25$perceivedF,rev(q75$perceivedF)),c(q25$trueF,rev(q75$trueF)),col=rgb(1,0,0,0.6),border=F)
#polygon(c(q05$perceivedF,rev(q95$perceivedF)),c(q05$trueF,rev(q95$trueF)),col=rgb(1,0,0,0.3),border=F)
#lines(q50$perceivedF,q50$trueF,lwd=2)
#
#
#
## plot for the Yields
#q95<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.95)})
#q75<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.75)})
#q50<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.50)})
#q25<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.25)})
#q05<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.05)})
#
#
#xrange  <- range(c(0,q50$perceivedF))
#yrange  <- c(0,range(pretty(c(q95$Yields)),na.rm=T)[2])
#
#
#plot(q50$perceivedF,q50$Yields,bty="n",type="l",col="red", lwd=2, xaxt="n",xlab="perceived Fbar",ylab="Yields",main="Yields",yaxt="n",xlim=xrange,ylim=yrange)
#
#axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000)
#axis(4,las=1,at=pretty(pretty(yrange)/max(pretty(yrange)))*max(pretty(yrange)),labels=pretty(pretty(yrange)/max(pretty(yrange))))
#axis(1,las=1,at=pretty(xrange),labels=pretty(xrange),pos=0)
#
#polygon(c(q25$perceivedF,rev(q75$perceivedF)),c(q25$Yields,rev(q75$Yields)),col=rgb(1,0,0,0.6),border=F)
#polygon(c(q05$perceivedF,rev(q95$perceivedF)),c(q05$Yields,rev(q95$Yields)),col=rgb(1,0,0,0.3),border=F)
#
##q50b<-aggregate(YieldsP~perceivedF,data=res,function(x) {quantile(x,0.50)})
##lines(q50b$perceivedF,q50b$YieldsP,lwd=2,lty="dotted",bty="n")
##legend(x="topright",c("true","perceived"),lty=c("solid","dotted"),bty="n")
##
#
##  MSY
#Fmsy<-q50$perceivedF[rev(order(q50$Yields))[1] ]
#
#risk<-data.frame(perceivedF=c(0,q50$perceivedF),RiskBpa,RiskBlim)
##lines(risk$perceivedF,risk$RiskBpa*max(q95$Yields),lty="dotted",col="green")
#risk$dif<-abs(risk$RiskBlim-0.05)
#Fmsyprec<-risk$perceivedF[order(risk$dif)[1]]
#MSYprec<-q50$Yields[q50$perceivedF==Fmsyprec]
#
#lines(q50$perceivedF,q50$Yields,lwd=2)
#segments(Fmsy,-2,Fmsy,max(q50$Yields))
#segments(Fmsyprec,-2,Fmsyprec,MSYprec,col="green")
#lines(risk$perceivedF[risk$RiskBlim>0],risk$RiskBlim[risk$RiskBlim>0]*max(pretty(yrange)),col="green",lwd=2)
#abline(h=0.05*max(pretty(yrange)),col="green")
#text(0.07,0.05*max(pretty(yrange)),"5% risk",pos=3,col="green")
#text(Fmsy,max(q50$Yields),Fmsy,pos=3)
#text(Fmsyprec,MSYprec,Fmsyprec,pos=3,col="green")
#abline(h=0)
#
## plot for the SSBs
#q95<-aggregate(SSBsT~perceivedF,data=res,function(x) {quantile(x,0.95)})
#q75<-aggregate(SSBsT~perceivedF,data=res,function(x) {quantile(x,0.75)})
#q50<-aggregate(SSBsT~perceivedF,data=res,function(x) {quantile(x,0.50)})
#q25<-aggregate(SSBsT~perceivedF,data=res,function(x) {quantile(x,0.25)})
#q05<-aggregate(SSBsT~perceivedF,data=res,function(x) {quantile(x,0.05)})
#
#plot(q50$perceivedF,q50$SSBsT,type="l",col="red", lwd=2, xlab="perceived Fbar",ylab="SSBs",main="Spawning Stock Biomass",ylim=range(c(q05,q95)))
#polygon(c(q25$perceivedF,rev(q75$perceivedF)),c(q25$SSBs,rev(q75$SSBsT)),col=rgb(1,0,0,0.6),border=F)
#polygon(c(q05$perceivedF,rev(q95$perceivedF)),c(q05$SSBs,rev(q95$SSBsT)),col=rgb(1,0,0,0.3),border=F)
#lines(q50$perceivedF,q50$SSBsT,lwd=2)
#
#q50$diffBpa<-abs(q50$SSBsT-Bpa)
#q50$diffBlim<-abs(q50$SSBsT-Blim)
#
#Fpa  <-  q50$perceivedF[order(q50$diffBpa)[1]]
#Flim <-  q50$perceivedF[order(q50$diffBlim)[1]]
#
#segments(0,Bpa,Fpa,Bpa)
#segments(0,Blim,Flim,Blim)
#segments(Fpa,Bpa,Fpa,0)
#segments(Flim,Blim,Flim,0)
#text(Fpa,Bpa,paste("Fpa=",Fpa),pos=3)
#text(Flim,Blim,paste("Flim=",Flim),pos=3)
#
#MSYBtrigger<-round(q05$SSBsT[q05$perceivedF==Fmsyprec],0)
#segments(0,MSYBtrigger,Fmsyprec,MSYBtrigger)
#segments(Fmsyprec,MSYBtrigger,Fmsyprec,0)
#text(Fmsyprec,MSYBtrigger,paste("MSYBtrigger=",MSYBtrigger),pos=3)
#
#q50<-aggregate(SSBsP~perceivedF,data=res[res$perceivedF<=0.6,],function(x) {quantile(x,0.50)})
#lines(q50$perceivedF,q50$SSBsP,lwd=2,lty="dotted",bty="n")
#legend(x="topright",c("true","perceived"),lty=c("solid","dotted"),bty="n")
#
#
#dev.off()
#








png(file=paste(outPath2,"plotsBaseCase/plot_equilibrium_and_risk",prlength2,"_",nYrMean,".png",sep=""),width=12,height=4,units="in",res=200)
par(mfrow=c(1,3),mar=c(4,4,4,4))


risk<-data.frame(perceivedF=c(q50$perceivedF),
                          Risk1ShortT,Risk1MidT,Risk1LongT,
                          Risk2ShortT,Risk2MidT,Risk2LongT,
                          Risk3ShortT,Risk3MidT,Risk3LongT)



# plot for the Yields
q95<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.95)})
q75<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.75)})
q50<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.50)})
q25<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.25)})
q05<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.05)})


xrange  <- range(c(0,q50$perceivedF))
yrange  <- c(0,range(pretty(c(q95$Yields)),na.rm=T)[2])


plot(q50$perceivedF,q50$Yields,bty="n",type="l",col="red", lwd=2,lty=0, xaxt="n",xlab="perceived Fbar",ylab="Yields",main="equ Yields and Risk1",yaxt="n",xlim=xrange,ylim=yrange)

axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000)
axis(4,las=1,at=pretty(pretty(yrange)/max(pretty(yrange)))*max(pretty(yrange)),labels=pretty(pretty(yrange)/max(pretty(yrange))))
axis(1,las=1,at=pretty(xrange),labels=pretty(xrange),pos=0)

polygon(c(q25$perceivedF,rev(q75$perceivedF)),c(q25$Yields,rev(q75$Yields)),col=rgb(1,0,0,0.6),border=F)
polygon(c(q05$perceivedF,rev(q95$perceivedF)),c(q05$Yields,rev(q95$Yields)),col=rgb(1,0,0,0.3),border=F)


YPR<-lowess(q50$Yields~q50$perceivedF,f=1/6)
Fmsy<-YPR$x[rev(order(YPR$y))[1] ]
lines(YPR,lwd=2)
segments(Fmsy,-2,Fmsy,max(YPR$y))
text(Fmsy,1.1*max(q50$Yields),Fmsy,pos=3)
segments(0,max(YPR$y),Fmsy,max(YPR$y))
text(0,max(YPR$y),pos=3,round(max(YPR$y)/1000))

lines(risk$perceivedF,risk$Risk1ShortT*max(pretty(yrange))/100,lty="dotted",col="green",lwd=2)
lines(risk$perceivedF,risk$Risk1MidT  *max(pretty(yrange))/100,lty="dashed",col="green",lwd=2)
lines(risk$perceivedF,risk$Risk1LongT *max(pretty(yrange))/100,lty="solid",col="green",lwd=2)
abline(h=0.05*max(pretty(yrange)),col="green")
text(0.07,0.05*max(pretty(yrange)),"5% risk",pos=3,col="green")

risk$difShortT<-abs(risk$Risk1ShortT-5)
risk$difMidT<-abs(risk$Risk1MidT-5)
risk$difLongT<-abs(risk$Risk1LongT-5)

FmsyprecShortT<-risk$perceivedF[order(risk$difShortT)[1]]
if(FmsyprecShortT==0.5) FmsyprecShortT<-1
FmsyprecMidT<-risk$perceivedF[order(risk$difMidT)[1]]
FmsyprecMidT1<-FmsyprecMidT
FmsyprecLongT<-risk$perceivedF[order(risk$difLongT)[1]]
MSYprecShortT<-q50$Yields[q50$perceivedF==FmsyprecShortT]
 if(FmsyprecShortT==1) MSYprecShortT<-1
MSYprecMidT<-q50$Yields[q50$perceivedF==FmsyprecMidT]
MSYprecLongT<-q50$Yields[q50$perceivedF==FmsyprecLongT]

segments(FmsyprecShortT,-2,FmsyprecShortT,MSYprecShortT,col="green",lty="dotted")
segments(FmsyprecMidT,-2,FmsyprecMidT,MSYprecMidT,col="green",lty="dashed")
segments(FmsyprecLongT,-2,FmsyprecLongT,MSYprecLongT,col="green",lty="solid")
text(FmsyprecShortT,MSYprecShortT,FmsyprecShortT,pos=3,col="green")
text(FmsyprecMidT,MSYprecMidT,FmsyprecMidT,pos=3,col="green")
text(FmsyprecLongT,MSYprecLongT,FmsyprecLongT,pos=3,col="green")
abline(h=0)
legend(x='topright',lty=c("dotted","dashed","solid"),c("ShortT","MidT","LongT"),col="green",bty="n")






# plot for the Yields
q95<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.95)})
q75<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.75)})
q50<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.50)})
q25<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.25)})
q05<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.05)})


xrange  <- range(c(0,q50$perceivedF))
yrange  <- c(0,range(pretty(c(q95$Yields)),na.rm=T)[2])


plot(q50$perceivedF,q50$Yields,bty="n",type="l",col="red", lwd=2,lty=0, xaxt="n",xlab="perceived Fbar",ylab="Yields",main="equ Yields and Risk2",yaxt="n",xlim=xrange,ylim=yrange)

axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000)
axis(4,las=1,at=pretty(pretty(yrange)/max(pretty(yrange)))*max(pretty(yrange)),labels=pretty(pretty(yrange)/max(pretty(yrange))))
axis(1,las=1,at=pretty(xrange),labels=pretty(xrange),pos=0)

polygon(c(q25$perceivedF,rev(q75$perceivedF)),c(q25$Yields,rev(q75$Yields)),col=rgb(1,0,0,0.6),border=F)
polygon(c(q05$perceivedF,rev(q95$perceivedF)),c(q05$Yields,rev(q95$Yields)),col=rgb(1,0,0,0.3),border=F)


YPR<-lowess(q50$Yields~q50$perceivedF,f=1/6)
Fmsy<-YPR$x[rev(order(YPR$y))[1] ]
lines(YPR,lwd=2)
segments(Fmsy,-2,Fmsy,max(YPR$y))
text(Fmsy,1.1*max(q50$Yields),Fmsy,pos=3)
segments(0,max(YPR$y),Fmsy,max(YPR$y))
text(0,max(YPR$y),pos=3,round(max(YPR$y)/1000))




lines(risk$perceivedF,risk$Risk2ShortT*max(pretty(yrange))/100,lty="dotted",col="green",lwd=2)
lines(risk$perceivedF,risk$Risk2MidT  *max(pretty(yrange))/100,lty="dashed",col="green",lwd=2)
lines(risk$perceivedF,risk$Risk2LongT *max(pretty(yrange))/100,lty="solid",col="green",lwd=2)
abline(h=0.05*max(pretty(yrange)),col="green")
text(0.07,0.05*max(pretty(yrange)),"5% risk",pos=3,col="green")

risk$difShortT<-abs(risk$Risk2ShortT-5)
risk$difMidT<-abs(risk$Risk2MidT-5)
risk$difLongT<-abs(risk$Risk2LongT-5)

FmsyprecShortT<-risk$perceivedF[order(risk$difShortT)[1]]
if(FmsyprecShortT==0.5) FmsyprecShortT<-1
FmsyprecMidT<-risk$perceivedF[order(risk$difMidT)[1]]
FmsyprecMidT2<-FmsyprecMidT
FmsyprecLongT<-risk$perceivedF[order(risk$difLongT)[1]]
MSYprecShortT<-q50$Yields[q50$perceivedF==FmsyprecShortT]
 if(FmsyprecShortT==1) MSYprecShortT<-1
MSYprecMidT<-q50$Yields[q50$perceivedF==FmsyprecMidT]
MSYprecLongT<-q50$Yields[q50$perceivedF==FmsyprecLongT]

segments(FmsyprecShortT,-2,FmsyprecShortT,MSYprecShortT,col="green",lty="dotted")
segments(FmsyprecMidT,-2,FmsyprecMidT,MSYprecMidT,col="green",lty="dashed")
segments(FmsyprecLongT,-2,FmsyprecLongT,MSYprecLongT,col="green",lty="solid")
text(FmsyprecShortT,MSYprecShortT,FmsyprecShortT,pos=3,col="green")
text(FmsyprecMidT,MSYprecMidT,FmsyprecMidT,pos=3,col="green")
text(FmsyprecLongT,MSYprecLongT,FmsyprecLongT,pos=3,col="green")
abline(h=0)
legend(x='topright',lty=c("dotted","dashed","solid"),c("ShortT","MidT","LongT"),col="green",bty="n")




# plot for the Yields
q95<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.95)})
q75<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.75)})
q50<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.50)})
q25<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.25)})
q05<-aggregate(Yields~perceivedF,data=res,function(x) {quantile(x,0.05)})


xrange  <- range(c(0,q50$perceivedF))
yrange  <- c(0,range(pretty(c(q95$Yields)),na.rm=T)[2])


plot(q50$perceivedF,q50$Yields,bty="n",type="l",col="red", lwd=2,lty=0, xaxt="n",xlab="perceived Fbar",ylab="Yields",main="equ Yields and Risk3",yaxt="n",xlim=xrange,ylim=yrange)

axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000)
axis(4,las=1,at=pretty(pretty(yrange)/max(pretty(yrange)))*max(pretty(yrange)),labels=pretty(pretty(yrange)/max(pretty(yrange))))
axis(1,las=1,at=pretty(xrange),labels=pretty(xrange),pos=0)

polygon(c(q25$perceivedF,rev(q75$perceivedF)),c(q25$Yields,rev(q75$Yields)),col=rgb(1,0,0,0.6),border=F)
polygon(c(q05$perceivedF,rev(q95$perceivedF)),c(q05$Yields,rev(q95$Yields)),col=rgb(1,0,0,0.3),border=F)


YPR<-lowess(q50$Yields~q50$perceivedF,f=1/6)
Fmsy<-YPR$x[rev(order(YPR$y))[1] ]
lines(YPR,lwd=2)
segments(Fmsy,-2,Fmsy,max(YPR$y))
text(Fmsy,1.1*max(q50$Yields),Fmsy,pos=3)
segments(0,max(YPR$y),Fmsy,max(YPR$y))
text(0,max(YPR$y),pos=3,round(max(YPR$y)/1000))




lines(risk$perceivedF,risk$Risk3ShortT*max(pretty(yrange))/100,lty="dotted",col="green",lwd=2)
lines(risk$perceivedF,risk$Risk3MidT  *max(pretty(yrange))/100,lty="dashed",col="green",lwd=2)
lines(risk$perceivedF,risk$Risk3LongT *max(pretty(yrange))/100,lty="solid",col="green",lwd=2)
abline(h=0.05*max(pretty(yrange)),col="green")
text(0.07,0.05*max(pretty(yrange)),"5% risk",pos=3,col="green")

risk$difShortT<-abs(risk$Risk3ShortT-5)
risk$difMidT<-abs(risk$Risk3MidT-5)
risk$difLongT<-abs(risk$Risk3LongT-5)

FmsyprecShortT<-risk$perceivedF[order(risk$difShortT)[1]]
if(FmsyprecShortT==0.5) FmsyprecShortT<-1
FmsyprecMidT<-risk$perceivedF[order(risk$difMidT)[1]]
FmsyprecMidT3<-FmsyprecMidT
FmsyprecLongT<-risk$perceivedF[order(risk$difLongT)[1]]
MSYprecShortT<-q50$Yields[q50$perceivedF==FmsyprecShortT]
 if(FmsyprecShortT==1) MSYprecShortT<-1
MSYprecMidT<-q50$Yields[q50$perceivedF==FmsyprecMidT]
MSYprecLongT<-q50$Yields[q50$perceivedF==FmsyprecLongT]

segments(FmsyprecShortT,-2,FmsyprecShortT,MSYprecShortT,col="green",lty="dotted")
segments(FmsyprecMidT,-2,FmsyprecMidT,MSYprecMidT,col="green",lty="dashed")
segments(FmsyprecLongT,-2,FmsyprecLongT,MSYprecLongT,col="green",lty="solid")
text(FmsyprecShortT,MSYprecShortT,FmsyprecShortT,pos=3,col="green")
text(FmsyprecMidT,MSYprecMidT,FmsyprecMidT,pos=3,col="green")
text(FmsyprecLongT,MSYprecLongT,FmsyprecLongT,pos=3,col="green")
abline(h=0)
legend(x='topright',lty=c("dotted","dashed","solid"),c("ShortT","MidT","LongT"),col="green",bty="n")


dev.off()

#############################







png(file=paste(outPath2,"plotsBaseCase/densityMSY_and_risks",prlength2,"_",nYrMean,".png",sep=""),width=8,height=8,units="in",res=200)






MSY<-data.frame(Bmsy=NA,Fmsy=NA,Msy=NA,iter=NA)
for (its in 1:nits)
{
sub<-res[res$iter==its,]
msy<-rev(order(sub$Yields))[1]
MSY<-rbind(MSY,data.frame(Bmsy=sub$SSBsP[msy],Fmsy=sub$perceivedF[msy],Msy=sub$Yields[msy],iter=its))
}
MSY<-MSY[-1,]

MSY<-MSY[order(MSY$Fmsy),]
MSY$count<-c(1:dim(MSY)[1])/max(MSY$iter)

medFmsy<-median(MSY$Fmsy)



plot(0,0,type="l", xlim=range(res$perceivedF) , ylim= c(0,1),xlab="perceived Fbar",ylab="",bty="n")
lines(MSY$Fmsy,MSY$count,col="red")
abline(v=medFmsy,col="red")
text(medFmsy,0.5,paste("med Fmsy",round(medFmsy,2)),pos=2,col="red")
lines(risk$perceivedF,risk$Risk1MidT/100,lty="dashed",col="green",lwd=2)
lines(risk$perceivedF,risk$Risk2MidT/100,lty="dotted",col="green",lwd=2)
lines(risk$perceivedF,risk$Risk3MidT/100,lty="solid",col="green",lwd=2)
abline(h=0.05,lty="dashed",col="green")
text(0.07,0.05,"5% risk",pos=3,col="green")
abline(h=0)
abline(v=FmsyprecMidT1,col="green",lty="dashed")
abline(v=FmsyprecMidT2,col="green",lty="dotted")
abline(v=FmsyprecMidT3,col="green",lty="solid")

legend(x="topleft",bty="n",col=c("red",rep("green",3)),lty=c("solid","dashed","dotted","solid"),c("Fmsy distrib","MidT risk1","MidT risk2","MidT risk3"))


dev.off()











png(file=paste(outPath2,"plotsBaseCase/otherRefPoints",prlength2,"_",nYrMean,".png",sep=""),width=8,height=8,units="in",res=200)

# plot for the SSBs
q95<-aggregate(SSBsT~perceivedF,data=res,function(x) {quantile(x,0.95)})
q75<-aggregate(SSBsT~perceivedF,data=res,function(x) {quantile(x,0.75)})
q50<-aggregate(SSBsT~perceivedF,data=res,function(x) {quantile(x,0.50)})
q25<-aggregate(SSBsT~perceivedF,data=res,function(x) {quantile(x,0.25)})
q05<-aggregate(SSBsT~perceivedF,data=res,function(x) {quantile(x,0.05)})

plot(q50$perceivedF,q50$SSBsT,type="l",col="red", lwd=2, xlab="perceived Fbar",ylab="SSBs",main="Spawning Stock Biomass",ylim=range(c(q05,q95)))
polygon(c(q25$perceivedF,rev(q75$perceivedF)),c(q25$SSBs,rev(q75$SSBsT)),col=rgb(1,0,0,0.6),border=F)
polygon(c(q05$perceivedF,rev(q95$perceivedF)),c(q05$SSBs,rev(q95$SSBsT)),col=rgb(1,0,0,0.3),border=F)
lines(q50$perceivedF,q50$SSBsT,lwd=2)

q50$diffBpa<-abs(q50$SSBsT-Bpa)
q50$diffBlim<-abs(q50$SSBsT-Blim)

Fpa  <-  q50$perceivedF[order(q50$diffBpa)[1]]
Flim <-  q50$perceivedF[order(q50$diffBlim)[1]]

segments(0,Bpa,Fpa,Bpa)
segments(0,Blim,Flim,Blim)
segments(Fpa,Bpa,Fpa,0)
segments(Flim,Blim,Flim,0)
text(Fpa,1.3*Bpa,paste("Fpa=",Fpa),pos=4)
text(Flim,Blim,paste("Flim=",Flim),pos=4)


Fmsy<-0.21
abline(v=Fmsy)
MSYBtrigger<-round(q05$SSBsT[q05$perceivedF==Fmsy],0)
segments(0,MSYBtrigger,Fmsy,MSYBtrigger)
segments(Fmsy,MSYBtrigger,Fmsy,0)
text(Fmsy,1.8*MSYBtrigger,paste("MSYBtrigger=",MSYBtrigger),pos=3)

q50<-aggregate(SSBsP~perceivedF,data=res[res$perceivedF<=0.6,],function(x) {quantile(x,0.50)})
lines(q50$perceivedF,q50$SSBsP,lwd=2,lty="dotted",bty="n")
legend(x="topright",c("true","perceived"),lty=c("solid","dotted"),bty="n")


dev.off()
