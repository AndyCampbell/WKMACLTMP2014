#-------------------------------------------------------------------------------
# WKHELP
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
# Performs an MSE of North Sea Herring under different TAC scenario's
#
# Date: 02-Sep-2012
#
# Build for R2.13.2
#-------------------------------------------------------------------------------

rm(list=ls())

library(FLSAM)
library(MASS)
library(msm)
wine <- F

library(FLCore)
#library(PBSadmb)
library(lattice)
library(MASS)

ac<-function(x) {return(as.character(x))}
an<-function(x) {return(as.numeric(x))}


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




source(paste(codePath,'functions.r',sep=""))

# load the true and observed stocks at the start of the simulation  from :
RecRegime <-  "srest"

#define year ranges
ShortT    <-  ac(2014:2018)
MidT      <-  ac(2019:2028)
LongT     <-  ac(2029:2052)



##-------------------------------------------------------------------------------
## Setup array to save results
##-------------------------------------------------------------------------------



diags<-data.frame(scenario=NA,RecRegime=NA,Iterations=NA,Btrigger=NA,Blim=NA,Ftarget=NA,
      Risk1ShortT=NA,Risk1MidT=NA,Risk1LongT=NA,
      Risk2ShortT=NA,Risk2MidT=NA,Risk2LongT=NA,
      Risk3ShortT=NA,Risk3MidT=NA,Risk3LongT=NA,
      SSBend=NA,
      meanSSBShortT=NA,meanSSBMidT=NA,meanSSBLongT=NA,
      Fend=NA,
      meanFShortT=NA,meanFMidT=NA,meanFLongT=NA,
      meanYieldShortTerm=NA,meanYieldMidTerm=NA,meanYieldLongTerm=NA,
      meanrelTACIAV=NA,noIAVrestrictup=NA,noIAVrestrictdown=NA,TACup=NA,TACdown=NA,
      SmeanAgeShortT=NA,SmeanAgeMidT=NA,SmeanAgeLongT=NA,
      SmeanWeightShortT=NA,SmeanWeightMidT=NA,SmeanWeightLongT=NA) 


##-------------------------------------------------------------------------------
## Load results
##-------------------------------------------------------------------------------


counter<-1 


for (sc in c("2.2mt","2.4mt"))
for (sc in c("2.4mt"))
{
for (TACvarlim in c(T,F))
{
for(Fvarlim in c(T,F))
{ 
for (opt in c(6,7,8,11,16,21))
{
scen          <- paste("LTMP",sc,sep="")         


BBscen        <- "noBB"   
LastRecOp     <- "geom"            
#
cat(scen,opt,"\n")
##
#
perm<-T
if (perm) cat("!!! scenario with permanent changes")
ifelse (perm,outPathp <- paste(outPath,"perm",sep=""), outPathp <- outPath)



sc.name<-paste(scen,opt,"_TACvarlim",TACvarlim,"_Fvarlim",Fvarlim,"_",BBscen,"_LastRec",LastRecOp,sep="")

outPath2<-paste(outPath,"HCR sensitivty varconst/",sc.name,"/",sep="")
source(paste(codePath,"07_scenarioDescription.r", sep=""))
mpPoints      <- get(scen)[[which(names(get(scen))==paste("opt",opt,sep=""))]]







load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalf.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalbiol.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalstocks.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalpercievedstocks.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_FinalTAC.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_FinalfSTF.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_FinalSSB.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalfishery.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_FinalmpPoints.RData",           sep=""))
source(paste(codePath,"functions.r",            sep=""))
source(paste(codePath,"04_forecastScenarios.r", sep=""))
load(file=paste(outPath,"settings.RData",       sep=""))
for(i in 1:length(settings)) assign(x=names(settings)[i],value=settings[[i]])




futureMaxYr<-histMinYr+length(biol@n[1,which(!is.na(biol@n[1,,,,,1])),,,,1])
projPeriod           <- 2014:(futureMaxYr-1)  
nits<-dim(f)[6]

print(counter)
Ref<-mpPoints
diags[counter,"scenario"]   <- opt
diags[counter,"RecRegime"]  <- RecType
diags[counter,"Iterations"] <- nits
diags[counter,"Btrigger"]   <- Ref$Btrigger
diags[counter,"Blim"]       <- Ref$Blim
diags[counter,"Ftarget"]    <- Ref$Ftarget
diags[counter,"permanent"]  <- perm
diags[counter,"TACvarlim"]  <- TACvarlim
diags[counter,"Fvarlim"]  <- Fvarlim


##-------------------------------------------------------------------------------
## Diagnostics on results
##-------------------------------------------------------------------------------

Btrg  <- Ref$Btrigger
Blim  <- Ref$Blim
Bpa  <- Ref$Bpa


Ssb<-ssbb(biol,f,stockstore)
percSsb<-ssb(stockstore)
Fbar<-quantMeans((f[ac(4:8),]))
Fbar2<-quantMeans((harvest(stockstore)[ac(4:8),]))

# risk related to Blim      .

risk<-  apply(Ssb<mpPoints$Blim,c(1:5),sum)/nits

diags[counter,"Risk1ShortT"]             <- mean(risk[,ShortT])*100       # percentage of iteration that reach Blim
diags[counter,"Risk1MidT"]               <- mean(risk[,MidT])  *100
diags[counter,"Risk1LongT"]              <- mean(risk[,LongT]) *100

diags[counter,"Risk2ShortT"]             <- length(unique(which(Ssb[,ShortT]<Blim,arr.ind=T)[,"dim6"]))/dims(Ssb)$iter  *100       # percentage of iteration that reach Blim
diags[counter,"Risk2MidT"]               <- length(unique(which(Ssb[,MidT]  <Blim,arr.ind=T)[,"dim6"]))/dims(Ssb)$iter  *100
diags[counter,"Risk2LongT"]              <- length(unique(which(Ssb[,LongT] <Blim,arr.ind=T)[,"dim6"]))/dims(Ssb)$iter  *100

diags[counter,"Risk3ShortT"]             <- max(risk[,ShortT])*100       # percentage of iteration that reach Blim
diags[counter,"Risk3MidT"]               <- max(risk[,MidT])  *100
diags[counter,"Risk3LongT"]              <- max(risk[,LongT]) *100


# risk related to Btrig      .

risk<-  apply(Ssb<mpPoints$Btrigger,c(1:5),sum)/nits

diags[counter,"RiskBT1ShortT"]             <- mean(risk[,ShortT])*100       # percentage of iteration that reach Blim
diags[counter,"RiskBT1MidT"]               <- mean(risk[,MidT])  *100
diags[counter,"RiskBT1LongT"]              <- mean(risk[,LongT]) *100

diags[counter,"RiskBT2ShortT"]             <- length(unique(which(Ssb[,ShortT]<Blim,arr.ind=T)[,"dim6"]))/dims(Ssb)$iter  *100       # percentage of iteration that reach Blim
diags[counter,"RiskBT2MidT"]               <- length(unique(which(Ssb[,MidT]  <Blim,arr.ind=T)[,"dim6"]))/dims(Ssb)$iter  *100
diags[counter,"RiskBT2LongT"]              <- length(unique(which(Ssb[,LongT] <Blim,arr.ind=T)[,"dim6"]))/dims(Ssb)$iter  *100

diags[counter,"RiskBT3ShortT"]             <- max(risk[,ShortT])*100       # percentage of iteration that reach Blim
diags[counter,"RiskBT3MidT"]               <- max(risk[,MidT])  *100
diags[counter,"RiskBT3LongT"]              <- max(risk[,LongT]) *100



# number of times SSB<Btrigger
trig<-  apply(Ssb<mpPoints$Btrigger,c(1:5),sum)
trig<-  Ssb<mpPoints$Btrigger

diags[counter,"bellowBtrigShortT"]      <- median(c(yearSums(trig[,ShortT])/length(ShortT)@.Data))
diags[counter,"bellowBtrigMidT"]       <- median(c(yearSums(trig[,MidT])/length(MidT)@.Data))
diags[counter,"bellowBtrigLongT"]     <- median(c(yearSums(trig[,LongT])/length(LongT)@.Data))



# stock and fishing mortality
diags[counter,"SSBend"]            <- round(median(c(apply(Ssb[,ac(futureMaxYr-1)],3:6,mean,na.rm=T))))   # or round(iterMeans(Ssb[,ac(futureMaxYr-1)]))  ?
diags[counter,"meanSSBLongT"]      <- round(median(c(apply(Ssb[,LongT],3:6,mean,na.rm=T))))
diags[counter,"meanSSBMidT"]       <- round(median(c(apply(Ssb[,MidT],3:6,mean,na.rm=T))))
diags[counter,"meanSSBShortT"]     <- round(median(c(apply(Ssb[,ShortT],3:6,mean,na.rm=T))))




diags[counter,"meanFShortT"]        <- round( median(   yearMeans(Fbar[,ShortT])@.Data) ,3)
diags[counter,"meanFMidT"]        <- round( median(   yearMeans(Fbar[,MidT])@.Data) ,3)
diags[counter,"meanFLongT"]        <- round( median(   yearMeans(Fbar[,LongT])@.Data) ,3)
diags[counter,"Fend"]                 <- round( median   (Fbar[,ac(futureMaxYr-1)]@.Data) ,3)

# difference between percieved and true stocks
diags[counter,"SSBabsBias"]           <-  round(apply( yearMeans(100*abs(percSsb[,ac(projPeriod)]-Ssb[,ac(projPeriod)])/Ssb[,ac(projPeriod)]),1:5,median,na.rm=T)@.Data,3)
diags[counter,"SSBBias"]              <-  round(apply( yearMeans(100*(percSsb[,ac(projPeriod)]-Ssb[,ac(projPeriod)])/Ssb[,ac(projPeriod)]),1:5,median,na.rm=T)@.Data,3)
diags[counter,"FbarabsBias"]          <-  round(apply( yearMeans(100*abs((Fbar2[,ac(projPeriod)])-(Fbar[,ac(projPeriod)]))/(Fbar[,ac(projPeriod)])),1:5,median,na.rm=T)@.Data,3)
diags[counter,"FbarBias"]             <-  round(apply( yearMeans(100*((Fbar2[,ac(projPeriod)])-(Fbar[,ac(projPeriod)]))/(Fbar[,ac(projPeriod)])),1:5,median,na.rm=T)@.Data,3)



# catches and quotas
diags[counter,"meanYieldShortTerm"]   <- round(median(c(yearMeans((computeLandings(fishery)[,ShortT])))))
diags[counter,"meanYieldMidTerm"]     <- round(median(c(yearMeans((computeLandings(fishery)[,MidT])))))
diags[counter,"meanYieldLongTerm"]    <- round(median(c(yearMeans((computeLandings(fishery)[,LongT])))))

# mean age and weight in the catches and in the SSB

      # mean age mature fish
sma<-quantSums(sweep((biol@n*biol@fec),c(1:6) ,c(0:12),"*"))/quantSums((biol@n*biol@fec))
      # mean weight mature fish
smw<-quantSums((biol@n*biol@wt*biol@fec)) / quantSums((biol@n*biol@fec))

diags[counter,"SmeanAgeShortT"]   <- round(median(c(yearMeans((sma[,ShortT])))),2)
diags[counter,"SmeanAgeMidT"]     <- round(median(c(yearMeans((sma[,MidT])))),2)
diags[counter,"SmeanAgeLongT"]    <- round(median(c(yearMeans((sma[,LongT])))),2)

diags[counter,"SmeanWeightShortT"]   <- round(1000*median(c(yearMeans((smw[,ShortT])))),0)
diags[counter,"SmeanWeightMidT"]     <- round(1000*median(c(yearMeans((smw[,MidT])))),0)
diags[counter,"SmeanWeightLongT"]    <- round(1000*median(c(yearMeans((smw[,LongT])))),0)


      # mean age based on the distribution of catches between age classes
yma<-quantSums(sweep((fishery@landings.n),c(1:6) ,c(0:12),"*"))/quantSums((fishery@landings.n))
      # mean weight based on the distribution of catches between age classes
ymw<-quantSums((fishery@landings.n*fishery@landings.wt))/quantSums((fishery@landings.n))

diags[counter,"YmeanAgeShortT"]   <- round(median(c(yearMeans((yma[,ShortT])))),2)
diags[counter,"YmeanAgeMidT"]     <- round(median(c(yearMeans((yma[,MidT])))),2)
diags[counter,"YmeanAgeLongT"]    <- round(median(c(yearMeans((yma[,LongT])))),2)

diags[counter,"YmeanWeightShortT"]   <- round(1000*median(c(yearMeans((ymw[,ShortT])))),0)
diags[counter,"YmeanWeightMidT"]     <- round(1000*median(c(yearMeans((ymw[,MidT])))),0)
diags[counter,"YmeanWeightLongT"]    <- round(1000*median(c(yearMeans((ymw[,LongT])))),0)







# quota variability
diags[counter,"meanrelTACIAV"]        <- round(median(c(apply(abs(TAC[,ac(projPeriod[2]:rev(projPeriod)[1])] - TAC[,ac(projPeriod[1]:rev(projPeriod)[2])]) /TAC[,ac(projPeriod[1]:rev(projPeriod)[2])] * 100,3:6,mean,na.rm=T))),3)
IAVUp   <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)])] == 1.2* TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],arr.ind=T)
IAVDown <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)])] == 0.8* TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],arr.ind=T)
#  #- Average number of times the IAV rule is applied upwards or downwards
diags[counter,"noIAVrestrictup"]<-0
if((nrow(IAVUp)) > 0 ){
  a <- IAVUp
  diags[counter,"noIAVrestrictup"]    <- max(0,median(aggregate(a[,"dim2"],by=list(a[,"dim6"]),function(x){length(x)})$x),na.rm=T)
}
 diags[counter,"noIAVrestrictdown"]<-0
if((nrow(IAVDown)) > 0 ){
  a <- IAVDown
  diags[counter,"noIAVrestrictdown"]  <- max(0,median(aggregate(a[,"dim2"],by=list(a[,"dim6"]),function(x){length(x)})$x),na.rm=T)
}
#
#  #- Which TAC of the runs go up and which go down
resUp   <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)])] > TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],arr.ind=T)
resDown <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)])] < TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],arr.ind=T)
#
#  #- Mean increase in TAC is TAC goes up, or mean decrease in TAC is TAC goes down
diags[counter,"TACup"]   <- round(mean((TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1])] - TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2])])@.Data[which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1])] > TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2])])]))
diags[counter,"TACdown"] <- round(mean((TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2])] - TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1])])@.Data[which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1])] < TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2])])]))
#




counter <- counter + 1
}
}
}
write.csv((diags),file=paste(outPath,"/tables_diags_HCR sensitivity varconst",BBscen,"_LastRec.csv",sep=""),row.names=T)

}

#write.csv(t(diags),file=paste(outPath,"/tables_diags",paste(scen,opt,"_TACvarlim",TACvarlim,"_Fvarlim",Fvarlim,"_",BBscen,"_LastRec",LastRecOp,sep=""),".csv",sep=""),row.names=T)

