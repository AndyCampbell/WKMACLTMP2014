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
library(msm)
library(FLFleet)
library(minpack.lm)

#for (opt_i in 1:21){
opt_i <- 1



wine <- F


path          <- getwd()      #"W:/IMARES/Data/ICES-WG/WKMACLTMP/"
inPath        <- "Data/"      #"W:/IMARES/Data/ICES-WG/WKMACLTMP/Data/"
codePath      <- "R code/"    #"W:/IMARES/Data/ICES-WG/WKMACLTMP/R code/"
outPath       <- "Results/"   #"W:/IMARES/Data/ICES-WG/WKMACLTMP/Results/"
                              # if(substr(R.Version()$os,1,3)== "lin"){
                              #   path        <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",path)
                              #   inPath      <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",inPath)
                              #   codePath    <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",codePath)
                              #   outPath     <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",outPath)
                              #   }
                              
                              # 
                              # home<-F
                              # if(home)
                              # {
                              # path          <- "D://MSE/"
                              # inPath        <- "D://MSE/Data/"
                              # codePath      <- "D://MSE/R code/"
                              # outPath       <- "D://MSE/Results/"
                              # }
                              # 


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
load(file=paste(outPath,"SRmod.RData",          sep=""))
load(file=paste(outPath,"resRFinal.RData",      sep=""))
for(i in 1:length(settings)) assign(x=names(settings)[i],value=settings[[i]])

source(paste(codePath,"functions_modified.r",            sep=""))
source(paste(codePath,"04_forecastScenarios.r", sep=""))


#- Settings
load(file=paste(outPath,"settings.RData",sep=""))
histMinYr   <- settings$histMinYr
histMaxYr   <- settings$histMaxYr
nyrs        <- settings$nyrs
futureMaxYr <- settings$futureMaxYr
histPeriod  <- settings$histPeriod
projPeriod  <- settings$projPeriod
recrPeriod  <- settings$recrPeriod
selPeriod   <- settings$selPeriod
fecYears    <- settings$fecYears
nits        <- settings$nits
RecType<-settings$RecType


# set manually the projection period
prlength    <-  40
projPeriod    <-ac((histMaxYr+1):(histMaxYr+prlength))
futureMaxYr   <- an(rev(projPeriod)[1]    )

settings$projPeriod<-projPeriod
settings$futureMaxYr<-futureMaxYr

nits <- 1000
short <- T
load(file=paste(outPath,"resNFinal_",ifelse(short,max(200,nits),nits),".RData",      sep=""))
load(file=paste(outPath,"resFFinal_",ifelse(short,max(200,nits),nits),".RData",      sep=""))


# run a shorter number of iterations : has to be 100

if (short) 
{
nits      <-  1000
biol      <-  biol[,,,,,1:nits]
fishery   <-  iter(fishery,1:nits)
ctch      <-  ctch[,,,,,1:nits]
Rindex    <-  Rindex[,,,,,1:nits]
stocks    <-  stocks[,,,,,1:nits]
SRmod     <-  SRmod[1:nits,]
devR      <-  devR[,,,,,1:nits]

}

# 
# 
 if (short) 
 {
# devN<-devN100
# devF<-devF100
 devN<-devN[,,,,,1:nits]
 devF<-devF[,,,,,1:nits]
 #rm(devN100)
 #rm(devF100)
 }

#
  #------------------------------------------------------------------------------#
  # 0) setup TACS & F's and Survivors and maximum change in effort
  #------------------------------------------------------------------------------#
maxEff                                <- 1000
TAC                                   <- FLQuant(NA,dimnames=list(age="all",year=histMinYr:(futureMaxYr+3),unit=c("A"),season="all",area=1,iter=1:nits))
TAC[,ac(2000:2014),"A"]               <- c(612,670,683,583,532,422,444,502,458,605,885,959,927,906,1396) *1000
#TACusage                              <- FLQuant(array(rep(c(rep(1,length(histMinYr:futureMaxYr)+3),rep(0.539755,length(histMinYr:futureMaxYr)+3),
#                                                   rep(1,length(histMinYr:futureMaxYr)+3),rep(1,length(histMinYr:futureMaxYr)+3)),nits),dim=c(1,length(histMinYr:futureMaxYr)+3,4,1,1,nits)),dimnames=dimnames(TAC[,ac(histMinYr:(futureMaxYr+3))]))
TACusage                              <- FLQuant(array(rep(c(rep(1,length(histMinYr:futureMaxYr)+3)),nits),dim=c(1,length(histMinYr:futureMaxYr)+3,1,1,1,nits)),dimnames=dimnames(TAC[,ac(histMinYr:(futureMaxYr+3))]))

HCRTAC                                <- TAC; HCRTAC[] <- NA; SSB <- HCRTAC[,,1]; HCRSSB <- SSB

stockstore <- stocks                    # object to store the perceived stocks. The object "stocks" being updated at each time step, it does keep track of the percieved stock

f                                     <- FLQuant(NA,dimnames=dimnames(fishery@landings.n))
for(iFsh in dimnames(f)$unit)
  f[,ac(histMinYr:histMaxYr),iFsh]    <- sweep(harvest(stocks)[,ac(histMinYr:histMaxYr)],c(1,3:5),propN[,,iFsh],"*")
fSTF                                  <- f; fSTF@.Data[] <- NA
survivors                             <- FLQuant(NA,dimnames=dimnames(n(biol)[,ac(2011)]))

  #------------------------------------------------------------------------------#
  # 1) Select Scenario's
  #------------------------------------------------------------------------------#

scen          <- c("LTMP")              # 
opt           <- opt_i                     # for multiple scenario combinations a counter
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
                                                     # "geom" = replace by geomean 1990/(TaY-1)
                                                     # "RCT3" = replace by RCT3 output


mpOptions<-list(opt=opt,TACvarlim=TACvarlim,Fvarlim=Fvarlim,BBscen=BBscen)

HCR.method='without.stf'   #stf, without.stf, fix.prop


sc.name<-paste(scen,opt,"_TACvarlim",TACvarlim,"_Fvarlim",Fvarlim,"_",BBscen,"_LastRec",LastRecOp,"_wstf",sep="")
outPath2<-paste(outPathp,sc.name,"/",sep="")
source(paste(codePath,"07_scenarioDescription_wstf.r", sep=""))
mpPoints      <- get(scen)[[which(names(get(scen))==paste("opt",opt,sep=""))]]

# mpPoints$Btrigger <- 2*10^6
# mpPoints$Ftarget <- 0.15
#mpPoints[['alpha']] <- 0.2
#mpPoints[['beta']] <- 0

dir.create(outPath2)

  #------------------------------------------------------------------------------#
  # 2) Start running
  #------------------------------------------------------------------------------#


start.time <- Sys.time()
for (iYr in an(projPeriod)){
  cat(iYr,"\n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

  #- Define mortality rates for iYr-1 to calculate survivors to iYr
  m           <- m(biol)[,ac(iYr-1),,]
  z           <- unitSums(f[,ac(iYr-1),,,,]) + m

   # - previous year recruitment 
if ((iYr-1)>histMaxYr)
{ 
     ssbtp<-ssbb(biol[,ac(iYr-1),,,,],f[,ac(iYr-1),,,,],stockstore[,ac(iYr-1),,,,])
     for (its in 1:nits)   n(biol)[1,ac(iYr-1),,,,its][] <- B2Rbis (ssbtp[,,,,,its],SRmod[its,],its,iYr-1,devR)
  }     
      

  #- Update biological model to iYr
    #- Survivors
  survivors   <- n(biol)[,ac(iYr-1)] * exp(-z)
  n(biol)[ac((range(biol,"min")+1):range(biol,"max")),ac(iYr),,] <- survivors[-dim(survivors)[1],,,,,]@.Data

    
     
     
           
    #- Plusgroup
  if (!is.na(range(biol,"plusgroup"))){
    n(biol)[ac(range(biol,"max")),ac(iYr),] <- n(biol)[ac(range(biol,"max")),ac(iYr),] + survivors[ac(range(biol,"max"))]
  }

  cat("\n Finished biology \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
 
 
 
 
 
  #- Update fishery to year iYr-1
  landings.n(fishery)[,ac(iYr-1)]     <- sweep(sweep(f[,ac(iYr-1),,,,],c(1:2,4:6),z,"/"),c(1:2,4:6),n(biol)[,ac(iYr-1)]*(1-exp(-z)),"*")

  #- Create stock object for assessment
  yrmin1      <- iYr -1
  TaY         <- yrmin1               #Terminal assessment year
  ImY         <- TaY+1                #Intermediate Year
  FcY         <- TaY+2                #Forecast year

 idxyrmin1   <- which(dimnames(biol@n)$year == yrmin1)
#  tmp_biol    <- biol[,1:idxyrmin1]   #Same but faster as window(biol,histMinYr,yrmin1)
#  tmp_fishery <- window(fishery,histMinYr,yrmin1)
#  tmp_stocks  <- stocks[,1:idxyrmin1] #Same but faster as window(stocks,histMinYr,yrmin1)
#
  #- Update stocks to year iYr -1
#  tmp_stocks  <- updateStocks(tmp_stocks,tmp_fishery,yrmin1,tmp_biol,ctch)


                 stocks@catch.n[,ac(yrmin1)]        <- unitSums(catch.n(fishery)[,ac(yrmin1)])   * ctch[,ac(yrmin1)]
                 stocks@landings.n[,ac(yrmin1)]     <- unitSums(landings.n(fishery)[,ac(yrmin1)])* ctch[,ac(yrmin1)]
                 stocks@landings.wt[,ac(yrmin1)]    <- stocks@catch.wt[,ac(yrmin1)]
                 stocks@discards.n[,ac(yrmin1)]     <- 0
                 stocks@discards.wt[,ac(yrmin1)]    <- 0
                 stocks@catch[,ac(yrmin1)]          <- computeCatch(stocks)[,ac(yrmin1)]
                 stocks@landings[,ac(yrmin1)]       <- computeLandings(stocks)[,ac(yrmin1)]
                 stocks@discards[,ac(yrmin1)]       <- computeDiscards(stocks)[,ac(yrmin1)]

                 



  #- Overwrite results from update to stock again (but not for 2011, as that result is already known)
#  if(iYr > an(projPeriod[1]))
#    stocks      <- tmp2stocks(stocks,tmp_stocks,TaY)
# #   
#                TaYtmp_stocks         <- tmp_stocks[,ac(TaY)]
#                TaYstocks             <- stocks[,ac(TaY)]
#                TaYtmp_stocks@stock   <- TaYstocks@stock
#                TaYtmp_stocks@stock.n <- TaYstocks@stock.n
#                TaYtmp_stocks@harvest <- TaYstocks@harvest
#                TaYtmp_stocks@name    <- TaYstocks@name
#                TaYtmp_stocks@desc    <- TaYstocks@desc
#                TaYtmp_stocks@range   <- TaYstocks@range
#                stocks[,ac(TaY)]       <- TaYtmp_stocks[,ac(TaY)] 
#            
#   
    
    
    

  cat("\n Finished update stocks\n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))


   Rindex@index[,ac(TaY)]  <-   exp ( log(n(biol)[1,ac(TaY)] *  Rindex@index.q[,ac(TaY)]) + rnorm(nits,0,c(Rindex@index.var[,ac(TaY)]@.Data)))
  cat("\n Finished update survey\n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))




  #-Do the assessment
  
    # length of the assessment period (should be 2013-1980)
  lassess<-dim(devN)[2]
  assessPeriod<-(TaY-lassess+1):TaY
  # for year ouside of period for which assessment errors are calculated (ie too old or future years) assume no error :
  stocks@stock.n<-biol@n
  stocks@harvest <- unitSums(f)
  
  # for the relevant years, apply the assessment errors
  stocks@stock.n[,ac(assessPeriod)] <- biol@n[,ac(assessPeriod)]        * exp(devN[,,ac(iYr),,,])
  stocks@harvest[,ac(assessPeriod)] <- unitSums((f[,ac(assessPeriod)])) * exp(devF[,,ac(iYr),,,])
  stocks@stock   <- computeStock(stocks)
  
  # overwrite the last estimated recruitment?
  if (LastRecOp == "geom")    stocks@stock.n[1,ac(TaY)] <-       exp(yearMeans(log(stock.n(stocks)[1,ac(1990:(TaY-1))])))
  
  
  if (LastRecOp == "RCT3")
  {
  iTer<-nits
      for (i in 1:iTer)
      {
      # prepare init file for RCT3
      R<-iter(stock.n(stocks)[ac(0),ac(1990:TaY)],iTer)
      R<-c(R@.Data)
      R[length(R)]<--11
      IBTS.index<-c(rep(-11,8),c(iter(Rindex@index[,ac(1998:TaY)],iTer)@.Data))
      years<-1990:TaY
      
      
      # remove files in the RCT3 folder !!!!
      file.name<-paste(outPath2,"RCT3init.txt",sep="")
      file.remove(file=file.name)
      write.table(data.frame("RCT3 for NEA Mackerel"),file=file.name,quote=F,col.names=FALSE,row.names=FALSE,append=TRUE,sep="\t")
      write.table(data.frame(1,length(R),2,"SAM","IBTS.index"),file=file.name,quote=F,col.names=FALSE,row.names=FALSE,append=TRUE,sep="\t")
      write.table(data.frame(years,R,IBTS.index),file=file.name,col.names=FALSE,quote=F,row.names=FALSE,append=TRUE,sep="\t")
      write.table(data.frame(c("SAM","IBTS.index")),file=file.name,col.names=FALSE,quote=F,row.names=FALSE,append=TRUE,sep="\t")
      
      source(paste(codePath,"RCT3v4a.r",sep=""))
      Rct3<-RCT3(file.name,logged=T)
      RCT3res<-Rct3$output()
     
      
      stocks@stock.n[1,ac(TaY),,,,iTer]   <-      RCT3res$Years$WAPred  
      }
  }
  
  
  
  # copy the perception of the stock in terminal assessment year to the stockstore object
  stockstore[,ac(TaY)]@stock.n  <-    stocks[,ac(TaY)]@stock.n
  stockstore[,ac(TaY)]@harvest  <-    stocks[,ac(TaY)]@harvest
  stockstore[,ac(TaY)]@catch.wt  <-   stocks[,ac(TaY)]@catch.wt
  
  
  # survivors for the short term forecast
   dimnames(survivors)$year<-ac(iYr)
   
  # recruitment for the first year in the short term forecast is the geometric mean of the historical time series 
  survivors[ac(0),]                  <- exp(yearMeans(log(stock.n(stocks)[1,ac(1990:(TaY-1))])))

  #Set plusgroup at 12 (which is true plusgroup - recruitment)
  survivors[-1,]    <- FLQuant(setPlusGroup(stocks[,ac(TaY)]@stock.n * exp(-stocks[,ac(TaY)]@harvest-stocks[,ac(TaY)]@m),11)@.Data,
                               dimnames=list(age=dimnames(stocks@stock.n)$age[-1],year=ac(TaY),unit=dimnames(stocks@stock.n)$unit,
                                             season=dimnames(stocks@stock.n)$season,area=dimnames(stocks@stock.n)$area,iter=1:nits))

  cat("\n Finished stock assessment \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
 
 
 
  #- Project 

  res <- HCR.TAC(HCR.method,stocks[,1:idxyrmin1],survivors,window(fishery,histMinYr,yrmin1),iYr,TAC,mpPoints$scen,histMaxYr,mpPoints,mpOptions)

      
  
#  res   <- projectMac(stocks[,1:idxyrmin1],survivors,window(fishery,histMinYr,yrmin1),iYr,TAC,mpPoints$scen,NULL,histMaxYr,mpPoints,mpOptions)
  
  TAC[,ac(FcY)]             <- res[["TAC"]]
  HCRTAC[,ac(FcY)]          <- res[["HCRTAC"]]
  HCRSSB[,ac(FcY)]          <- res[["SSB"]][["HCRSSB"]][,ac(FcY)]
  SSB[,ac(FcY)] <- res[["SSB"]][["SSB"]][,ac(FcY)]
  
  if(iYr != rev(projPeriod)[1]) fSTF[,ac(FcY)]            <- res[["fSTF"]]

  cat("\n Finished forecast \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  #-Calculate effort accordingly (assuming constant catchability)
  f[,ac(ImY)]               <- sweep(catch.sel(fishery)[,ac(ImY)],c(3,6),pmin(maxEff,f31tF(TAC*TACusage,biol,as.numeric(ImY),fishery)),"*")
  cat("\n Finished effort calc \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

  #- Save each year the output
  #save.image(file=paste(outPath,scen,"_",opt,"_",mpPoints$FadultA,"_",iYr,".RData",sep=""))
  #save.image(file=paste("/home/hintz001/WKHELP_test2_",iYr,".RData",sep=""))
  #save(file=paste("D:/WKHELP_test3_",iYr,".RData",sep=""))
}

stockstore@landings.n   <- stockstore@harvest * stockstore@stock.n * (1- exp(-stockstore@harvest - stockstore@m)) / (stockstore@harvest + stockstore@m)
stockstore@landings.wt<-stockstore@catch.wt
  #-------------------------------------------------------------------------------
  # 3): Save the objects
  #-------------------------------------------------------------------------------

save(biol          ,file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalbiol.RData",        sep=""))
save(fishery       ,file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalfishery.RData",     sep=""))
save(Rindex        ,file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalsurveys.RData",     sep=""))
save(stocks        ,file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalstocks.RData",      sep=""))
save(stockstore    ,file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalpercievedstocks.RData",      sep=""))
save(TAC           ,file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_FinalTAC.RData",         sep=""))
save(HCRTAC        ,file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_FinalHCRTAC.RData",      sep=""))
save(f             ,file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalf.RData",           sep=""))
save(fSTF          ,file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_FinalfSTF.RData",        sep=""))
save(mpPoints      ,file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_FinalmpPoints.RData",    sep=""))
save(settings      ,file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalsettings.RData",    sep=""))


}
