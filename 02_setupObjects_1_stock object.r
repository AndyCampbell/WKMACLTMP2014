#-------------------------------------------------------------------------------
# Mackerel management plan evaluation
#
# Author: Thomas Brunel    (based on Niels Hintzen's code for North Sea herring management plan evaluation)
#         IMARES, The Netherland
#
# Performs an MSE of NEA Mackerel under different TAC scenario's
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


#-------------------------------------------------------------------------------
#
#  WARNING    : the mcmc function generate samples of the stock uses a fixed
#   seed, so that if this code is run again, the output of the mcmc procedure is 
#   the same, and matches the SR parameters estimated outside this code
#
#-------------------------------------------------------------------------------




rm(list=ls())
library(FLCore)
#library(FLAssess)
library(FLSAM)
library(MASS)
library(msm)

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

setwd(path)

home<-F
if(home)
{
path          <- "D://MSE/"
inPath        <- "D://MSE/Data/"
codePath      <- "D://MSE/R code/"
outPath       <- "D://MSE/Results/"
}

#- Load objects
load(file=paste(outPath,"Mac.RData",     sep=""))
load(file=paste(outPath,"Macsam.RData",  sep=""))
load(file=paste(outPath,"Mactun.RData",  sep=""))
load(file=paste(outPath,"/Macctrl.RData", sep=""))

#- Settings
assess.name <-  "NEA-Mac-WGWIGE-2014-V2"
histMinYr   <- 1980
histMaxYr   <- 2013
nyrs        <- 100
futureMaxYr <- histMaxYr + nyrs
histPeriod  <- ac(histMinYr:histMaxYr)
projPeriod  <- ac((histMaxYr+1):futureMaxYr)
recrPeriod  <- ac(1990:(histMaxYr-1))
selPeriod   <- ac(1990:2013)
fecYears    <- ac(1990:2013)
nits        <- 1000
                    
                    
RecType     <-  "srest"                # chose from
BiolType    <-  "AUTOCORperm"                 
                                          # chose from    ARMArev      : ARMA process assuming reversibility of recent changes
                                          #               ARMAperm     : ARMA process assuming permantent recent changes
                                          #               AUTOCORperm  : simple autocorrelated process assuming permantent recent changes
                                          #               MEAN         :  a simple, mean,  constant value (useless)


settings    <- list(histMinYr=histMinYr,histMaxYr=histMaxYr,futureMaxYr=futureMaxYr,
                    histPeriod=histPeriod,projPeriod=projPeriod,recrPeriod=recrPeriod,
                    nyrs=nyrs,nits=nits,fecYears=fecYears,RecType=RecType,BiolType=BiolType)


source(paste(codePath,"functions.r",sep=""))


  #-------------------------------------------------------------------------------
  # 1): Create stock object & use vcov for new realisations
  #-------------------------------------------------------------------------------
Mac.sam@control<-Mac.ctrl
run.dir<-paste(inPath,"/",assess.name,"/run",sep="")

stocks                            <- monteCarloStock2(Mac,Mac.sam,nits,run.dir=run.dir,seed_number=floor(pi*10000))
stocks                            <- window(stocks,start=histMinYr,end=futureMaxYr)
stocks@catch.n                    <- stocks@stock.n * stocks@harvest / (stocks@harvest + stocks@m) * (1 - exp(-stocks@harvest - stocks@m))
stocks@landings.n                 <- stocks@catch.n

stocks@harvest.spwn[,projPeriod]  <- stocks@harvest.spwn[,ac(histMaxYr)]
stocks@m.spwn[,projPeriod]        <- stocks@m.spwn[,ac(histMaxYr)]



    
 #------------------------------------------------------------------------------
  # 2): prepare the SR pairs for the SR parameters estimation using ADMB by Einar
  #-------------------------------------------------------------------------------
library(mac)    
#### create a table with the SR pairs from the assessment

ssbest <-ssb(Mac)[,ac(recrPeriod)]
r<-rec(Mac)[,ac(recrPeriod)]
year<-recrPeriod
sr_data<-data.frame(year,r=c(r@.Data),ssb=c(ssbest@.Data))



#### combined SR pairs of all replicates in a single table    
x <- melt(ssb(stocks))[,c("iter","year","value")]
names(x)[3] <- "ssb"
y <- melt(rec(stocks))[,c("iter","year","value")]
names(y)[3] <- "rec"
sr_data_mcmc <- join(y,x)
sr_data_mcmc <- sr_data_mcmc[sr_data_mcmc$year %in% 1990:2012,]
    

  #-------------------------------------------------------------------------------
  # 3): Save the objects
  #-------------------------------------------------------------------------------

saveobjects <- F
if (saveobjects)
{
  save(sr_data       ,file=paste(outPath,"sr_data.RData",       sep=""))    
  save(sr_data_mcmc  ,file=paste(outPath,"sr_data_mcmc.RData",  sep=""))    
  save(stocks        ,file=paste(outPath,"stocks_save.RData",   sep=""))
  save(settings      ,file=paste(outPath,"settings.RData",      sep=""))
}


  #-------------------------------------------------------------------------------
  # 4): Continue setup objects in file 02_setupObjects_2_other objects
  #-------------------------------------------------------------------------------
