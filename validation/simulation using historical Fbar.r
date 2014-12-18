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
require(minpack.lm)

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


outPathp  <-  outPath
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

source(paste(codePath,"functions.r",            sep=""))
source(paste(codePath,"04_forecastScenarios.r", sep=""))
source(paste(codePath,"validation/FourierSurrogates.r", sep=""))

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
nits      <-  1000

# set manually the projection period
prlength    <-  length(histPeriod)                              # max 100
projPeriod    <-ac((histMaxYr+1):(histMaxYr+prlength))
futureMaxYr   <- an(rev(projPeriod)[1]    )

settings$projPeriod<-projPeriod
settings$futureMaxYr<-futureMaxYr



# run a shorter number of iterations : has to be 100
short<-T
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




load(file=paste(outPath,"resNFinal_simple.RData",      sep=""))
load(file=paste(outPath,"resFFinal_simple.RData",      sep=""))



if (short) 
{
devN<-devN[,,,,,1:nits]
devF<-devF[,,,,,1:nits]
}

  
 # generate fourrier surrowgates of the historical Fbar time series 
 #get Fbar 
fbar<-quantMeans(harvest(Mac)[ac(4:8),histPeriod])   
fbarsim<-FourierSurrogates(c(fbar@.Data),100)  
fbar<-expand(fbar,year=1980:futureMaxYr)  

# finding one surrogate which starts close to Fbar2013
welke<-order(abs(fbarsim[1,]-c(fbar[,ac(histMaxYr)]@.Data)))[11]
  
fbar[,projPeriod]<-fbarsim[,welke]  
plot(fbar,col=c(rep("blue",length(histPeriod)),
              rep("red",length(projPeriod))),
              main="historic (blue) and surrogate (red) Fbar",
              ylab= "Fbar")  
  
  
  
  
  
  
  
  
               


        
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
              
              stockstore <- stocks # object to store the perceived stocks. The object "stocks" being updated at each time step, it does keep track of the percieved stock
              
              f                                     <- FLQuant(NA,dimnames=dimnames(fishery@landings.n))
              for(iFsh in dimnames(f)$unit)
                f[,ac(histMinYr:histMaxYr),iFsh]    <- sweep(harvest(stocks)[,ac(histMinYr:histMaxYr)],c(1,3:5),propN[,,iFsh],"*")
              fSTF                                  <- f; fSTF@.Data[] <- NA
              survivors                             <- FLQuant(NA,dimnames=dimnames(n(biol)[,ac(2011)]))

                #------------------------------------------------------------------------------#
                # 1) Select Scenario's   not relevant  here
                #------------------------------------------------------------------------------#
              
              outPath2       <- paste(outPath,"/validation/",sep="")
              
              scen          <- "surrogateFbar"              # 
              opt           <- ""                      # for multiple scenario combinations a counter
              if(perm==F) opt <- ""
              mpOptions<-list()

              mpPoints      <- list(Fequ="FS")
              #
            
              
                #------------------------------------------------------------------------------#
                # 2) Start running
                #------------------------------------------------------------------------------#
              

              for (iYr in an(projPeriod)){

                cat(iYr,"\n")
                #- Define mortality rates for iYr-1 to calculate survivors to iYr
                m           <- m(biol)[,ac(iYr-1),,]
                z           <- (f[,ac(iYr-1),,,,]) + m
              
                              
               # - previous year recruitment 
                if ((iYr-1)>histMaxYr)
                { 
                     ssbtp<-ssbb(biol[,ac(iYr-1),,,,],f[,ac(iYr-1),,,,],stockstore[,ac(iYr-1),,,,])
                     n(biol)[1,ac(iYr-1)]<-B2RF(ssbtp,SRmod,iYr-1,devR)
                }    
              
                #- Update biological model to iYr
                  #- Survivors
                survivors   <- n(biol)[,ac(iYr-1)] * exp(-z)
                n(biol)[ac((range(biol,"min")+1):range(biol,"max")),ac(iYr),,] <- survivors[-dim(survivors)[1],,,,,]@.Data
                
              
                #- Plusgroup
                if (!is.na(range(biol,"plusgroup"))){
                  n(biol)[ac(range(biol,"max")),ac(iYr),] <- n(biol)[ac(range(biol,"max")),ac(iYr),] + survivors[ac(range(biol,"max"))]
                }
              
              
               
                #- Update fishery to year iYr-1
                landings.n(fishery)[,ac(iYr-1)]     <- sweep(sweep(f[,ac(iYr-1),,,,],c(1:2,4:6),z,"/"),c(1:2,4:6),n(biol)[,ac(iYr-1)]*(1-exp(-z)),"*")
              
                #- Create stock object for assessment
                yrmin1      <- iYr -1
                TaY         <- yrmin1               #Terminal assessment year
                ImY         <- TaY+1                #Intermediate Year
                FcY         <- TaY+2                #Forecast year
              
                idxyrmin1   <- which(dimnames(biol@n)$year == yrmin1)
                stocks@catch.n[,ac(yrmin1)]        <- unitSums(catch.n(fishery)[,ac(yrmin1)])   * ctch[,ac(yrmin1)]
                stocks@landings.n[,ac(yrmin1)]     <- unitSums(landings.n(fishery)[,ac(yrmin1)])* ctch[,ac(yrmin1)]
                stocks@landings.wt[,ac(yrmin1)]    <- stocks@catch.wt[,ac(yrmin1)]
                stocks@discards.n[,ac(yrmin1)]     <- 0
                stocks@discards.wt[,ac(yrmin1)]    <- 0
                stocks@catch[,ac(yrmin1)]          <- computeCatch(stocks)[,ac(yrmin1)]
                stocks@landings[,ac(yrmin1)]       <- computeLandings(stocks)[,ac(yrmin1)]
                stocks@discards[,ac(yrmin1)]       <- computeDiscards(stocks)[,ac(yrmin1)]
              
                #-Do the assessment
                # for the relevant years, apply the assessment errors
                stocks@stock.n[,ac(iYr-1)] <- biol@n[,ac(iYr-1)]     * exp(devN[,1,ac(iYr),,,])
                stocks@harvest[,ac(iYr-1)] <- f[,ac(iYr-1)]          * exp(devF[,1,ac(iYr),,,])
                stocks@stock   <- computeStock(stocks)

                # overwrite the last estimated recruitment?
                 stocks@stock.n[1,ac(TaY)] <-    exp(yearMeans(log(stocks@stock.n[1,ac(1990:(TaY-1))])))
                
                
                f48<-ac(4:8)
                # copy the perception of the stock in terminal assessment year to the stockstore object
                stockstore[,ac(TaY)]@stock.n  <-    stocks[,ac(TaY)]@stock.n
                stockstore[,ac(TaY)]@harvest  <-    sweep(stocks[,ac(TaY)]@harvest * c(fbar[,ac(ImY)]@.Data),c(2:6)  ,quantMeans(stocks@harvest[f48,ac(TaY)]),"/")
                stockstore[,ac(TaY)]@catch.wt  <-   stocks[,ac(TaY)]@catch.wt
                
                
                # survivors for the short term forecast
                 dimnames(survivors)$year<-ac(iYr)
                    
                # recruitment for the first year in the short term forecast is the geometric mean of the historical time series 
                survivors[ac(0),]                  <-   exp(yearMeans(log(stock.n(stocks)[1,ac(1990:(TaY-1))])))
                                                                                         #(rep("!!!! no assessment error implemented for the moment!!!",10))
                #Set plusgroup at 12 (which is true plusgroup - recruitment)
                survivors[-1,]    <- FLQuant(setPlusGroup(stocks[,ac(TaY)]@stock.n * exp(-stocks[,ac(TaY)]@harvest-stocks[,ac(TaY)]@m),11)@.Data,
                                             dimnames=list(age=dimnames(stocks@stock.n)$age[-1],year=ac(TaY),unit=dimnames(stocks@stock.n)$unit,
                                                           season=dimnames(stocks@stock.n)$season,area=dimnames(stocks@stock.n)$area,iter=1:nits))
              
                
                #- Project 4-fleet setup
                f48<-ac(4:8)
                fm<- sweep(stocks[,ac(TaY)]@harvest * c(fbar[,ac(ImY)]@.Data),c(2:6)   ,quantMeans(stocks@harvest[f48,ac(TaY)]),"/")
                TAC[,ac(ImY)] <-  quantSums(catch.wt(stocks)[,ac(ImY)] *  fm/(0.15+fm) * survivors*(1-exp(- (fm+m))))
                   
                 
                #-Calculate effort accordingly (assuming constant catchability)
                f[,ac(ImY)]               <- sweep(catch.sel(fishery)[,ac(ImY)],c(3,6),pmin(maxEff,f31tF(TAC*TACusage,biol,ImY,fishery)),"*")
              
              }
              
              stockstore@landings.n   <- stockstore@harvest * stockstore@stock.n * (1- exp(-stockstore@harvest - stockstore@m)) / (stockstore@harvest + stockstore@m)
              stockstore@landings.wt<-stockstore@catch.wt
                #-------------------------------------------------------------------------------
                # 3): Save the objects
                #-------------------------------------------------------------------------------
              
  
             
             
save(biol          ,file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalbiol.RData",        sep=""))
save(fishery       ,file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalfishery.RData",     sep=""))
save(stockstore    ,file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalpercievedstocks.RData",      sep=""))
save(f             ,file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),nits,"its",prlength,"yrs_Finalf.RData",           sep=""))
             
            
             

