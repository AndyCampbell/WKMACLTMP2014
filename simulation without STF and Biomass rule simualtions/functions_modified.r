#-------------------------------------------------------------------------------
#- Conversions
#-------------------------------------------------------------------------------
an <- function(x){return(as.numeric(x))}
ac <- function(x){return(as.character(x))}



#-------------------------------------------------------------------------------
#- simple plot function for checking input
#-------------------------------------------------------------------------------

plquant<-function(quant)
{
library(lattice)
xyplot(data~year|age,groups=iter,data=quant,type="l")
}

#-------------------------------------------------------------------------------
#- Draw random numbers from a log-normal distribution, but with truncation possibility
#-------------------------------------------------------------------------------
rtlnorm <- function (n, mean = 0, sd = 1, lower = -Inf, upper = Inf)
{
   ret <- numeric()
   if (length(n) > 1)
       n <- length(n)
   while (length(ret) < n) {
       y <- rlnorm(n - length(ret), mean, sd)
       y <- y[y >= lower & y <= upper]
       ret <- c(ret, y)
   }
   stopifnot(length(ret) == n)
   ret
}

#-------------------------------------------------------------------------------
#- Fit vonBertalanffy curve to weight data and draw new values based on sd
#-------------------------------------------------------------------------------
randomWeightBertalanffy <- function(histWeight,ages,histPeriod,projPeriod,iters,estim=list(Winf,K,t0),givenMeanWeight=NULL){

                              #- Make object to store new weights in
                              storeWt           <- array(NA,dim=c(length(ages),length(projPeriod),1,1,1,iters))
                              #- Get historic observed weights
                              wtsdat            <- as.data.frame(histWeight[,histPeriod])
                              colnames(wtsdat)  <- c(colnames(wtsdat)[1:6],"weight")
                              #- Fit von Bertalanffy to data
                              wvBertalanffy     <- nls(weight ~   (Winf*(1-exp(-K*(age-t0)))^3.038),data=wtsdat,start=list(Winf=estim$Winf,K=estim$K,t0=estim$t0),control=list(minFactor=1/2^16,maxiter=1000,warnOnly=T))

                              #- Get mean and sd of fit
                              meanAge <- log(predict(wvBertalanffy,new=data.frame(age=ages)))
                              if(is.null(givenMeanWeight)==F) meanAge <- log(givenMeanWeight)
                              sdAge   <- apply(log(histWeight),1,sd,na.rm=T)@.Data
                              #- Determine 2sd bound
                              bounds  <- cbind(exp(meanAge)-2*apply(histWeight,1,sd)@.Data,
                                               exp(meanAge)+2*apply(histWeight,1,sd)@.Data)
                              bounds[bounds<=0] <- 0.01

                              #- Draw new weights
                              for(iAge in 1:length(ages)){
                                for(iYr in 1:length(projPeriod)){
                                  for(iTer in 1:iters){
                                    set.seed(iAge*iYr*iTer)
                                    if(iAge == 1){ storeWt[iAge,iYr,,,,iTer] <- rtlnorm(1,mean=meanAge[iAge],sd=sdAge[iAge,,,,,],lower=bounds[iAge,1],upper=bounds[iAge,2])
                                    } else {       storeWt[iAge,iYr,,,,iTer] <- rtlnorm(1,mean=meanAge[iAge],sd=sdAge[iAge,,,,,],lower=storeWt[(iAge-1),iYr,,,,iTer],upper=bounds[iAge,2])
                                      }
                                  }
                                }
                              }
                            return(storeWt)}
                            
#-------------------------------------------------------------------------------
#- Update stock objects with 1 year of fisheries and biological data
#-------------------------------------------------------------------------------
updateStocks <- function(iStocks,iFishery,yr,iBiol,iCatchResids){


# iStocks        <- tmp_stocks
# iFishery       <- tmp_fishery
# yr             <- yrmin1
# iBiol          <- tmp_biol
# iCatchResids   <- ctch 
#



                  iStocks@catch.n[,ac(yr)]        <- unitSums(catch.n(iFishery)[,ac(yr)])   * iCatchResids[,ac(yr)]
                  iStocks@landings.n[,ac(yr)]     <- unitSums(landings.n(iFishery)[,ac(yr)])* iCatchResids[,ac(yr)]
                  iStocks@discards.n[,ac(yr)]     <- unitSums(discards.n(iFishery)[,ac(yr)])
                  iStocks@landings.wt[,ac(yr)]    <- unitSums(landings.n(iFishery)[,ac(yr)] * landings.wt(iFishery)[,ac(yr)]) / unitSums(landings.n(iFishery)[,ac(yr)])
                  iStocks@discards.wt[,ac(yr)]    <- 0
                  iStocks@catch.wt[,ac(yr)]       <- iStocks@landings.wt[,ac(yr)]
                  iStocks@catch[,ac(yr)]          <- computeCatch(iStocks)[,ac(yr)]
                  iStocks@landings[,ac(yr)]       <- computeLandings(iStocks)[,ac(yr)]
                  iStocks@discards[,ac(yr)]       <- computeDiscards(iStocks)[,ac(yr)]
                  
                  iStocks@mat[,ac(yr)]            <- iBiol@fec[,ac(yr)]
                  iStocks@stock.wt[,ac(yr)]       <- iBiol@wt[,ac(yr)]
                  iStocks@m[,ac(yr)]              <- iBiol@m[,ac(yr)]
                return(iStocks)}
                
#-------------------------------------------------------------------------------
#- Calculate the harvest by fleet given a TAC for a stock
#-------------------------------------------------------------------------------
fleet.harvest <- function(stk,iYr,TACS){
                    nUnits                          <- dims(stk)$unit
                    nIter                           <- dims(TACS)$iter
                    res                             <- matrix(NA,ncol=nIter,nrow=nUnits,dimnames=list(units=dimnames(stk@stock.n)$unit,iter=1:nIter))
                    for(iTer in 1:nIter) res[,iTer] <- nls.lm(par=rep(1,nUnits),rescaleF,stk=iter(stk,iTer),iYr=iYr,TACS=c(iter(TACS,iTer)),nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL)$par
                    stk@harvest[,iYr]               <- sweep(stk@harvest[,iYr],c(3,6),res,"*")
                 return(stk@harvest[,iYr])}

  #faster version   ImY
fleet.harvestF<- function(stk,iYr,TACS){
                    nUnits                          <- dims(stk)$unit
                    nIter                           <- dims(TACS)$iter
                    res                             <- matrix(NA,ncol=nIter,nrow=nUnits,dimnames=list(units=dimnames(stk@stock.n)$unit,iter=1:nIter))
                    for(iTer in 1:nIter){
                      Ns      = stk@stock.n[,iYr,1,,,iTer]@.Data
                      Fs      = stk@harvest[,iYr,,,,iTer]@.Data
                      Cwts    = stk@catch.wt[,iYr,,,,iTer]@.Data
                      Ms      = stk@m[,iYr,1,,,iTer]@.Data
                      res[,iTer] <- nls.lm(par=rep(1,nUnits),rescaleFF,Ns=Ns,Fs=Fs,Cwts=Cwts,Ms=Ms,TACS=TACS[,,,,,iTer],nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL)$par
                    }
                    stk@harvest[,iYr]               <- sweep(stk@harvest[,iYr],c(3,6),res,"*")
                 return(stk@harvest[,iYr])}
                 
  #faster version
fleet.harvestFF<- function(stk,iYrr,TACS){
# stk <- stf ;  iYrr<-  ImY  ; TACS <- iTAC[,ImY]

                    lin     <- substr(R.Version()$os,1,3)== "lin"
                    nUnits                          <- dims(stk)$unit
                    nIter                           <- dims(TACS)$iter
                    res                             <- matrix(NA,ncol=nIter,nrow=nUnits,dimnames=list(units=dimnames(stk@stock.n)$unit,iter=1:nIter))
                    for(iTer in 1:nIter){
                      Ns      = c(stk@stock.n[,iYrr,1,,,iTer]@.Data)
                      Fs      = stk@harvest[,iYrr,,,,iTer,drop=T]@.Data
                      Cwts    = stk@catch.wt[,iYrr,,,,iTer,drop=T]@.Data
                      Ms      = c(stk@m[,iYrr,1,,,iTer]@.Data)
                      if(lin) res[,iTer] <- nls.lm(par=rep(1,nUnits),rescaleFFF,Ns=Ns,Fs=Fs,Cwts=Cwts,Ms=Ms,TACS=c(TACS[,,,,,iTer]),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits),nls.lm.control(ftol = (.Machine$double.eps)))$par
                      if(!lin)res[,iTer] <- nls.lm(par=rep(1,nUnits),rescaleFFF,Ns=Ns,Fs=Fs,Cwts=Cwts,Ms=Ms,TACS=c(TACS[,,,,,iTer]),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits),nls.lm.control(ftol = (.Machine$double.eps)))$par
                    }
                    stk@harvest[,iYrr]               <- sweep(stk@harvest[,iYrr],c(3,6),res,"*")
                 return(stk@harvest[,iYrr])}
                 
#-------------------------------------------------------------------------------
#- Function to scale the F pattern for stock
#-------------------------------------------------------------------------------
rescaleF      <- function(mult,stk.=stk,iYr.=iYr,TACS.=TACS){
                    stk.@harvest[,iYr.] <- sweep(stk.@harvest[,iYr.],3,mult,"*")
                    stkZ                <- unitSums(stk.@harvest[,iYr.]) + stk.@m[,iYr.,1]
                    res                 <- sqrt(c((TACS. - c(apply(sweep(stk.@stock.n[,iYr.] * stk.@catch.wt[,iYr.] * sweep(stk.@harvest[,iYr.],c(1:2,4:6),stkZ,"/"),c(1:2,4:6),(1-exp(-stkZ)),"*"),3:6,sum,na.rm=T)))^2))
                 return(res)}

  #Faster version
rescaleFF     <- function(mult,Ns=Ns,Fs=Fs,Cwts=Cwts,Ms=Ms,TACS.=TACS){
                    Fs    <- sweep(Fs,3,mult,"*")
                    stkZ  <- apply(Fs,1,sum) + Ms
                    res   <- sqrt(c((TACS. - c(apply(sweep(sweep(Fs,c(1:2,4:6),stkZ,"/")* Cwts,c(1:2,4:6),Ns * (1-exp(-stkZ)),"*") ,3:6,sum))))^2)
                 return(res)}
                 
  #Faster version
rescaleFFF     <- function(mult,Ns.=Ns,Fs.=Fs,Cwts.=Cwts,Ms.=Ms,TACS.=TACS){
                    Fs.    <- t(t(Fs.) * mult)
                    stkZ  <- rowSums(Fs.) + Ms.
                    res   <- sqrt(c((TACS. - c(colSums(sweep(sweep(Fs.,1,stkZ,"/")* Cwts.,1,Ns. * (1-exp(-stkZ)),"*")))))^2)
                 return(res)}
                 
#-------------------------------------------------------------------------------
#- Calculate the catch by fleet
#-------------------------------------------------------------------------------
harvestCatch  <-  function(stk.,iYr){
                    stkZ      <- unitSums(stk.@harvest[,iYr]) + stk.@m[,iYr,1]
                    res       <- apply(sweep(stk.@stock.n[,iYr] * stk.@catch.wt[,iYr] * sweep(stk.@harvest[,iYr],c(1:2,4:6),stkZ,"/"),c(1:2,4:6),(1-exp(-stkZ)),"*"),3:6,sum,na.rm=T)
                  return(res)}
                  
#-------------------------------------------------------------------------------
#- Management plan: calculate A and B TAC
#-------------------------------------------------------------------------------

                 
find.FABF     <- function(mult,stk.,mpPoints.){

#mult<-1
#stk.<- iter(stf[,FcY],iTer)
#mpPoints.<-mpPoints
                    Fs      <- sweep(stk.@harvest@.Data,3,mult,"*")
                    Ns      <- stk.@stock.n@.Data[,,1,,,]
                    Ms      <- stk.@m@.Data[,,1,,,]
                    Hspwns  <- stk.@harvest.spwn@.Data[,,1,,,]
                    Mspwns  <- stk.@m.spwn@.Data[,,1,,,]
                    Mats    <- stk.@mat@.Data[,,1,,,]
                    Swghts  <- stk.@stock.wt@.Data[,,1,,,]
                    Cwghts  <- stk.@catch.wt@.Data

                    bigF              <- apply(Fs,1,sum)
                    ssb               <- sum(Ns * Swghts * exp(-bigF*Hspwns - Ms*Mspwns) * Mats)
                    if(ssb < mpPoints.$Btrigger) resA <- mpPoints.$Ftarget* ssb / mpPoints.$Btrigger
                    if(ssb > mpPoints.$Btrigger) resA <- mpPoints.$Ftarget   
                                 
                    fbarA     <- mean(bigF[ac(4:8)])
                    
                    ret       <- (fbarA-resA)^2
               
                 return(ret)}

#-------------------------------------------------------------------------------
#- Function to calculate the change in effort
#-------------------------------------------------------------------------------
f21t <- function(fmult,Q1,morts,Fsq,landsel,nums,land_wts){
    bigF <- Fsq * fmult
    return(abs(Q1 - sum( ((bigF)/(bigF + morts))*nums*(1-exp(-bigF-morts))* land_wts)))
  }

#-------------------------------------------------------------------------------
#- Function to return the change in effort
#-------------------------------------------------------------------------------
f31t <- function(Q1,iBiol,ImY,iFishery,iFailRuns=NULL){
  lan_wt  <- landings.wt( iFishery)[,ac(ImY)]
  fpattern<- catch.sel(iFishery)[,ac(ImY)]
  iters   <- dims(Q1)$iter
  nUnits  <- dims(Q1)$unit
  fm      <- matrix(NA,nrow=dims(iFishery)$unit,ncol=iters,dimnames=list(dimnames(iFishery@landings.n)$unit,dimnames(iFishery@landings.n)$iter))
  if(is.null(iFailRuns)==T){
    for(i in 1:iters){
      iMort    <- iBiol@m[,ac(ImY),,,,i]@.Data
      iNum     <- iBiol@n[,ac(ImY),,,,i]@.Data
      iFsq     <- fpattern[,,,,,i]@.Data
      iLanwt   <- lan_wt[,,,,,i]@.Data
      iQ1      <- Q1[,ac(ImY),,,,i]@.Data
      fm[,i]   <- nls.lm(par=rep(1,nUnits),rescaleFBiol,iFsq.=iFsq,iMort.=iMort,iNum.=iNum,iLanwt.=iLanwt,TACS.=iQ1,nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits))$par
    }
  } else {
    for(i in (1:iters)[-iFailRuns]){
      iMort    <- iBiol@m[,ac(ImY),,,,i]@.Data
      iNum     <- iBiol@n[,ac(ImY),,,,i]@.Data
      iFsq     <- fpattern[,,,,,i]@.Data
      iLanwt   <- lan_wt[,,,,,i]@.Data
      iQ1      <- Q1[,ac(ImY),,,,i]@.Data
      fm[,i]   <- nls.lm(par=rep(1,nUnits),rescaleFBiol,iFsq.=iFsq,iMort.=iMort,iNum.=iNum,iLanwt.=iLanwt,TACS.=iQ1,nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits))$par
    }
  }
  return(fm)
}






f31tF <- function(Q1,iBiol,ImY,iFishery,iFailRuns=NULL){

 #  Q1=TAC*TACusage ; iBiol=biol ; iFishery =fishery


  lin     <- substr(R.Version()$os,1,3)== "lin"
  lan_wt  <- landings.wt( iFishery)[,ac(ImY)]
  fpattern<- catch.sel(iFishery)[,ac(ImY)]
  iters   <- dims(Q1)$iter
  nUnits  <- dims(Q1)$unit
  fm      <- matrix(NA,nrow=dims(iFishery)$unit,ncol=iters,dimnames=list(dimnames(iFishery@landings.n)$unit,dimnames(iFishery@landings.n)$iter))
  if(is.null(iFailRuns)==T){
    for(i in 1:iters){
      iMort    <- c(iBiol@m[,ac(ImY),,,,i]@.Data)
      iNum     <- c(iBiol@n[,ac(ImY),,,,i]@.Data)
      iNum[1] <- c(iBiol@n[1,ac(ImY-1),,,,i]@.Data)
      iFsq     <- fpattern[,,,,,i,drop=T]@.Data
      iLanwt   <- lan_wt[,,,,,i,drop=T]@.Data
      iQ1      <- c(Q1[,ac(ImY),,,,i]@.Data)
      if(lin) fm[,i]   <- nls.lm(par=rep(1,nUnits),rescaleFBiolF,iFsq.=iFsq,iMort.=iMort,iNum.=iNum,iLanwt.=iLanwt,TACS.=iQ1,nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits))$par
      if(!lin)fm[,i]   <- nls.lm(par=rep(1,nUnits),rescaleFBiolF,iFsq.=iFsq,iMort.=iMort,iNum.=iNum,iLanwt.=iLanwt,TACS.=iQ1,nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits))$par
    }
  } else {
    for(i in (1:iters)[-iFailRuns]){
      iMort    <- c(iBiol@m[,ac(ImY),,,,i]@.Data)
      iNum     <- c(iBiol@n[,ac(ImY),,,,i]@.Data)
      iNum[1] <- c(iBiol@n[1,ac(ImY-1),,,,i]@.Data)
      iFsq     <- fpattern[,,,,,i,drop=T]@.Data
      iLanwt   <- lan_wt[,,,,,i,drop=T]@.Data
      iQ1      <- c(Q1[,ac(ImY),,,,i]@.Data)
      if(lin) fm[,i]   <- nls.lm(par=rep(1,nUnits),rescaleFBiolF,iFsq.=iFsq,iMort.=iMort,iNum.=iNum,iLanwt.=iLanwt,TACS.=iQ1,nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits))$par
      if(!lin)fm[,i]   <- nls.lm(par=rep(1,nUnits),rescaleFBiolF,iFsq.=iFsq,iMort.=iMort,iNum.=iNum,iLanwt.=iLanwt,TACS.=iQ1,nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits))$par
    }
  }
  return(fm)
}


#-------------------------------------------------------------------------------
#- Function to scale the F pattern for biol
#-------------------------------------------------------------------------------
rescaleFBiol <- function(mult,iFsq.,iMort.,iNum.,iLanwt.,TACS.){
                    iF          <- sweep(iFsq.,3,mult,"*")
                    stkZ        <- apply(iF,1,sum) + iMort.
                    res         <- sqrt(c((c(TACS.) - c(apply(sweep(sweep(iLanwt. * iF,c(1:2,4:6),stkZ,"/"),c(1:2,4:6),iNum. * (1-exp(-stkZ)),"*"),3:6,sum,na.rm=T)))^2))
                 return(res)}

rescaleFBiolF <- function(mult,iFsq.,iMort.,iNum.,iLanwt.,TACS.){
                    iF          <- t(t(iFsq.) * mult)
                    stkZ        <- rowSums(iF) + iMort.
                    res         <- sqrt(c((TACS. -    c(colSums(sweep(sweep(iF,1,stkZ,"/")* iLanwt.,1,iNum. * (1-exp(-stkZ)),"*")))))^2)
                 return(res)}


#-------------------------------------------------------------------------------
#- Function to calculate the ssb based on biological numbers at age
#-------------------------------------------------------------------------------

ssbb <- function(iBiol,iHarvest,iStck){

#  iBiol= biol[,ac(yrs)] ; iHarvest =   ;iStck=
          res <- quantSums(iBiol@n * iBiol@wt * exp(-iHarvest * iStck@harvest.spwn - iBiol@m * iStck@m.spwn) * iBiol@fec)
        return(res)}

#-------------------------------------------------------------------------------
#- Function to estimate the mean and sd of the lognormal distribution that fits
#  the observed recruitment best
#-------------------------------------------------------------------------------

optimRecDistri <- function(pars,recs){
                    sim <- plnorm(recs,meanlog=pars[1],sdlog=pars[2])
                    lik <- sum(abs(log(sim) - log(seq(0,1,length.out=length(recs)+2)[-c(1,length(recs)+2)])))
                  return(lik)}

#-------------------------------------------------------------------------------
#- Function for paralelisation
#-------------------------------------------------------------------------------

startCluster <- function(snow.cores){
                  #-Start up cluster
                  library(snow)
                  library(Rmpi)
                  num.cores       <- snow.cores
                  snow.cluster    <- c(rep("localhost",num.cores))
                  cl              <- makeMPIcluster(count=length(snow.cluster))
                  options(timeout=14*86400)
                  clusterEvalQ(cl,library(FLSAM))
                  clusterEvalQ(cl,library(MASS))
                  clusterEvalQ(cl,library(msm))
                return(cl)}
                
exportCluster <- function(cl,lst){
                  clusterExport(cl,lst)
                 }

retroLB                     <- function(x){
                                cat(paste("\n\n\n\ THIS IS ITERATION NUMBER",x,"\n\n\n"))
                                stck  <- window(iter(stocks,x),end=histMaxYr)
                                tun   <- FLIndices()
                                for(iSurv in names(surveys))
                                  tun[[iSurv]] <- iter(window(surveys[[iSurv]],end=range(NSH.tun[[iSurv]])["maxyear"]),x)
                                res <- retro(stock=stck,indices=tun,
                                             control=ctrl,retro=10,base.assess=base.assess)
                               return(res)}

#-------------------------------------------------------------------------------
#- Function to overwrite the updated tmp_stocks to stocks (but do not overwrite
#   the elements listed below from tmp to stocks, hence they are taken from
#   the stocks object)
#-------------------------------------------------------------------------------


tmp2stocks <- function(stocks2,tmp_stocks2,yr){

  
                TaYtmp_stocks         <- tmp_stocks2[,ac(yr)]
                TaYstocks             <- stocks2[,ac(yr)]
                TaYtmp_stocks@stock   <- TaYstocks@stock
                TaYtmp_stocks@stock.n <- TaYstocks@stock.n
                TaYtmp_stocks@harvest <- TaYstocks@harvest
                TaYtmp_stocks@name    <- TaYstocks@name
                TaYtmp_stocks@desc    <- TaYstocks@desc
                TaYtmp_stocks@range   <- TaYstocks@range
                stocks2[,ac(yr)]       <- TaYtmp_stocks[,ac(yr)]
              return(stocks2)}
              
              
              
#-----------------------------------------------------------------------------------------------------------
# Recruitment function => gives the R accoring to the SSB and SR parameters
#-----------------------------------------------------------------------------------------------------------

B2R<-function (Biom, SRp,ite,yr,recreg){

#Biom<- ssb(biol)[,ac(iYr-1),,,,its]
#SRp<-  SRmod[its,]
#ite<-  its
#yr<-   iYr
#recreg<- RecType
#
      if( recreg  !=  "SGRec"){
          mod <-  SRp$mod
          A     <-  SRp$A
          B     <-  SRp$B
          sig   <-  SRp$sigma  # this rescaling should have been done in the Bayesian protocol

          if (mod=="segreg" ) mu <-  A*(Biom>=B)+A*Biom*(Biom<B)/B
          if (mod=="ricker" ) mu <-  A*Biom*exp(-B*Biom)
          if (mod=="bevholt") mu <-  A*Biom/(B+Biom)
          
          REC<-exp(log(mu)+rnorm(1,0,sig)) 
          
          REC<- max(REC,10)
        }
        if ( recreg == "SGRec")
          REC   <- SRp[1,yr]

    return(REC)}              


    
B2Rdet<-function (Biom, SRp,ite,yr,recreg){

#Biom<- ssb(biol)[,ac(iYr-1),,,,its]
#SRp<-  SRmod[its,]
#ite<-  its
#yr<-   iYr
#recreg<- RecType
#
      if( recreg  !=  "SGRec"){
          mod <-  SRp$mod
          A     <-  SRp$A
          B     <-  SRp$B
          sig   <-  SRp$sigma  # this rescaling should have been done in the Bayesian protocol

          if (mod=="segreg" ) mu <-  A*(Biom>=B)+A*Biom*(Biom<B)/B
          if (mod=="ricker" ) mu <-  A*Biom*exp(-B*Biom)
          if (mod=="bevholt") mu <-  A*Biom/(B+Biom)
          
          REC<-mu
          
          REC<- max(REC,10)
        }
        if ( recreg == "SGRec")
          REC   <- SRp[1,yr]

    return(REC)}      



    
B2Rbis<-function (Biom, SRp,its,yr,dev){

# Biom  <-  ssb(biol)[,ac(iYr-1),,,,its]
# SRp   <-  SRmod[its,]
# ite   <-  its
# dev   <-  devR
# yr    <-  iYr-1

          mod   <-  SRp$mod
          A     <-  SRp$A
          B     <-  SRp$B
          sig   <-  SRp$sigma  # this rescaling should have been done in the Bayesian protocol

          if (mod=="segreg" ) mu <-  A*(Biom>=B)+A*Biom*(Biom<B)/B
          if (mod=="ricker" ) mu <-  A*Biom*exp(-B*Biom)
          if (mod=="bevholt") mu <-  A*Biom/(B+Biom)
          
          REC<- mu * exp(dev[,ac(yr),,,,its])
          
          REC<- max(REC,10)
     

    return(REC)}              




#-----------------------------------------------------------------------------------------------------------
# modified mcmc function to resample stock from the SAM assessment.
#-----------------------------------------------------------------------------------------------------------
    
monteCarloStock2 <- 
function (stck, sam, realisations, run.dir = tempdir(),seed_number) 
{
  require(MASS)
  ctrl <- sam@control
  mcstck <- propagate(stck, iter = realisations)
  mcstck <- window(mcstck, start = range(sam)["minyear"], end = range(sam)["maxyear"])
  mcstck@stock.n[] <- NA
  mcstck@harvest[] <- NA
  set.seed(seed_number)
  random.param <- mvrnorm(realisations, sam@params$value, sam@vcov)
  save(random.param, file = file.path(run.dir, "random.param.RData"))
  n.states <- length(unique(ctrl@states[names(which(ctrl@fleets == 
                                                      0)), ]))
  yrs <- dims(sam)$minyear:dims(sam)$maxyear
  ages <- dims(sam)$age
  u <- random.param[, which(colnames(random.param) == "U")]
  ca <- random.param[, which(colnames(random.param) == "logCatch")]
  idxNs <- c(mapply(seq, from = seq(1, ncol(u), ages + n.states), 
                    to = seq(1, ncol(u), ages + n.states) + ages - 1, by = 1))
  idxFs <- c(mapply(seq, from = seq(1, ncol(u), ages + n.states) + 
                      ages, to = seq(1, ncol(u), ages + n.states) + n.states + 
                      ages - 1, by = 1))
  mcstck@stock.n[] <- exp(aperm(array(u[, idxNs], dim = c(realisations, 
                                                          ages, length(yrs))), perm = c(2, 3, 1)))
  mcstck@harvest[] <- exp(aperm(array(u[, idxFs], dim = c(realisations, 
                                                          n.states, length(yrs))), perm = c(2, 3, 1)))[ctrl@states[names(which(ctrl@fleets == 
                                                                                                                                 0)), ], , ]
  mcstck@catch[] <- exp(aperm(array(ca, dim = c(realisations, 
                                                length(yrs))), perm = c(2, 1)))
  return(mcstck)
}

#stk1 <- monteCarloStock2(Mac,Mac.sam,nits,run.dir="tmp",seed_number = 555)
#stk2 <- monteCarloStock2(Mac,Mac.sam,nits,run.dir="tmp",seed_number = 555)
#identical(stk1,stk2)
#



#-----------------------------------------------------------------------------------------------------------
# Calculating TAC with different HCR methods
#-----------------------------------------------------------------------------------------------------------

# stocks <- stocks[,1:idxyrmin1]
# fishery <- window(fishery,histMinYr,yrmin1)
# scen <- mpPoints$scen

HCR.TAC <- function(HCR.method,stocks,survivors,fishery,iYr,TAC,scen,histMaxYr,mpPoints,mpOptions){
  switch(HCR.method,
         stf         =projectMac(stocks,survivors,fishery,iYr,TAC,scen,NULL,histMaxYr,mpPoints,mpOptions),
         without.stf = HCR.without.stf(stocks,survivors,fishery,iYr,TAC,scen,NULL,histMaxYr,mpPoints,mpOptions),
         fix.prop    = HCR.fix.prop(stocks, survivors,fishery,iYr, TAC, scen, NULL, histMaxYr,mpPoints,mpOptions))
}

#-----------------------------------------------------------------------------------------------------------
# Calculating TAC with HCR without SHORT TERM FORECAST
#-----------------------------------------------------------------------------------------------------------



HCR.without.stf <- function(iStocks,iSurvivors,iFishery,iYr,iTAC,iScenario,iFailRuns=NULL,iHistMaxYr,mpPoints,mpoptions){
  
  iStocks<-stocks[,1:idxyrmin1]
  iSurvivors<-survivors                                 
  iFishery<-window(fishery,histMinYr,yrmin1)
  iTAC<-TAC
  iScenario<-scen
  iFailRuns=NULL 
  iHistMaxYr<-histMaxYr
  mpoptions<-mpOptions
   
  require(minpack.lm)
  
  lin     <- substr(R.Version()$os,1,3)== "lin"
  stk     <- iStocks
  stk.sur <- iSurvivors
  #===============================================================================
  # Setup control file
  #===============================================================================
  
  DtY         <- ac(range(stk)["maxyear"]) #Data year
  ImY         <- ac(an(DtY)+1) #Intermediate year
  FcY         <- ac(an(DtY)+2) #Forecast year
  CtY         <- ac(an(DtY)+3) #Continuation year
  CtY1        <- ac(an(DtY)+4)
  FuY         <- c(ImY,FcY,CtY,CtY1)#Future years
  pyears      <- an(DtY) - iHistMaxYr #Years away from original historic max year
  
  #- We substract pyears from the starting point, because with every year forward, we move the year ranges one forward too.
  RECS        <- FLQuants("ImY"=iSurvivors[1,],"FcY"=exp(apply(log(rec(stk)[,ac((range(stk)["maxyear"]-12-pyears):(range(stk)["maxyear"]))]),3:6,mean,na.rm=T)),
                          "CtY"=exp(apply(log(rec(stk)[,ac((range(stk)["maxyear"]-12-pyears):(range(stk)["maxyear"]))]),3:6,mean,na.rm=T)))
  
  yrs1        <- list("m.spwn","harvest.spwn","stock.wt")
  yrs3        <- list("mat")
  yrs5        <- list("m")
  
  dsc         <- "NEa mackerel  "
  nam         <- "Mac"
  dms         <- dimnames(stk@m)
  dms$year    <- c(rev(rev(dms$year)[1:3]),ImY,FcY,CtY,CtY1)
  dms$unit    <- dimnames(iFishery@landings.n)$unit
  
  
  f48         <- ac(4:8)
  
  #===============================================================================
  # Setup stock file
  #===============================================================================
  
  stf         <- FLStock(name=nam,desc=dsc,m=FLQuant(NA,dimnames=dms))
  for(i in dms$unit) stf[,,i] <- window(stk,start=an(dms$year)[1],end=rev(an(dms$year))[1])
  units(stf)  <- units(stk)
  # Fill slots that are the same for all fleets
  for(i in c(unlist(yrs1),unlist(yrs3),unlist(yrs5))){
    if(i %in% unlist(yrs1)){ for(j in dms$unit){ slot(stf,i)[,FuY,j] <- slot(stk,i)[,DtY]}}
    if(i %in% unlist(yrs3)){ for(j in dms$unit){ slot(stf,i)[,FuY,j] <- yearMeans(slot(stk,i)[,ac((an(DtY)-2):an(DtY))])}}
    if(i %in% unlist(yrs5)){ for(j in dms$unit){ slot(stf,i)[,FuY,j] <- yearMeans(slot(stk,i)[,ac((an(DtY)-4):an(DtY))])}}
  }
  
  # Fill slots that are unique for the fleets
  for(iFuY in FuY){
    stf@harvest[,iFuY]       <- stk@harvest[,DtY]
    stf@catch.wt[,iFuY]      <- iFishery@landings.wt[,ac(DtY)]
    stf@landings.wt[,iFuY]   <- iFishery@landings.wt[,ac(DtY)]
  }
  
  # Fill slots that have no meaning for NSAS
  stf@discards.n[]        <- 0
  stf@discards[]          <- 0
  stf@discards.wt[]       <- 0
  stf@m[]                 <- 0.15
  #===============================================================================
  # Intermediate year
  #===============================================================================
  
  for(i in dms$unit)        stf@stock.n[,ImY,i]     <- stk.sur
  stf@harvest[,ImY]         <- fleet.harvestFF(stk=stf,iYrr=ImY,TACS=iTAC[,ImY])
  for(i in dms$unit){
    stf@catch.n[,ImY,i]     <- stf@stock.n[,ImY,i]*(1-exp(-unitSums(stf@harvest[,ImY])-stf@m[,ImY,i]))*(stf@harvest[,ImY,i]/(unitSums(stf@harvest[,ImY])+stf@m[,ImY,i]))
    stf@catch[,ImY,i]       <- computeCatch(stf[,ImY,i])
    stf@landings.n[,ImY,i]  <- stf@catch.n[,ImY,i]
    stf@landings[,ImY,i]    <- computeLandings(stf[,ImY,i])
  }
  
  #===============================================================================
  # Forecast year
  #===============================================================================
  
  for(i in dms$unit) stf[,FcY,i]    <- stf[,ImY,i] 
  #  for(i in dms$unit) stf@stock.n[2:(dims(stf)$age-1),FcY,i]  <- (stf@stock.n[,ImY,1]*exp(-unitSums(stf@harvest[,ImY])-stf@m[,ImY,1]))[ac(range(stf)["min"]:(range(stf)["max"]-2)),]
  #  for(i in dms$unit) stf@stock.n[dims(stf)$age,FcY,i]        <- apply((stf@stock.n[,ImY,1]*exp(-unitSums(stf@harvest[,ImY])-stf@m[,ImY,1]))[ac((range(stf)["max"]-1):range(stf)["max"]),],2:6,sum,na.rm=T)
  
  ###--- Management options ---###
  
  #- First calculate the TAC according to the HCR (not the LTMP)
  res <- matrix(NA,nrow=dims(stf)$unit,ncol=dims(stf)$iter,dimnames=list(dimnames(stf@stock.n)$unit,dimnames(stf@stock.n)$iter))
  if(is.null(iFailRuns) == T){
    for(iTer in 1:dims(stf)$iter){
      if(lin)  res[,iTer]            <- optimize(find.FABF,stk.=iter(stf[,FcY],iTer), mpPoints.=mpPoints,
                                                 lower=rep(0,dims(stf)$unit),upper=rep(1e5,dims(stf)$unit))$minimum
      if(!lin) res[,iTer]            <- optimize(find.FABF,stk.=iter(stf[,FcY],iTer), mpPoints.=mpPoints,
                                                 lower=rep(0,dims(stf)$unit),upper=rep(1e5,dims(stf)$unit))$minimum
    }
    
  } else {
    for(iTer in (1:dims(stf)$iter)[-iFailRuns]){
      if(lin)  res[,iTer]            <- optimize(find.FABF,stk.=iter(stf[,FcY],iTer), mpPoints.=mpPoints,
                                                 lower=rep(0,dims(stf)$unit),upper=rep(1e5,dims(stf)$unit))$minimum
      if(!lin) res[,iTer]            <- optimize(find.FABF,stk.=iter(stf[,FcY],iTer), mpPoints.=mpPoints,
                                                 lower=rep(0,dims(stf)$unit),upper=rep(1e5,dims(stf)$unit))$minimum
    }
  }
  
  stf@harvest[,FcY]         <- sweep(stf@harvest[,FcY],c(3,6),res,"*")
  for(i in dms$unit){
    stf@catch.n[,FcY,i]     <- stf@stock.n[,FcY,i]*(1-exp(-unitSums(stf@harvest[,FcY])-stf@m[,FcY,i]))*(stf@harvest[,FcY,i]/(unitSums(stf@harvest[,FcY])+stf@m[,FcY,i]))
    stf@catch[,FcY,i]       <- computeCatch(stf[,FcY,i])
    stf@landings.n[,FcY,i]  <- stf@catch.n[,FcY,i]
    stf@landings[,FcY,i]    <- computeLandings(stf[,FcY,i])
  }
  HCRTAC                    <- harvestCatch(stf,FcY)
  SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
  HCRSSB                    <- SSB
  F48                       <- apply(unitSums(stf@harvest[f48,FcY]),2:6,mean)
  
  print(round(c(apply(F48,1,mean)),5))
  
  
  #-----------------------------------------------------------------------------
  #- Scenario no IAV
  #-----------------------------------------------------------------------------
  if(mpoptions$TACvarlim == F | ImY=="2014" ){
    #- No change to TAC calculation
    iTAC[,FcY]              <- HCRTAC
    
  }
  
  #-----------------------------------------------------------------------------
  #- Scenario current Long Term Management Plan
  #-----------------------------------------------------------------------------
  if(mpoptions$TACvarlim == T & ImY!="2014"){
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  #-----------------------------------------------------------------------------
  #- Scenario: if HCR TAC results in F outside of 15% of target F, remove TAC contraint, else use it
  #-----------------------------------------------------------------------------
  if(mpoptions$Fvarlim == T){
    tmpTAC                    <- iTAC[,FcY,]
    tmpTAC[]                  <- HCRTAC
    iTAC[,FcY]                <- HCRTAC
    #- First calculate what the F would be if the TAC would be constraint by 15%
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")])
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")])
    if(length(bidx)>0) tmpTAC[,,c("A"),,,bidx]   <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) tmpTAC[,,c("A"),,,sidx]   <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=tmpTAC)
    
    for(i in dms$unit){
      stf@catch.n[,FcY,i]     <- stf@stock.n[,FcY,i]*(1-exp(-unitSums(stf@harvest[,FcY])-stf@m[,FcY,i]))*(stf@harvest[,FcY,i]/(unitSums(stf@harvest[,FcY])+stf@m[,FcY,i]))
      stf@catch[,FcY,i]       <- computeCatch(stf[,FcY,i])
      stf@landings.n[,FcY,i]  <- stf@catch.n[,FcY,i]
      stf@landings[,FcY,i]    <- computeLandings(stf[,FcY,i])
    }
    #- Determine which iters are outside the +/- 15% from target F or inside this range
    outidx                    <- which(apply(unitSums(stf@harvest[f48,FcY]),2:6,mean,na.rm=T) <  0.9 * F48 |
                                         apply(unitSums(stf@harvest[f48,FcY]),2:6,mean,na.rm=T) >  1.1 * F48)
    inidx                     <- which(apply(unitSums(stf@harvest[f48,FcY]),2:6,mean,na.rm=T) >= 0.9 * F48 &
                                         apply(unitSums(stf@harvest[f48,FcY]),2:6,mean,na.rm=T) <= 1.1 * F48)
    #- If outside this range, do not apply TAC constraint and use 15% cap on F
    if(length(outidx)>0){
      sidx                    <- which((apply(unitSums(stf@harvest[f48,FcY,,,,outidx]),2:6,mean,na.rm=T)) <  0.9 * F48[,,,,,outidx])
      bidx                    <- which((apply(unitSums(stf@harvest[f48,FcY,,,,outidx]),2:6,mean,na.rm=T)) >  1.1 * F48[,,,,,outidx])
      
      if(length(sidx)>0) stf@harvest[,FcY,"A",,,outidx[sidx]]       <- sweep(stf@harvest[,FcY,"A",,,outidx[sidx]],2:6,
                                                                             (F48[,,,,,outidx[sidx]] * 0.9)  /
                                                                               apply(unitSums(stf@harvest[f48,FcY,"A",,,outidx[sidx]]),2:6,mean,na.rm=T),"*")
      if(length(bidx)>0) stf@harvest[,FcY,"A",,,outidx[bidx]]       <- sweep(stf@harvest[,FcY,"A",,,outidx[bidx]],2:6,
                                                                             (F48[,,,,,outidx[bidx]] * 1.1)  /
                                                                               apply(unitSums(stf@harvest[f48,FcY,"A",,,outidx[bidx]]),2:6,mean,na.rm=T),"*")
      
      iTAC[,FcY,c("A"),,,outidx]                <-  quantSums(stf@catch.wt[,FcY,"A",,,outidx] * stf@stock.n[,FcY,"A",,,outidx] *
                                                                (1-exp(-unitSums(stf@harvest[,FcY,,,,outidx])-stf@m[,FcY,"A",,,outidx])) *
                                                                (stf@harvest[,FcY,"A",,,outidx]/(unitSums(stf@harvest[,FcY,,,,outidx])+stf@m[,FcY,"A",,,outidx])))
    }
    if(length(inidx)>0){
      bidx                    <- which(HCRTAC[,,c("A"),,,inidx] > 1.20*iTAC[,ImY,c("A"),,,inidx])
      sidx                    <- which(HCRTAC[,,c("A"),,,inidx] < 0.80*iTAC[,ImY,c("A"),,,inidx])
      if(length(bidx)>0) iTAC[,FcY,c("A"),,,inidx[bidx]] <- 1.20*iTAC[,ImY,c("A"),,,inidx[bidx]]
      if(length(sidx)>0) iTAC[,FcY,c("A"),,,inidx[sidx]] <- 0.80*iTAC[,ImY,c("A"),,,inidx[sidx]]
    }
    #- If SSB below Blim, then no cap at all
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, bank, other years repay + borrow
  #-----------------------------------------------------------------------------
  if(mpoptions$BBscen == "BB"){
    
    #- First HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    #- Second banking and borrowing
    if(iYr == "2013"){             #Bank
      iTAC[,FcY,c("A")]         <- 0.9 * iTAC[,FcY,"A"]
    }
    if(iYr == "2014"){             #Repay banked part                               Borrow
      iTAC[,FcY,c("A")]         <- (iTAC[,ImY,c("A")] / 0.9 - iTAC[,ImY,c("A")]) + (1.1 * iTAC[,FcY,c("A")])
    }
    if(an(iYr) > 2014){            #Repay borrowed part                             Borrow
      iTAC[,FcY,c("A")]         <- (iTAC[,ImY,c("A")] / 1.1 - iTAC[,ImY,c("A")]) + (1.1 * iTAC[,FcY,c("A")])
    }
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, bank, other years : bank
  #-----------------------------------------------------------------------------
  if(mpoptions$BBscen == "Banking"){
    
    #-  HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    #- first year of  banking
    if(iYr == "2013"){             #Bank
      iTAC[,FcY,c("A")]         <- 0.9 * iTAC[,FcY,"A"]
    }
    if(an(iYr) >= 2013){             #Repay banked part                               Bank
      iTAC[,FcY,c("A")]         <- (0.1 * iTAC[,ImY,c("A")] ) + (0.9 * iTAC[,FcY,c("A")])
    }
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, borrow, other years : borrow
  #-----------------------------------------------------------------------------
  if(mpoptions$BBscen == "Borrowing"){
    
    #-  HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    #- first year of  banking
    if(iYr == "2013"){             #Borrow
      iTAC[,FcY,c("A")]         <- 1.1 * iTAC[,FcY,"A"]
    }
    if(an(iYr) >= 2013){             #Repay                   and          Borrow again
      iTAC[,FcY,c("A")]         <- (- 0.1 * iTAC[,ImY,c("A")] ) + (1.1 * iTAC[,FcY,c("A")])
    }
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, borrow, then alternate bank and borrow
  #-----------------------------------------------------------------------------
  if(mpoptions$BBscen == "AlternateBorrow"){
    
    #-  HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    #- first year of  banking
    
    if(iYr == "2013"){             #Borrow
      iTAC[,FcY,c("A")]         <- 1.1 * iTAC[,FcY,"A"]
    }
    if(an(iYr) >= 2013){ 
      
      al<-(an(iYr)%%2-0.5)/5          # in 2014 al = - 0.1
      #al = -0.1 => bank  al=0.1 => borrow
      
      #Repay                   and          Borrow again
      iTAC[,FcY,c("A")]         <- ( al * iTAC[,ImY,c("A")] ) + ((1+al) * iTAC[,FcY,c("A")])
    }
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  
  
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, bank, then alternate borrow and bank
  #-----------------------------------------------------------------------------
  if(mpoptions$BBscen == "AlternateBorrow"){
    
    #-  HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    #- first year of  banking
    
    if(iYr == "2013"){             #Borrow
      iTAC[,FcY,c("A")]         <- 0.9 * iTAC[,FcY,"A"]
    }
    if(an(iYr) >= 2013){ 
      
      al<- - (an(iYr)%%2-0.5)/5          # in 2014 al = - 0.1
      #al = -0.1 => bank  al=0.1 => borrow
      
      #Repay                   and          Borrow again
      iTAC[,FcY,c("A")]         <- ( al * iTAC[,ImY,c("A")] ) + ((1+al) * iTAC[,FcY,c("A")])
    }
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, bank, then alternate borrow and bank
  #-----------------------------------------------------------------------------
  if(mpoptions$BBscen == "MinVar"){
    
    #           -  Scenario 6 (yield stability): 
    #  (y) = sign(Y(y-1) +   (y-1) TAChcr(y-1) - TAChcr(y)) 
    #  		* min{|(Y(y-1) +   (y-1) TAChcr(y-1))/TAChcr(y) - 1|, 0.1},
    #for all years y? 2013, where the "sign" function is defined as:
    #sign(x) = 1 if x>0;  = 0 if x=0; = -1 if x>0.
    #
    
    
    #-  HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    #- first year of  banking
    
    if(iYr == "2013"){             #Borrow
      iTAC[,FcY,c("A")]         <- 0.9 * iTAC[,FcY,"A"]
    }
    if(an(iYr) >= 2013){ 
      
      al<- - (an(iYr)%%2-0.5)/5          # in 2014 al = - 0.1
      #al = -0.1 => bank  al=0.1 => borrow
      
      #Repay                   and          Borrow again
      iTAC[,FcY,c("A")]         <- ( al * iTAC[,ImY,c("A")] ) + ((1+al) * iTAC[,FcY,c("A")])
    }
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  res<-list(TAC=iTAC[,ac(FcY)],HCRTAC=HCRTAC[,ac(FcY)],SSB=list(HCRSSB=HCRSSB[,ac(FcY)],SSB=SSB[,ac(FcY)]),fSTF=unitSums(stf@harvest[,FcY]))
  return(res)
  
}


#----------------------------------------------------------------------------------------------#
#
#  Calculate TAC with a fix proportion of a harvestable biomass
#
#-----------------------------------------------------------------------------------------------#


HCR.fix.prop <- function(iStocks,iSurvivors,iFishery,iYr,iTAC,iScenario,iFailRuns=NULL,iHistMaxYr,mpPoints,mpoptions){
  
  #iStocks<-stocks[,1:idxyrmin1]
  #iSurvivors<-survivors                                 
  #iFishery<-tmp_fishery
  #iTAC<-TAC
  #iScenario<-scen
  #iFailRuns=NULL 
  #iHistMaxYr<-histMaxYr
  # mpoptions<-mpOptions
  #  
  require(minpack.lm)
  
  lin     <- substr(R.Version()$os,1,3)== "lin"
  stk     <- iStocks
  stk.sur <- iSurvivors
  #===============================================================================
  # Setup control file
  #===============================================================================
  
  DtY         <- ac(range(stk)["maxyear"]) #Data year
  ImY         <- ac(an(DtY)+1) #Intermediate year
  FcY         <- ac(an(DtY)+2) #Forecast year
  CtY         <- ac(an(DtY)+3) #Continuation year
  CtY1        <- ac(an(DtY)+4)
  FuY         <- c(ImY,FcY,CtY,CtY1)#Future years
  pyears      <- an(DtY) - iHistMaxYr #Years away from original historic max year
  
  #- We substract pyears from the starting point, because with every year forward, we move the year ranges one forward too.
  RECS        <- FLQuants("ImY"=iSurvivors[1,],"FcY"=exp(apply(log(rec(stk)[,ac((range(stk)["maxyear"]-12-pyears):(range(stk)["maxyear"]))]),3:6,mean,na.rm=T)),
                          "CtY"=exp(apply(log(rec(stk)[,ac((range(stk)["maxyear"]-12-pyears):(range(stk)["maxyear"]))]),3:6,mean,na.rm=T)))
  
  yrs1        <- list("m.spwn","harvest.spwn","stock.wt")
  yrs3        <- list("mat")
  yrs5        <- list("m")
  
  dsc         <- "NEa mackerel  "
  nam         <- "Mac"
  dms         <- dimnames(stk@m)
  dms$year    <- c(rev(rev(dms$year)[1:3]),ImY,FcY,CtY,CtY1)
  dms$unit    <- dimnames(iFishery@landings.n)$unit
  
  
  f48         <- ac(4:8)
  
  #===============================================================================
  # Setup stock file
  #===============================================================================
  
  stf         <- FLStock(name=nam,desc=dsc,m=FLQuant(NA,dimnames=dms))
  for(i in dms$unit) stf[,,i] <- window(stk,start=an(dms$year)[1],end=rev(an(dms$year))[1])
  units(stf)  <- units(stk)
  # Fill slots that are the same for all fleets
  for(i in c(unlist(yrs1),unlist(yrs3),unlist(yrs5))){
    if(i %in% unlist(yrs1)){ for(j in dms$unit){ slot(stf,i)[,FuY,j] <- slot(stk,i)[,DtY]}}
    if(i %in% unlist(yrs3)){ for(j in dms$unit){ slot(stf,i)[,FuY,j] <- yearMeans(slot(stk,i)[,ac((an(DtY)-2):an(DtY))])}}
    if(i %in% unlist(yrs5)){ for(j in dms$unit){ slot(stf,i)[,FuY,j] <- yearMeans(slot(stk,i)[,ac((an(DtY)-4):an(DtY))])}}
  }
  
  # Fill slots that are unique for the fleets
  for(iFuY in FuY){
    stf@harvest[,iFuY]       <- stk@harvest[,DtY]
    stf@catch.wt[,iFuY]      <- iFishery@landings.wt[,ac(DtY)]
    stf@landings.wt[,iFuY]   <- iFishery@landings.wt[,ac(DtY)]
  }
  
  # Fill slots that have no meaning for NSAS
  stf@discards.n[]        <- 0
  stf@discards[]          <- 0
  stf@discards.wt[]       <- 0
  stf@m[]                 <- 0.15
  #===============================================================================
  # Intermediate year
  #===============================================================================
  
  for(i in dms$unit)        stf@stock.n[,ImY,i]     <- stk.sur
  stf@harvest[,ImY]         <- fleet.harvestFF(stk=stf,iYrr=ImY,TACS=iTAC[,ImY])
  for(i in dms$unit){
    stf@catch.n[,ImY,i]     <- stf@stock.n[,ImY,i]*(1-exp(-unitSums(stf@harvest[,ImY])-stf@m[,ImY,i]))*(stf@harvest[,ImY,i]/(unitSums(stf@harvest[,ImY])+stf@m[,ImY,i]))
    stf@catch[,ImY,i]       <- computeCatch(stf[,ImY,i])
    stf@landings.n[,ImY,i]  <- stf@catch.n[,ImY,i]
    stf@landings[,ImY,i]    <- computeLandings(stf[,ImY,i])
  }
  
  #===============================================================================
  # Calculating TAC
  #===============================================================================
  for(i in dms$unit)        stf[,FcY,i]     <- stf[,ImY,i]
  B2      <- quantSums(stf@stock.n[3:13,FcY,1]*stf@catch.wt[3:13,FcY,1])
  SSB     <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                         exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
  
  HCRTAC <- HCRTAC[,FcY]
  
  minFun <- HCRTAC
  minFun[] <- NA
  for(i in 1: dim(minFun)[6]){
    iter(minFun,i)[] <- min(1, as.numeric(iter(SSB/mpPoints$Btrigger,i)))
  }
  HCRTAC[]                    <- mpPoints$alpha*minFun*B2 
  
  #harvestCatch(stf,FcY)
  #  SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
  HCRSSB                    <- SSB
  # F48                       <- apply(unitSums(stf@harvest[f48,FcY]),2:6,mean)
  
 # print(round(c(apply(F48,1,mean)),5))
  
  
  #-----------------------------------------------------------------------------
  #- Scenario no IAV
  #-----------------------------------------------------------------------------
  if(mpoptions$TACvarlim == F | ImY=="2014" ){
    #- No change to TAC calculation
    iTAC[,FcY]              <- HCRTAC
  }
  
  #-----------------------------------------------------------------------------
  #- Scenario current Long Term Management Plan
  #-----------------------------------------------------------------------------
  if(mpoptions$TACvarlim == T & ImY!="2014"){
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  #-----------------------------------------------------------------------------
  #- Scenario: if HCR TAC results in F outside of 15% of target F, remove TAC contraint, else use it
  #-----------------------------------------------------------------------------
#   if(mpoptions$Fvarlim == T){
#     tmpTAC                    <- iTAC[,FcY,]
#     tmpTAC[]                  <- HCRTAC
#     iTAC[,FcY]                <- HCRTAC
#     #- First calculate what the F would be if the TAC would be constraint by 15%
#     bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")])
#     sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")])
#     if(length(bidx)>0) tmpTAC[,,c("A"),,,bidx]   <- 1.20*iTAC[,ImY,c("A"),,,bidx]
#     if(length(sidx)>0) tmpTAC[,,c("A"),,,sidx]   <- 0.80*iTAC[,ImY,c("A"),,,sidx]
#     stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=tmpTAC)
#     
#     for(i in dms$unit){
#       stf@catch.n[,FcY,i]     <- stf@stock.n[,FcY,i]*(1-exp(-unitSums(stf@harvest[,FcY])-stf@m[,FcY,i]))*(stf@harvest[,FcY,i]/(unitSums(stf@harvest[,FcY])+stf@m[,FcY,i]))
#       stf@catch[,FcY,i]       <- computeCatch(stf[,FcY,i])
#       stf@landings.n[,FcY,i]  <- stf@catch.n[,FcY,i]
#       stf@landings[,FcY,i]    <- computeLandings(stf[,FcY,i])
#     }
#     #- Determine which iters are outside the +/- 15% from target F or inside this range
#     outidx                    <- which(apply(unitSums(stf@harvest[f48,FcY]),2:6,mean,na.rm=T) <  0.9 * F48 |
#                                          apply(unitSums(stf@harvest[f48,FcY]),2:6,mean,na.rm=T) >  1.1 * F48)
#     inidx                     <- which(apply(unitSums(stf@harvest[f48,FcY]),2:6,mean,na.rm=T) >= 0.9 * F48 &
#                                          apply(unitSums(stf@harvest[f48,FcY]),2:6,mean,na.rm=T) <= 1.1 * F48)
#     #- If outside this range, do not apply TAC constraint and use 15% cap on F
#     if(length(outidx)>0){
#       sidx                    <- which((apply(unitSums(stf@harvest[f48,FcY,,,,outidx]),2:6,mean,na.rm=T)) <  0.9 * F48[,,,,,outidx])
#       bidx                    <- which((apply(unitSums(stf@harvest[f48,FcY,,,,outidx]),2:6,mean,na.rm=T)) >  1.1 * F48[,,,,,outidx])
#       
#       if(length(sidx)>0) stf@harvest[,FcY,"A",,,outidx[sidx]]       <- sweep(stf@harvest[,FcY,"A",,,outidx[sidx]],2:6,
#                                                                              (F48[,,,,,outidx[sidx]] * 0.9)  /
#                                                                                apply(unitSums(stf@harvest[f48,FcY,"A",,,outidx[sidx]]),2:6,mean,na.rm=T),"*")
#       if(length(bidx)>0) stf@harvest[,FcY,"A",,,outidx[bidx]]       <- sweep(stf@harvest[,FcY,"A",,,outidx[bidx]],2:6,
#                                                                              (F48[,,,,,outidx[bidx]] * 1.1)  /
#                                                                                apply(unitSums(stf@harvest[f48,FcY,"A",,,outidx[bidx]]),2:6,mean,na.rm=T),"*")
#       
#       iTAC[,FcY,c("A"),,,outidx]                <-  quantSums(stf@catch.wt[,FcY,"A",,,outidx] * stf@stock.n[,FcY,"A",,,outidx] *
#                                                                 (1-exp(-unitSums(stf@harvest[,FcY,,,,outidx])-stf@m[,FcY,"A",,,outidx])) *
#                                                                 (stf@harvest[,FcY,"A",,,outidx]/(unitSums(stf@harvest[,FcY,,,,outidx])+stf@m[,FcY,"A",,,outidx])))
#     }
#     if(length(inidx)>0){
#       bidx                    <- which(HCRTAC[,,c("A"),,,inidx] > 1.20*iTAC[,ImY,c("A"),,,inidx])
#       sidx                    <- which(HCRTAC[,,c("A"),,,inidx] < 0.80*iTAC[,ImY,c("A"),,,inidx])
#       if(length(bidx)>0) iTAC[,FcY,c("A"),,,inidx[bidx]] <- 1.20*iTAC[,ImY,c("A"),,,inidx[bidx]]
#       if(length(sidx)>0) iTAC[,FcY,c("A"),,,inidx[sidx]] <- 0.80*iTAC[,ImY,c("A"),,,inidx[sidx]]
#     }
#     #- If SSB below Blim, then no cap at all
#     stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
#     SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
#                                              exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
#     idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
#     if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
#   }
  
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, bank, other years repay + borrow
  #-----------------------------------------------------------------------------
  if(mpoptions$BBscen == "BB"){
    
    #- First HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    #- Second banking and borrowing
    if(iYr == "2013"){             #Bank
      iTAC[,FcY,c("A")]         <- 0.9 * iTAC[,FcY,"A"]
    }
    if(iYr == "2014"){             #Repay banked part                               Borrow
      iTAC[,FcY,c("A")]         <- (iTAC[,ImY,c("A")] / 0.9 - iTAC[,ImY,c("A")]) + (1.1 * iTAC[,FcY,c("A")])
    }
    if(an(iYr) > 2014){            #Repay borrowed part                             Borrow
      iTAC[,FcY,c("A")]         <- (iTAC[,ImY,c("A")] / 1.1 - iTAC[,ImY,c("A")]) + (1.1 * iTAC[,FcY,c("A")])
    }
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, bank, other years : bank
  #-----------------------------------------------------------------------------
  if(mpoptions$BBscen == "Banking"){
    
    #-  HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    #- first year of  banking
    if(iYr == "2013"){             #Bank
      iTAC[,FcY,c("A")]         <- 0.9 * iTAC[,FcY,"A"]
    }
    if(an(iYr) >= 2013){             #Repay banked part                               Bank
      iTAC[,FcY,c("A")]         <- (0.1 * iTAC[,ImY,c("A")] ) + (0.9 * iTAC[,FcY,c("A")])
    }
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, borrow, other years : borrow
  #-----------------------------------------------------------------------------
  if(mpoptions$BBscen == "Borrowing"){
    
    #-  HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    #- first year of  banking
    if(iYr == "2013"){             #Borrow
      iTAC[,FcY,c("A")]         <- 1.1 * iTAC[,FcY,"A"]
    }
    if(an(iYr) >= 2013){             #Repay                   and          Borrow again
      iTAC[,FcY,c("A")]         <- (- 0.1 * iTAC[,ImY,c("A")] ) + (1.1 * iTAC[,FcY,c("A")])
    }
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, borrow, then alternate bank and borrow
  #-----------------------------------------------------------------------------
  if(mpoptions$BBscen == "AlternateBorrow"){
    
    #-  HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    #- first year of  banking
    
    if(iYr == "2013"){             #Borrow
      iTAC[,FcY,c("A")]         <- 1.1 * iTAC[,FcY,"A"]
    }
    if(an(iYr) >= 2013){ 
      
      al<-(an(iYr)%%2-0.5)/5          # in 2014 al = - 0.1
      #al = -0.1 => bank  al=0.1 => borrow
      
      #Repay                   and          Borrow again
      iTAC[,FcY,c("A")]         <- ( al * iTAC[,ImY,c("A")] ) + ((1+al) * iTAC[,FcY,c("A")])
    }
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  
  
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, bank, then alternate borrow and bank
  #-----------------------------------------------------------------------------
  if(mpoptions$BBscen == "AlternateBorrow"){
    
    #-  HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    #- first year of  banking
    
    if(iYr == "2013"){             #Borrow
      iTAC[,FcY,c("A")]         <- 0.9 * iTAC[,FcY,"A"]
    }
    if(an(iYr) >= 2013){ 
      
      al<- - (an(iYr)%%2-0.5)/5          # in 2014 al = - 0.1
      #al = -0.1 => bank  al=0.1 => borrow
      
      #Repay                   and          Borrow again
      iTAC[,FcY,c("A")]         <- ( al * iTAC[,ImY,c("A")] ) + ((1+al) * iTAC[,FcY,c("A")])
    }
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, bank, then alternate borrow and bank
  #-----------------------------------------------------------------------------
  if(mpoptions$BBscen == "MinVar"){
    
    #           -  Scenario 6 (yield stability): 
    #  (y) = sign(Y(y-1) +   (y-1) TAChcr(y-1) - TAChcr(y)) 
    #  		* min{|(Y(y-1) +   (y-1) TAChcr(y-1))/TAChcr(y) - 1|, 0.1},
    #for all years y? 2013, where the "sign" function is defined as:
    #sign(x) = 1 if x>0;  = 0 if x=0; = -1 if x>0.
    #
    
    
    #-  HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.20*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.80*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.20*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.80*iTAC[,ImY,c("A"),,,sidx]
    
    #- first year of  banking
    
    if(iYr == "2013"){             #Borrow
      iTAC[,FcY,c("A")]         <- 0.9 * iTAC[,FcY,"A"]
    }
    if(an(iYr) >= 2013){ 
      
      al<- - (an(iYr)%%2-0.5)/5          # in 2014 al = - 0.1
      #al = -0.1 => bank  al=0.1 => borrow
      
      #Repay                   and          Borrow again
      iTAC[,FcY,c("A")]         <- ( al * iTAC[,ImY,c("A")] ) + ((1+al) * iTAC[,FcY,c("A")])
    }
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                             exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }

  
  res<-list(TAC=iTAC[,ac(FcY)],HCRTAC=HCRTAC[,ac(FcY)],SSB=list(HCRSSB=HCRSSB[,ac(FcY)],SSB=SSB[,ac(FcY)]),fSTF=unitSums(stf@harvest[,FcY]))
 
 return(res)
  
}

