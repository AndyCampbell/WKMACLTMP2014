  #-------------------------------------------------------------------------------
  # 5): compute the deviation from the "true" stock (Mac) for each replicate
  #-------------------------------------------------------------------------------


# get the CVs on the Ns and Fs from the latests assessment
devN<- log(sweep(stock.n(stocks)[,ac(histPeriod)],1:5,stock.n(Mac)[,ac(histPeriod)]  ,"/"))      # this is actually cheating, I should trake the CVs from the output of SAM, and not recalculate them from the MCMC realisiations
devF<- log(sweep(harvest(stocks)[,ac(histPeriod)],1:5,harvest(Mac)[,ac(histPeriod)]  ,"/"))

cvN<-(iterVars(devN))^0.5
cvF<-(iterVars(devF))^0.5


# create objects which have the same dimension as the assessment error objects ie (age, years, projection years,1,1,iterations)
cvN<-propagate(cvN,1000)
dnms<-dimnames(cvN)
dnms$unit<-ac(c(2013,projPeriod))
cvN2<-FLQuant(NA,dimnames=dnms)
for (un in dnms$unit)  cvN2[,,un,,,]<-cvN

cvF<-propagate(cvF,1000)
dnms<-dimnames(cvF)
dnms$unit<-ac(c(2013,projPeriod))
cvF2<-FLQuant(NA,dimnames=dnms)
for (un in dnms$unit)  cvF2[,,un,,,]<-cvF



rhoN<-0.5  # assume an autocorrelation for the moment
rhoF<-0.5

#create an array which repeats the standard FLquant as many times as there are projection years
dnms<-dimnames(stock.n(stocks)[,histPeriod])
dnms$unit<-ac(c(2013,projPeriod))
dnms$iter<-1:1000
En<-FLQuant(NA,dimnames=dnms)
Ef<-En

# generate the autocorrelated errors
En<-           (1-rhoN^2)^0.5*rnorm(length(c(En[]@.Data)),0,1)*cvN2
Ef<-           (1-rhoF^2)^0.5*rnorm(length(c(Ef[]@.Data)),0,1)*cvF2


En[-1,-1,-1] <-   rhoN* En[-dim(En)[1],-dim(En)[2],-dim(En)[3]] + En[-1,-1,-1]
Ef[-1,-1,-1] <-   rhoF* Ef[-dim(Ef)[1],-dim(Ef)[2],-dim(Ef)[3]] + Ef[-1,-1,-1]


