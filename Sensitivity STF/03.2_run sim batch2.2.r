
rm(list=ls())
library(FLCore)
#library(FLFleet)
#library(FLAssess)
#library(FLSAM)
library(MASS)
#library(msm)


for (LastRecOp in c("geom","RCT3","SAM"))
{

# permanent changes?
perm<-T
if (perm)cat("!!! scenario with permanent changes")
# numer of years for the projections
prlength    <- 40
# reduce the number of iterations to
short<-T
nits2<-100
opt   <-9
BBscen<-"noBB"
scen          <- c("LTMP2.2mt")              # 
TACvarlim     <- T                      # whether or not to apply the 20% limit on TAC change
Fvarlim       <- T                      # whether or not to apply the 10% limit on Fbar change
                                                     # "Banking"          : always bank 
                                                     # "Borrowing"        : always borrow
                                                     # "AlternateBank"    : bank first and then alternate
                                                     # "AlternateBorrow"  : borrow first and then alternate
                                                     # "MinVar"           : use BB to minise TAC variability
                                                     



# run the MSE


codePath      <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/R code/Sensitivity STF/"
if(substr(R.Version()$os,1,3)== "lin")
{
  codePath    <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",codePath)
  }
source(paste(codePath,'03.1_run MSE for sim batch.r',sep=""))
}

