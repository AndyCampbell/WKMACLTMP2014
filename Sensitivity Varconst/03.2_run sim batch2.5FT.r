
rm(list=ls())
library(FLCore)
#library(FLFleet)
#library(FLAssess)
#library(FLSAM)
library(MASS)
#library(msm)


for (opt in c(6,7,8,11,16,21))
{
# permanent changes?
perm<-T
if (perm)cat("!!! scenario with permanent changes")
# numer of years for the projections
prlength    <- 40
# reduce the number of iterations to
short<-T
nits2<-200
TACvarlim<-F
Fvarlim<-T
scen          <- c("LTMP2.5mt")              # 
BBscen        <- "noBB"
                                                     # "Banking"          : always bank 
                                                     # "Borrowing"        : always borrow
                                                     # "AlternateBank"    : bank first and then alternate
                                                     # "AlternateBorrow"  : borrow first and then alternate
                                                     # "MinVar"           : use BB to minise TAC variability
                                                     
LastRecOp     <- "geom"                  # option to replace the last estimated recruitment : "SAM", "geom", "RCT3"
                                                     # "SAM"  = don't overwrite the SAM estimmate
                                                     # "geom" = replace by geomean 1990/(TaY-1)
                                                     # "RCT3" = replace by RCT3 output


# run the MSE


codePath      <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/R code/Sensitivity Varconst/"
if(substr(R.Version()$os,1,3)== "lin")
{
  codePath    <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",codePath)
  }
source(paste(codePath,'03.1_run MSE for sim batch.r',sep=""))
}

