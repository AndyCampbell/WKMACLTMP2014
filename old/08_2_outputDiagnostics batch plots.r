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

TACvarlim     <- T                 
Fvarlim       <- T                 
BBscen        <- "noBB"   
LastRecOp     <- "geom"          


diags<-read.csv(file=paste(outPath,"/tables_diags_screening","_TACvarlim",TACvarlim,"_Fvarlim",Fvarlim,"_",BBscen,"_LastRec.csv",sep=""))



#write.csv(t(diags),file=paste(outPath,"/tables_diags",paste(scen,opt,"_TACvarlim",TACvarlim,"_Fvarlim",Fvarlim,"_",BBscen,"_LastRec",LastRecOp,sep=""),".csv",sep=""),row.names=T)

