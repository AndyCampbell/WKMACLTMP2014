
#
#install.packages("devtool")
#install.packages("PBSadmb")
#install.packages("R2admb")
#
#library(devtools)
#install_github("einarhjorleifsson/mac")
#
library(mac)
library(knitr)
library(RColorBrewer)
library(PBSadmb)


compile_admb2(windose=T,compile=F)

compile_admb2 <- function(tpl_name="srest.tpl",windose=F,compile=TRUE) {
  #if(!(compile | windose)) stop("Must compile in non-Windows operating systems")
  #if(compile)
  #{
    file.copy(paste(path.package("mac"),"/extdata/",tpl_name,sep=''),tpl_name)
    tpl <- str_replace(tpl_name,".tpl","")
    compile_admb(tpl,verbose=T)
    clean_admb(tpl)
  #} else {
  #  stop(".exe file not yet available")
    #file.copy(paste(path.package("mac"),"/bin/",tpl,".exe",sep=""),paste(tpl,".exe","")
  #}
}

sam <- read.fit(paste(inPath,"NEA-Mac-WGWIGE-2014-V2/run/sam",sep=""))
rby <- read.rby(sam)
sr <- rby[rby$year < 2013 &rby$year >= 1990 ,c("year","r","ssbest")]
names(sr) <- c("year","r","ssb")
mac_srest_cat("Mackerel",
              y1=min(sr$year),y2=max(sr$year),
              aR=0,aP=12,
              opt_sr_model = 3,opt_pen = 1,
              r = sr$r,ssb=sr$ssb,year=sr$year)



