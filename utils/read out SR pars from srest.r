
rm(list=ls())
load("C:\\Users\\brune001\\Dropbox\\WKMACLTMP (1)\\data\\sr_fit.RData")

nits<-length(table(fit$iter))
nmods<-length(table(fit$model))

SRmod<-expand.grid(model=c(1:nmods),iter=c(1:nits))

for (i in 1:dim(SRmod)[1])
{

mod  <- SRmod$model[i]
iter <- SRmod$iter[i]

SRmod$ap[i]              <-  fit$value[fit$model==mod & fit$iter==iter & fit$name=="ap"]
SRmod$a[i]               <-  fit$value[fit$model==mod & fit$iter==iter & fit$name=="a"]
SRmod$bp[i]              <-  fit$value[fit$model==mod & fit$iter==iter & fit$name=="bp"]
SRmod$b[i]               <-  fit$value[fit$model==mod & fit$iter==iter & fit$name=="b"]
SRmod$sigR[i]            <-  fit$value[fit$model==mod & fit$iter==iter & fit$name=="sigR"]
SRmod$scor[i]            <-  fit$value[fit$model==mod & fit$iter==iter & fit$name=="scor"]
SRmod$g[i]               <-  fit$value[fit$model==mod & fit$iter==iter & fit$name=="g"]
SRmod$neglnL[i]          <-  fit$value[fit$model==mod & fit$iter==iter & fit$name=="neglnL"]

}

write.csv(SRmod,file="W:/IMARES/Data/ICES-WG/WKMACLTMP/Data/SR pars from srest.csv")

