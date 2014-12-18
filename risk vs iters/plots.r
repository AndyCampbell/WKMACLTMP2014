

outPath       <- "W:/IMARES/Data/ICES-WG/WKMACLTMP/Results/risk vs iters/"
res<-read.csv(file=paste(outPath,"res50.csv",sep=""))


vars<-c("Risk1MidT","Risk1LongT","Risk2MidT","Risk2LongT","Risk3MidT","Risk3LongT")
par(mfrow=c(3,2))
for (vv in vars)  boxplot(res[,vv]~res$sample,main=vv)
