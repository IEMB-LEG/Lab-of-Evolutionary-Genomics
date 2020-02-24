setwd("")
library(zoo)
file<- read.table("coverage.matrix",header=F,sep=' ')
x<-rollapply(file,1000,mean,by=500) #500=step,1000=window 
write.table(x,"coverage.1kb.window")
