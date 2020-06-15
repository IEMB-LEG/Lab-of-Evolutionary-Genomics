file=read.table("D://Input.txt",header=T)

randomizer=function(file) {
  x=NULL
  for (i in 1:length(file$Lines)){
    x=c(x,rep(file$Lines[i],file$Replicates[i]))}
  z=matrix(sample(x,length(x),replace=FALSE),8,length(x)/8)
  write.csv(z,"D://file2.txt")
}

randomizer(file)