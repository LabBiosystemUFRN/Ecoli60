rm(list = ls())



timeline<-function(file="bSumHigh.csv",
                   quant = 7, 
                   normalize = 0,
                   type = "codon"){
  setwd("/home/clovis/Doutorado/Projetos/Ecoli60/data_files/")
  somat<-read.csv(file = file,
                  stringsAsFactors = F)
  somat<-na.exclude(somat)
  
  codonUsage<-read.csv("dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  deltaW<-read.csv("cDeltaW.csv",
                   header = T,
                   stringsAsFactors = F)
  
  #quant=7
  #normalize=0
  somat<-merge(somat,codonUsage[,c(1,3)],by.x = "codon",by.y=c("ref"))
  p<-somat
  if(normalize){
    somat[,3:124]<-somat[,3:124]/somat$factor
  }
  
  
  quantil<-quantile(deltaW$Yoon,c(seq(from = 0,to = 1,by = 1/quant)),type = 1)
  vdw<-data.frame(min=round(quantil[seq(1,quant)],2),max=round(quantil[seq(2,quant+1)],2))
  quantil[1]<-quantil[1]-0.1
  
  for(i in 1:(length(quantil)-1)){
    somat$faixa[somat$dW> quantil[i] & somat$dW<=quantil[i+1]]<-i
  }
  
  library(ggplot2)
  p<-ggplot()+
    theme_bw()+
    xlab("Generations (x10,000)")
  if(normalize){
    p<-p+ylab("Normalized Frequency")
  }else{
    p<-p+ylab("Frequency")
  }
  i=1
  colors<-rainbow(n=quant)
  max<-0
  
  for(i in 1:quant){
    q1<-somat[somat$faixa == i,]
    qt1<-data.frame(x=seq(0,60.5,0.5),
                    y=apply(X = q1[,3:124],MARGIN = 2, FUN = sum))
    max<-max(max,max(qt1$y))
  }
  min<-max/3
  #colors<-c("red","blue","green","orange")
  legenda<-data.frame(x=rep(0,quant),
                      y=seq(max,min,length.out = quant),
                      min=vdw$min,
                      max=vdw$max,
                      color=colors)
  for(i in 1:quant){
    q1<-somat[somat$faixa == i,]
    qt1<-data.frame(x=seq(0,60.5,0.5),
                    y=apply(X = q1[,3:124],MARGIN = 2, FUN = sum))
    p<-p+geom_line(data = qt1, aes(x,y),col=colors[i])
    p<-p+geom_point(data = legenda[i,],aes(x,y),col=colors[i])
    p<-p+geom_text(data = legenda[i,],aes(x+5,y,label=paste(min,"to",max)),cex=3)
  }
  
  
  plot(p)
}

timeline(file="bSumHigh.csv",normalize = T)
