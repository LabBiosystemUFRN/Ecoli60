rm(list = ls())

type="High"
pval=0.01
quant=7
normalize = T
plot = T
rank=100
tail="L"
TsTv = "all"
fix = T
graph = "box"
workdir = "/home/clovis/Dropbox/Ecoli60/data_files/"

enrichFixated<- function(type="High",
                  pval=0.01,
                  quant=7, 
                  normalize = F, 
                  rank=100, 
                  tail="L",
                  TsTv ,
                  fix = T,
                  graph = "box",
                  workdir = "/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
  if(tail%in%c("L","H")){
    tail<- (tail == "L")
  }else{
    cat('Use values "L" (Low) or "H" (High) for tail')
    return(0)
  }
  if(!TsTv%in%c("all","Ts","Tv")){
    cat('Use values "all","Ts" (transitions), or "Tv" (transversions) for TsTv')
    return(0)
  }
  vTsTv<-TsTv
  
  setwd(workdir)
  deltaW<-read.csv("cDeltaW.csv",
                   header = T,
                   stringsAsFactors = F)
  codonUsage<-read.csv("dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  
  
  mutAllel<-read.csv(paste0("a",type,"MutSynFreqDW.csv"),
                     header = T,
                     stringsAsFactors = F)
  mutAllelAll<-read.csv(paste0("a",type,"MutAllDW.csv"),
                     header = T,
                     stringsAsFactors = F)
  
  #só syn
  mutAllel<-mutAllelAll[mutAllelAll$Annotation == "synonymous",]
  
  #remove stop codons
  mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"),]
  #filtra Ts ou Tv
  if(TsTv != "all"){
    mutAllel<-subset(mutAllel, TsTv == vTsTv)
  }
  if(fix){
    mutAllel<-mutAllel[!is.na(mutAllel$fixed),]
  }else{
    mutAllel<-mutAllel[is.na(mutAllel$fixed),]
  }
  library(reshape2)

  quantil<-quantile(deltaW$LeGac,c(seq(from = 0,to = 1,by = 1/quant)),type = 1)
  quantil[1]<-quantil[1]-0.1
  
  for(pop in unique(mutAllel$pop)){
    tmp<-mutAllel[mutAllel$pop==pop,c("Gene","dw")]
    
    for(i in 1:(length(quantil)-1)){
      if(i == 1){
        tmp$faixa[tmp$dw == quantil[1]]<-1
      }
      tmp$faixa[tmp$dw> quantil[i] & tmp$dw<=quantil[i+1]]<-i
    }
    
    tmp$dw<-NULL
    tmp<-as.data.frame(table(tmp))
    tmp2<-dcast(tmp,Gene~faixa)
    #rownames(tmp2)<-tmp2$Gene
    # tmp2$sum<-apply(X = tmp2[,2:ncol(tmp2)],
    #                 MARGIN = 1,sum)
    # tmp3<-tmp2[tmp2$sum>=2,]
    # tmp3$sum<-NULL
    library(ggplot2)
    library("pheatmap")
    pheatmap(tmp2[,2:ncol(tmp2)],cluster_cols = F,main=pop)
  }
  p<-ggplot()+theme_bw()+
    labs(title = ifelse(fix,"Fixed Mutations","Never Fixed Mutations"))+
    xlab("Transit time (x1000)")+
    ylab(expression(paste(Delta,"w")))
  i=2
  for(i in 2:(ncol(tmp3))){
    tmp4<-data.frame(y=tmp3[,i])
    p<-p+geom_boxplot(data = tmp4, aes(x=i, y= y))
  }
  p
  tmp2$Gene<-NULL
  sum(tmp2$sum)
  length(unique(tmp$Gene))
  nrow(tmp[tmp$Freq>2,])
  
  
  
  if(graph == "box"){
    library(ggplot2)
    interval=1/quant
    quantil<-quantile(mutAllel$transit,seq(0,1,interval))
    for(i in 1:(length(quantil)-1)){
      if(i == 1){
        mutAllel$nfaixa[mutAllel$transit == quantil[1]]<-paste(quantil[i]/1000,"to",quantil[i+1]/1000)
        mutAllel$faixa[mutAllel$transit == quantil[1]]<-1#paste(quantil[i]/1000,"to",quantil[i+1]/1000)
      }
      mutAllel$nfaixa[mutAllel$transit> quantil[i] & mutAllel$transit<=quantil[i+1]]<-paste(quantil[i]/1000,"to",quantil[i+1]/1000)
      mutAllel$faixa[mutAllel$transit> quantil[i] & mutAllel$transit<=quantil[i+1]]<-i#paste(quantil[i]/1000,"to",quantil[i+1]/1000)
    }
    
    
      
    p<-ggplot()+theme_bw()+
      labs(title = ifelse(fix,"Fixed Mutations","Never Fixed Mutations"))+
      xlab("Transit time (x1000)")+
      ylab(expression(paste(Delta,"w")))+
      geom_boxplot(data = mutAllel, aes(x= factor(nfaixa), y= dw))
    return(p)
  }
}


p<-enrichFixated(type="High",
               pval=0.001,
               quant=8, 
               normalize = T, 
               rank=100,
               graph = "box",
               fix = F,
               TsTv = "all",
               workdir = "/home/clovis/Dropbox/Ecoli60/data_files/")
  
p
