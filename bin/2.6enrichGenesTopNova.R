rm(list = ls())

type="HighLow"
pval=0.01
quant=5
normalize = T
plot = T
rank=0.125
tail="L"
TsTv = "all"
fix = "all"
graph = "box"
workdir = "/home/clovis/Dropbox/Ecoli60/data_files/"


# enrichTopGenes<- function(type="High",
#                   pval=0.01,
#                   quant=7, 
#                   normalize = F, 
#                   rank=100, 
#                   tail="L",
#                   TsTv ,
#                   fix = T,
#                   graph = "box",
#                   workdir = "/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
  if(rank<0){
    cat("Rank must be an integer or a number between 0 and 1.")
    return(0)
  }
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
  deltaW<-read.csv("AuxFiles/cDeltaW.csv",
                   header = T,
                   stringsAsFactors = F)
  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  
  
  # mutAllel<-read.csv(paste0("a",type,"MutSynFreqDW.csv"),
  #                    header = T,
  #                    stringsAsFactors = F)
  mutAllelAll<-read.csv(paste0("base/a",type,"MutAllDW.csv"),
                     header = T,
                     stringsAsFactors = F)
  
  #só syn
  mutAllel<-mutAllelAll[mutAllelAll$Annotation == "synonymous",]
  
  hiGenes<-read.csv("AuxFiles/highlyEGLeGac",
                    header = F,
                    stringsAsFactors = F)
  colnames(hiGenes)<-"gene"
  #hiGenes<-data.frame(gene=hiGenes[-(grep(pattern = "#",hiGenes$V1)),])
  
  #remove stop codons
  mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"),]
  mutAllelAll<-mutAllelAll[!(mutAllelAll$ref == "TAA"| mutAllelAll$ref == "TGA"),]
  #filtra Ts ou Tv
  if(TsTv != "all"){
    mutAllel<-subset(mutAllel, TsTv == vTsTv)
    mutAllelAll<-subset(mutAllelAll, TsTv == vTsTv)
  }
  
  if(fix !="all"){
    if(fix){
      mutAllel<-mutAllel[!is.na(mutAllel$fixed),]
      mutAllelAll<-mutAllelAll[!is.na(mutAllelAll$fixed),]
    }else{
      mutAllel<-mutAllel[is.na(mutAllel$fixed),]
      mutAllelAll<-mutAllelAll[is.na(mutAllelAll$fixed),]
    }
  }

  #acrescenta fator
  mutAllel<-merge(mutAllel,codonUsage[,c(1,3)],by=c("ref"))
  colnames(mutAllel)[136]<-"factor"
  
  #normaliza valores
  tmp<-mutAllel[,7:128]/mutAllel$factor
  tmp[is.na(tmp)]<- 0
  mutAllelNorm<-cbind(mutAllel[,1:6],tmp,mutAllel[129:136])
  
  rm(tmp)
  
  #calcula a frequencia maxima atingida em cada mutação
  mutAllelNorm$maxFreq<-apply(X = mutAllelNorm[,7:128],
                         MARGIN = 1,
                         max) 
  # acha max e min
  # minNorm<-min(mutAllelNorm$maxFreq)
  # maxNorm<-max(mutAllelNorm$maxFreq)
  # mutHigh<-mutAllelNorm[mutAllelNorm$Gene%in%hiGenes$gene,]
  # mutOthers<-mutAllelNorm[!mutAllelNorm$Gene%in%hiGenes$gene,]
  # 
  # min(mutHigh$maxFreq)
  # min(mutOthers$maxFreq)
  # 
  # median<-median(mutAllelNorm$maxFreq)
  
  # quantil<-quantile(deltaW$LeGac,c(seq(from = 0,to = 1,by = 1/8)),type = 1)
  # quantil[1]<-quantil[1]-0.1
  # 
  mutAllel<-mutAllelNorm
  rm(mutAllelNorm,mutAllelAll)
  mutAllel<-mutAllel[order(mutAllel$Position),]

  #teste estatistico de ocorrencia de higenes
  #total de high genes e outros genes
  if(rank>=0 & rank <=1){
    rankT<- round(as.numeric(nrow(mutAllel))*as.numeric(rank),0)
  }else{
    rankT<- round(rank,0)
  }
  
  top100<-mutAllel
  top<-top100$maxFreq[order(mutAllel$maxFreq, decreasing = F)]
  if(rankT > length(top))
    rankT <- length(top)
  if(top[rankT] == 0){
    top100<-top100[top100$maxFreq<top[rankT],]
  }else{
    top100<-top100[top100$maxFreq<=top[rankT],]
  }
  
  genesWhite<-length(unique(mutAllel$Gene[mutAllel$Gene%in%hiGenes$gene]))
  genesBlack<-length(unique(mutAllel$Gene[!mutAllel$Gene%in%hiGenes$gene]))
  genesFreq<-length(unique(top100$Gene[top100$Gene%in%hiGenes$gene]))
  genesNFreq<-length(unique(top100$Gene[!top100$Gene%in%hiGenes$gene]))
  genesDrawn<-length(unique(top100$Gene))
  
  phyper(genesFreq,genesWhite,genesBlack,genesDrawn,lower.tail = F)
  dbinom(genesFreq,genesDrawn,genesWhite/genesBlack)
  chisq.test(c(genesFreq,genesNFreq)
             ,c(genesWhite/genesBlack,
             1-(genesWhite/genesBlack)))
  
  library(DescTools)
  GTest(c(genesFreq,genesNFreq)
             ,c(genesWhite/genesBlack,
                1-(genesWhite/genesBlack)))
  
  
  quantil<-quantile(deltaW$dw,c(seq(from = 0,to = 1,by = 1/quant)),type = 1)
  vdw<-data.frame(min=round(quantil[seq(1,quant)],2),max=round(quantil[seq(2,quant+1)],2))
  quantil[1]<-quantil[1]-0.01
  
  for(i in 1:(length(quantil)-1)){
    if(i == 1){
      mutAllel$faixa[mutAllel$dw == quantil[1]]<-1
    }
    mutAllel$faixa[mutAllel$dw> quantil[i] & mutAllel$dw<=quantil[i+1]]<-i
  }
  
  #ordem<-unique(mutAllel$faixa[order(mutAllel$faixa)]) 
  tmp<-mutAllel[mutAllel$Gene%in%hiGenes$gene,c(4,137,138)]
  mutCount<-length(unique(tmp$Gene))
  
  tmp[70,1]<-"nada"
  tmp[70,2]<-0.0
  tmp[70,3]<-7
  tmp[71,1]<-c("nada",0,7)
  tmp$faixa<-factor(tmp$faixa)
  p<-ggplot()+theme_bw()+
    labs(title = paste("Top expressed genes:",
                       mutCount,
                       "of",
                       nrow(hiGenes)))+
    xlab(expression(paste(Delta,"w levels")))+
    ylab("Max Frequency")+
    #ylim(-1.9,1.9)+
    geom_boxplot(data = tmp, aes(x= faixa, y= maxFreq,fill=faixa))+
    scale_color_manual(expression(paste0(Delta,"w")) )
  
  p

  tmp<-mutAllel[!mutAllel$Gene%in%hiGenes$gene,c(4,137,138)]
  mutCount<-length(unique(tmp$Gene))
  
  tmp$faixa<-factor(tmp$faixa)
  p<-ggplot()+theme_bw()+
    labs(title = paste("Genes:",
                       mutCount,
                       "of",
                       length(unique(mutAllel$Gene))))+
    xlab(expression(paste(Delta,"w levels")))+
    ylab("Max Frequency")+
    #ylim(-1.9,1.9)+
    geom_boxplot(data = tmp, aes(x= faixa, y= maxFreq,fill=faixa))+
    scale_color_manual(expression(paste0(Delta,"w")) )
  
  
  p

  p<-ggplot()+theme_bw()+
    labs(title = paste("Genes:",
                       mutCount,
                       "of",
                       length(unique(mutAllel$Gene))))+
    xlab(expression(paste(Delta,"w levels")))+
    ylab("Max Frequency")+
    ylim(0,0.2)+
    #ylim(-1.9,1.9)+
    geom_boxplot(data = tmp, aes(x= faixa, y= maxFreq,fill=faixa),
                 outlier.shape = NA)+
    scale_color_manual(expression(paste0(Delta,"w")) )
  
  
  p
  
    
  codons<-unique(mutAllel[,c("ref","mut")])
  
  
  
  
  top100Base<-mutAllel[,c(1:5,135,137)]#,2,127,128,3)]

  quantil<-quantile(deltaW$LeGac,c(seq(from = 0,to = 1,by = 1/quant)),type = 1)
  vdw<-data.frame(min=round(quantil[seq(1,quant)],2),max=round(quantil[seq(2,quant+1)],2))
  quantil[1]<-quantil[1]-0.1
  
  for(i in 1:(length(quantil)-1)){
    if(i == 1){
      top100Base$faixa[top100Base$dw == quantil[1]]<-1
    }
    top100Base$faixa[top100Base$dw> quantil[i] & top100Base$dw<=quantil[i+1]]<-i
  }
  #garante a mesma ordem me mutAllel e top100Base
  top100Base<-top100Base[order(top100Base$Position),c("Position","Gene","Allele","ref","mut","dw","maxFreq","faixa")]

  if(exists("result")){rm(result)}
  result<-data.frame(faixa=numeric(),
                     freq=numeric(),
                     white=numeric(),
                     black=numeric(),
                     drawn=numeric(),
                     hypDrawn=numeric(),
                     hypWhite=numeric())
  
  i=1
  for(i in 1:(length(quantil)-1)){
    top100<-subset(top100Base,faixa == i)
    result[i,1]<-i
    if(nrow(top100)==0){
      next()
    }
    
    result$white[i]<-nrow(top100[top100$Gene%in%hiGenes$gene,])
    result$black[i]<-nrow(top100[!top100$Gene%in%hiGenes$gene,])

    if(rank>=0 & rank <=1){
      rankT<- round(as.numeric(nrow(top100))*as.numeric(rank),0)
    }else{
      rankT<- round(rank,0)
    }
    
    top<-top100$maxFreq[order(top100$maxFreq, decreasing = F)]
    if(rankT > length(top))
      rankT <- length(top)
    if(top[rankT] == 0){
      top100<-top100[top100$maxFreq<top[rankT],]
    }else{
      top100<-top100[top100$maxFreq<=top[rankT],]
    }
    result$drawn[i]<-nrow(top100)
    result$freq[i]<-nrow(top100[top100$Gene%in%hiGenes$gene,])
    result$hypDrawn[i]<-phyper(result$freq[i],
                        result$white[i],
                        result$black[i],
                        result$drawn[i],
                        lower.tail = tail)

  }
#   
#   
#   
#     result1<- data.frame(all=numeric(), top=numeric())
#     cat("Fisher Low frequency:\n")
#     i=1
#     library(DescTools)
#     for(i in 1:(length(quantil)-1)){
#       
#     result1[nrow(result1)+1,1]<- length(mutOthers[mutOthers$dw> quantil[i] & 
#                                            mutOthers$dw<=quantil[i+1] & 
#                                            mutOthers$maxFreq<median,"maxFreq"])
#     result1[nrow(result1),2]<- length(mutHigh[mutHigh$dw> quantil[i] & 
#                                               mutHigh$dw<=quantil[i+1] & 
#                                               mutHigh$maxFreq<median,"maxFreq"])
#   }
#   fisher.test(result1,
#               alternative="t")
#   GTest(result1)
#   result2<- data.frame(all=numeric(), top=numeric())
#   cat("Fisher High frequency:\n")
#   i=1
#   for(i in 1:(length(quantil)-1)){
#     result2[nrow(result2)+1,1]<- length(mutOthers[mutOthers$dw> quantil[i] & 
#                                                   mutOthers$dw<=quantil[i+1] & 
#                                                   mutOthers$maxFreq>=median,"maxFreq"])
#     result2[nrow(result2),2]<- length(mutHigh[mutHigh$dw> quantil[i] & 
#                                               mutHigh$dw<=quantil[i+1] & 
#                                               mutHigh$maxFreq>=median,"maxFreq"])
#   }
#   fisher.test(result2,
#               alternative="t")
#   GTest(result2)
#   sum(result1$all)+sum(result2$all)
#   sum(result1$top)+sum(result2$top)
#   hist(mutHigh$dw[mutHigh$maxFreq < median(mutOthers$maxFreq)])
#   hist(mutOthers$dw)
#   
#   
#   non<-mutFiltered[mutFiltered$Annotation != "synonymous",]
#   syn<-mutFiltered[mutFiltered$Annotation == "synonymous",]
#   # hist(non$transit)
#   # hist(syn$transit)
#   # hist(non$fixed,breaks = 6)
#   # hist(syn$fixed,breaks = 6)
#   
#   length(unique(mutFiltered$Gene))
#   
#   quantil<-quantile(deltaW$LeGac,c(seq(from = 0,to = 1,by = 1/quant)),type = 1)
#   vdw<-data.frame(min=round(quantil[seq(1,quant)],2),max=round(quantil[seq(2,quant+1)],2))
#   quantil[1]<-quantil[1]-0.1
#   
#   for(i in 1:(length(quantil)-1)){
#     if(i == 1){
#       syn$faixa[syn$dw == quantil[1]]<-1
#       mutAllel$faixa[mutAllel$dw == quantil[1]]<-1
#     }
#     syn$faixa[syn$dw> quantil[i] & syn$dw<=quantil[i+1]]<-i
#     mutAllel$faixa[mutAllel$dw> quantil[i] & mutAllel$dw<=quantil[i+1]]<-i
#   }
#   
#   countsTot<-as.data.frame(table(mutAllel[,c("faixa")]))
#   colnames(countsTot)<- c("faixa","white")
#   countsTot$black<-nrow(mutAllel)-countsTot$white
#   sum(countsTot$white)
#   
#   counts100<-as.data.frame(table(syn[,c("faixa")]),stringsAsFactors = F)
#   #counts100<-counts100[counts100$Freq!=0,]
#   colnames(counts100)<-c("faixa","Freq")
#   #completa a tabela de counts com zeros
#   faixas<-seq(1:quant)
#   for (padrao in (faixas[!faixas%in%counts100$faixa])) {
#     counts100[nrow(counts100)+1,]<-c(padrao,0)
#   }
#   counts100<-merge(counts100,countsTot,by=c("faixa"))
#   counts100$hyp<-phyper(counts100$Freq,
#                         counts100$white,
#                         counts100$black,
#                         nrow(mutAllel),
#                         lower.tail = T)
#   counts100<-counts100[, c("faixa","hyp")]
#   
#   
#   
#   library(reshape2)
# 
#   quantil<-quantile(deltaW$LeGac,c(seq(from = 0,to = 1,by = 1/quant)),type = 1)
#   quantil[1]<-quantil[1]-0.1
#   
#   for(pop in unique(mutAllel$pop)){
#     tmp<-mutAllel[mutAllel$pop==pop,c("Gene","dw")]
#     
#     for(i in 1:(length(quantil)-1)){
#       if(i == 1){
#         tmp$faixa[tmp$dw == quantil[1]]<-1
#       }
#       tmp$faixa[tmp$dw> quantil[i] & tmp$dw<=quantil[i+1]]<-i
#     }
#     
#     tmp$dw<-NULL
#     tmp<-as.data.frame(table(tmp))
#     tmp2<-dcast(tmp,Gene~faixa)
#     #rownames(tmp2)<-tmp2$Gene
#     # tmp2$sum<-apply(X = tmp2[,2:ncol(tmp2)],
#     #                 MARGIN = 1,sum)
#     # tmp3<-tmp2[tmp2$sum>=2,]
#     # tmp3$sum<-NULL
#     library(ggplot2)
#     library("pheatmap")
#     pheatmap(tmp2[,2:ncol(tmp2)],cluster_cols = F,main=pop)
#   }
#   p<-ggplot()+theme_bw()+
#     labs(title = ifelse(fix,"Fixed Mutations","Never Fixed Mutations"))+
#     xlab("Transit time (x1000)")+
#     ylab(expression(paste(Delta,"w")))
#   i=2
#   for(i in 2:(ncol(tmp3))){
#     tmp4<-data.frame(y=tmp3[,i])
#     p<-p+geom_boxplot(data = tmp4, aes(x=i, y= y))
#   }
#   p
#   tmp2$Gene<-NULL
#   sum(tmp2$sum)
#   length(unique(tmp$Gene))
#   nrow(tmp[tmp$Freq>2,])
#   
#   
#   
#   if(graph == "box"){
#     library(ggplot2)
#     interval=1/quant
#     quantil<-quantile(mutAllel$transit,seq(0,1,interval))
#     for(i in 1:(length(quantil)-1)){
#       if(i == 1){
#         mutAllel$nfaixa[mutAllel$transit == quantil[1]]<-paste(quantil[i]/1000,"to",quantil[i+1]/1000)
#         mutAllel$faixa[mutAllel$transit == quantil[1]]<-1#paste(quantil[i]/1000,"to",quantil[i+1]/1000)
#       }
#       mutAllel$nfaixa[mutAllel$transit> quantil[i] & mutAllel$transit<=quantil[i+1]]<-paste(quantil[i]/1000,"to",quantil[i+1]/1000)
#       mutAllel$faixa[mutAllel$transit> quantil[i] & mutAllel$transit<=quantil[i+1]]<-i#paste(quantil[i]/1000,"to",quantil[i+1]/1000)
#     }
#     
#     
#       
#     p<-ggplot()+theme_bw()+
#       labs(title = ifelse(fix,"Fixed Mutations","Never Fixed Mutations"))+
#       xlab("Transit time (x1000)")+
#       ylab(expression(paste(Delta,"w")))+
#       geom_boxplot(data = mutAllel, aes(x= factor(nfaixa), y= dw))
#     return(p)
#   }
# }


p<-enrichTopGenes(type="High",
               pval=0.001,
               quant=8, 
               normalize = T, 
               rank=100,
               graph = "box",
               fix = T,
               TsTv = "all",
               workdir = "/home/clovis/Dropbox/Ecoli60/data_files/")
  
p
