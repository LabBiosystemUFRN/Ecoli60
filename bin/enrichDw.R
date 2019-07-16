enrich<- function(type="High",
                  pval=0.01,
                  quant=7, 
                  normalize = F, 
                  plot = T, 
                  rank=100, 
                  tail="L"){
  if(tail%in%c("L","H")){
    tail<- (tail == "L")
  }else{
    cat('Use values "L" (Low) or "H" (High) for tail')
    return(0)
  }
  #rm(list = ls())
  setwd("/home/clovis/Doutorado/Projetos/Ecoli60/data_files/")
  deltaW<-read.csv("cDeltaW.csv",
                   header = T,
                   stringsAsFactors = F)
  codonUsage<-read.csv("dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  
  
  # class(deltaW$Yoon)
  # hist(deltaW$Yoon,breaks = 31)
  # hist(deltaW$SharpYoon,breaks = 31)
  
  
  
  mutAllel<-read.csv(paste0("a",type,"MutAllDW.csv"),
                     header = T,
                     stringsAsFactors = F)
  #remove stop codons
  mutAllel<-mutAllel[!(mutAllel$codon == "TAA"| mutAllel$codon == "TGA"),]
  #ordena por posição
  mutAllel<-mutAllel[order(mutAllel$Position),]
  
  codons<-unique(mutAllel[,c("codon","mut")])
  
  top100Base<-mutAllel[,c(1,2,127,128,3)]
  colnames(top100Base)<-c("Position","Gene","ref","mut","Allele")
  top100Base<-merge(top100Base,deltaW[,c(1,2,4)],by=c("ref","mut"))
  top100Base<-merge(top100Base,codonUsage[,c(1,3)],by=c("ref"))
  colnames(top100Base)<-c("ref","mut","Position","Gene","Allele","dw", "factor")
  #quant=7
  quantil<-quantile(deltaW$Yoon,c(seq(from = 0,to = 1,by = 1/quant)),type = 1)
  vdw<-data.frame(min=round(quantil[seq(1,quant)],2),max=round(quantil[seq(2,quant+1)],2))
  quantil[1]<-quantil[1]-0.1
  #quantil<-c(-2,-0.5,0,0.5,2)
  #contagem
  # for (i in 1:(length(quantil)-1)) {
  #   cat("Qualtil ",i,": ",sum(deltaW$Yoon>quantil[i] & deltaW$Yoon<=quantil[i+1]),"\n")  
  # }
  
  for(i in 1:(length(quantil)-1)){
    if(i == 1){
      top100Base$faixa[top100Base$dw == quantil[1]]<-1
    }
    top100Base$faixa[top100Base$dw> quantil[i] & top100Base$dw<=quantil[i+1]]<-i
  }
  # top100Base$faixa[top100Base$dw<= -.5]<-1
  # top100Base$faixa[top100Base$dw> -.5 & top100Base$dw<=0]<-2
  # top100Base$faixa[top100Base$dw> 0 & top100Base$dw<=0.5]<-3
  # top100Base$faixa[top100Base$dw> .5 ]<-4
  #garante a mesma ordem me mutAllel e top100Base
  top100Base<-top100Base[order(top100Base$Position),c("Position","Gene","Allele","ref","mut","dw","factor","faixa")]
  countsTot<-as.data.frame(table(top100Base[,c("faixa")]))
  countsTot<-countsTot[countsTot$Freq!=0,]
  colnames(countsTot)<- c("faixa","white")
  countsTot$black<-nrow(top100Base)-countsTot$white
  
  
  if(exists("result")){rm(result)}
  col="X33500"
  faixas<-seq(1:quant)
  for( col in colnames(mutAllel[5:126])){
    top100<-cbind(top100Base,mutAllel[,col])
    colnames(top100)<-c("Position","Gene","Allele","ref","mut","dw","factor","faixa","count")
    
    if(normalize){
    #normalize observations
    top100$cNorm<-apply(X = top100[,c("factor","count")],
                        MARGIN = 1,
                        FUN = function(x){
                          return(x[2]/x[1])
                        })
  }else{
    top100$cNorm<-top100$count
  }
    #colnames(top100)<-c("Position","Gene","Allele","ref","mut","dw","factor","faixa","count")
    top<-top100$cNorm[order(top100$cNorm, decreasing = T)]
    if(top[rank] == 0){
      top100<-top100[top100$cNorm>top[rank],]
    }else{
      top100<-top100[top100$cNorm>=top[rank],]
    }
    counts100<-as.data.frame(table(top100[,c("faixa")]),stringsAsFactors = F)
    #counts100<-counts100[counts100$Freq!=0,]
    colnames(counts100)<-c("faixa","Freq")
    for (padrao in (faixas[!faixas%in%counts100$faixa])) {
      counts100[nrow(counts100)+1,]<-c(padrao,0)
    }
    if(sum(counts100$Freq)-nrow(top100) !=0){
      cat(col," fail!")
      next
    }
    counts100<-merge(counts100,countsTot,by=c("faixa"))
    counts100$hyp<-phyper(counts100$Freq,
                          counts100$white,
                          counts100$black,
                          nrow(top100),
                          lower.tail = tail)
    #    counts100<-counts100[counts100$hyp<=pval, c("faixa","hyp")]
    # print(counts100)
    # Sys.sleep(2)
    counts100<-counts100[, c("faixa","hyp")]
    colnames(counts100)<-c("faixa",col)
    if(exists("result")){
      result<-merge(result,counts100,by=c("faixa"), all = T)
    }else{
      result<-counts100
    }
  }
  resCorrigido<-result[0,2:123]
  linha=1
  for( linha in faixas) {
    resCorrigido[nrow(resCorrigido)+1,]<-c(p.adjust(result[linha,2:123]))
  }
  resCorrigido[resCorrigido>pval]<-NA
  
  if(plot){
    library(ggplot2)
    resPlot<-resCorrigido
    resPlot[is.na(resPlot)]<-0
    #resPlot<-log10(res)
    p<-ggplot()+
      theme_bw()+
      xlab("Generations (x10,000)")+
      scale_y_continuous(trans='sqrt')
    if(normalize){
      p<-p+ylab("p-value (Normalized set)")
    }else{
      p<-p+ylab("p-value")
    }
    i=1
    colors<-rainbow(n=quant)
    max<-max(resPlot)
    
    min<-max/3
    #colors<-c("red","blue","green","orange")
    legenda<-data.frame(x=rep(0,quant),
                        y=seq(max,min,length.out = quant),
                        min=vdw$min,
                        max=vdw$max,
                        color=colors)
    for(i in 1:quant){
      qt1<-data.frame(x=seq(0,60.5,0.5),
                      y=t(resPlot[i,]))
      colnames(qt1)<-c("x","y")
      p<-p+geom_point(data = qt1[qt1$y!=0,], aes(x,y),col=colors[i],pch=15)
      p<-p+geom_point(data = legenda[i,],aes(x,y),col=colors[i])
      p<-p+geom_text(data = legenda[i,],aes(x+5,y,label=paste(min,"to",max)),cex=3)
    }
    plot(p)
    
  }
  
  
  return(resCorrigido)
  #log(1)
  
  #elimina zeros
  # result[result==0]<-NA
  # faixas<-result[,1]
  # result<-result[,-1]
  # #remove linhas que só contém NA
  # tmp1<-data.frame(t(!is.na(t(result))))
  # tmp1<-apply(X = tmp1,MARGIN = 2,sum)
  # sum(tmp1!=0)-1
  # tmp2<-result[tmp1!=0]
  # final<-data.frame(age=as.numeric(sub(pattern = "X",
  #                           replacement = "",
  #                           x = colnames(tmp2))),
  #                   padj=tmp2[1,], stringsAsFactors = F)
  # #tmp<-data.frame(t(result[,-1]))
  # 
  # #result$limites<-quantil[1:length(quantil)-1]
  # teste<-data.frame(age=seq(0,60500,500))
  # final<-merge(final,teste,by="age",all=T)
  # final$padj[is.na(final$padj)]<-0
  # 
  # library(ggplot2)
  # ggplot()+
  #   scale_x_continuous(name = "Ages (x10000)",breaks = seq(0,60000,10000) )+ 
  #   ylab("p-value")+
  #   ylim(-0.001,max(final$padj))+
  #   geom_point(data = final[final$padj!=0,],
  #              mapping = aes(x = age,y = padj),
  #              col="red", pch=20)+
  #   geom_point(data = final[final$padj ==0,],
  #              mapping = aes(x = age,y = padj),
  #              col="blue", pch=15)+
  #   theme_bw()
  # 
  # #result<-data.frame(t(result))
  # 
  # # #result<-merge(result,deltaW[,c(1,2,4)],by=c("faixa"))
  # # outW<-deltaW[,c(1,2,4)]
  # # outW$concat<-paste0(outW$ref,":",outW$mut)
  # # outW<-outW[!outW$concat%in%paste0(result$ref,":",result$mut),1:3]
  # # cat("dw negativos não significativos: ", nrow(outW[outW$Yoon<=0,]))
  # # cat("dw positivos não significativos: ", nrow(outW[outW$Yoon>0,]))
  # # 
}

res<-enrich(normalize = T, 
            plot = T, 
            pval = 0.001,
            rank = 100,
            type = "High", 
            tail = "L")
