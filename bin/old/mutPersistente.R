rm(list = ls())

type="High"
pval=0.01
quant=7
normalize = T
plot = T
rank=100
tail="L"


plotCorr<- function(type="High"){
  setwd("/home/clovis/Doutorado/Projetos/Ecoli60/data_files/")
  deltaW<-read.csv("cDeltaW.csv",
                   header = T,
                   stringsAsFactors = F)
  codonUsage<-read.csv("dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  
  tmp<-merge(codonUsage[,c("ref","count")],deltaW[,c("ref","mut","LeGac")],by="ref")
  tmp<-merge(tmp,codonUsage[,c("ref","count")],by.x="mut",by.y="ref")
  colnames(tmp)<-c("mut","ref","countRef","dw","countMut")
  tmp<-tmp[,c("ref","mut","countRef","countMut","dw")]
  tmp$codonRatio<-log10(tmp$countMut/tmp$countRef)
  lm<-lm(tmp$codonRatio~tmp$dw)
  corr<-cor.test(tmp$codonRatio,tmp$dw)
  p<-ggplot()+theme_bw()+
    xlab(expression(paste(Delta,"w")))+
    ylab("Codon Usage Ratio (log10)")+
    geom_point(data = tmp,
               aes(x=dw,y=codonRatio), 
               col="red", 
               pch=1)+
    geom_abline(slope = lm$coefficients[2], 
                intercept = lm$coefficients[1],
                lty=2,
                col="blue")+
    geom_text(aes(x=-1, y=1,
                  label=paste("Correlation:",
                              round(corr$estimate,3),
                              "\npvalue:",
                              format(corr$p.value,digits = 3, scientific=T))))
  p
  return(p)
}

mutPersist<- function(type="High",
                  pval=0.01,
                  quant=7, 
                  normalize = F, 
                  rank=100, 
                  tail="L"){
  if(tail%in%c("L","H")){
    tail<- (tail == "L")
  }else{
    cat('Use values "L" (Low) or "H" (High) for tail')
    return(0)
  }
  setwd("/home/clovis/Doutorado/Projetos/Ecoli60/data_files/")
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
  
  
  #remove stop codons
  mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"),]
  #ordena por posição
  mutAllel<-mutAllel[order(mutAllel$Position),]
  
  codons<-unique(mutAllel[,c("ref","mut")])
  
  # #transforma zeros em NA
  # m<-data.frame((apply(X = mutAllel,
  #       MARGIN = 2, FUN = function(x){
  #         x[x==0]<-NA
  #         return(x)
  #       }
  # )))
    
  genes<-unique(mutAllel$Gene[order(mutAllel$Gene)])
  gene<-genes[1]
  tmp<-data.frame(age=seq(1,122,1),freq=t(mutAllel[mutAllel$Gene==gene,5:126] ))
  tmp<-data.frame((apply(X = tmp,
        MARGIN = 2, FUN = function(x){
          x[x==-1]<-0
          return(x)
        }
  )))
  tmp$mean<-apply(X = tmp,
                 MARGIN = 1, FUN = function(x){
                   soma<-mean(x[2:length(x)])
                   return(soma)
                   }
                 )
  
  tmp<-tmp[tmp$mean!=0,]
  
  top100Base<-mutAllel[,c(1,2,127,128,3)]
  colnames(top100Base)<-c("Position","Gene","ref","mut","Allele")
  top100Base<-merge(top100Base,deltaW[,c(1,2,4)],by=c("ref","mut"))
  top100Base<-merge(top100Base,codonUsage[,c(1,4)],by=c("ref"))
  colnames(top100Base)<-c("ref","mut","Position","Gene","Allele","dw", "factor")

  #garante a mesma ordem me mutAllel e top100Base
  top100Base<-top100Base[order(top100Base$Position),c("Position","Gene","Allele","ref","mut","dw","factor")]

  if(exists("result")){rm(result)}
  col="X33500"
  for( col in colnames(mutAllel[5:126])){
    top100<-cbind(top100Base,mutAllel[,col])
    colnames(top100)<-c("Position","Gene","Allele","ref","mut","dw","factor","count")
    
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
  return(list(resCorrigido,vdw))
}

plotEnrDepl<-function(type="High",
                      pval=0.01,
                      quant=7, 
                      normalize = F, 
                      rank=100){
  library(plotly)
  
  resCorrigido<-enrich(normalize = normalize, 
                       pval = 0.001,
                       rank = 100,
                       type = "High", 
                       tail = "L")
  vdw<-resCorrigido[[2]]
  resCorrigido<-resCorrigido[[1]]
  
  resPlot<-resCorrigido[,1:121]
  resPlot[is.na(resPlot)]<-0
  colors<-rainbow(n=quant)
  titulo<-paste("Frequency of codon mutation",
                ifelse(normalize,"- normalized","- without normalization"))
  
  legenda<-data.frame(x=rep(0,quant+1),
                      y=rep(0,quant+1),
                      DW=c(vdw$min,vdw$max[quant]),
                      max=c(vdw$max,vdw$max[quant]),
                      color=c(colors,colors[quant]))
    p<-plot_ly()
    p<-add_trace(p=p,
                 type="scatter",
                 mode="markers",
                 data = legenda,
                 x=~x,  
                 y=~DW,  
                 name = "m",
                 color = ~DW,
                 colors = colors,
                 visible = "legendonly")%>%
      layout(xaxis = list(title = "k Generatios"),
             yaxis = list(title = "",
                          showticklabels=F,
                          showline=F,
                          showgrid=F,
                          range = c(-1,1)),
             title = titulo)
    i=1
    linha = c(-0.5,-0.6)
    contLin = 0
    for(i in 1:quant){
      qt1<-data.frame(x=seq(0,60,0.5),
                      y=t(resPlot[i,]))
      colnames(qt1)<-c("x","y")
      if(nrow(qt1[qt1$y!=0,])>1){
        contLin<-contLin + 1
        qt1$y[qt1$y!=0] <- linha[contLin]
        p<-add_trace(p = p, 
                     y = qt1$y[qt1$y!=0],                  
                     x = qt1$x[qt1$y!=0] , 
                     type="scatter", 
                     mode="markers",
                     marker=list(color = colors[i], width = 0.5, symbol = "square") )%>%
          layout(showlegend = FALSE)
      }
    }

    resCorrigido<-enrich(normalize = normalize, 
                         pval = 0.001,
                         rank = 100,
                         type = "High", 
                         tail = "H")
    resCorrigido<-resCorrigido[[1]]
    
    resPlot<-resCorrigido[,1:121]
    resPlot[is.na(resPlot)]<-0
    i=1
    linha = c(0.5,0.6)
    contLin = 0
    for(i in 1:quant){
      qt1<-data.frame(x=seq(0,60,0.5),
                      y=t(resPlot[i,]))
      colnames(qt1)<-c("x","y")
      if(nrow(qt1[qt1$y!=0,])>1){
        contLin<-contLin + 1
        qt1$y[qt1$y!=0] <- linha[contLin]
        p<-add_trace(p = p, 
                     y = qt1$y[qt1$y!=0],                  
                     x=qt1$x[qt1$y!=0] , 
                     type="scatter", 
                     mode="markers",
                     marker=list(color = colors[i], size = 5, symbol = "square") )%>%
          layout(showlegend = FALSE)
      }
    }
    
    texto <- c('Depleted', 'Enriched')
    x <- c(5, 5)
    y <- c(-0.4, 0.4)
    data <- data.frame(texto, x, y)
    
    p <- add_trace(p=p, 
                   data = data, 
                   x = ~x, 
                   y = ~y, 
                   type = 'scatter',
                   mode = 'text', 
                   text = ~text, 
                   textposition = 'middle right',
                   textfont = list(color = '#000000', size = 12)) 
  #p
  return(p)
}

p<-plotCorr(type="High",
                         pval=0.001,
                         quant=7, 
                         normalize = T, 
                         rank=100)
  p

  plotCorr(type="High")
  