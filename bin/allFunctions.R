# !diagnostics off
library(dplyr)
calcDwOld<-function(top=86,
                    workdir){
  setwd(workdir)
  deltaW<-readDeltaW(top,workdir)
  deltaW$dw<-NA
  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  aminos<-read.csv("AuxFiles/dAminoCodon.csv",
                       header = T,
                       stringsAsFactors = F)
  codonUsage<-merge(codonUsage[,c("ref","count")],
                    aminos, 
                    by="ref")
  #deltaW<-merge(deltaW,codonUsage,by="ref")
  amino="!"
  codonUsage$w<-NA
  for (amino in unique(aminos$amino[!aminos$amino%in%c("!","M","W")])){
    tmp<-codonUsage[codonUsage$amino==amino,]
    codonUsage<-codonUsage[!codonUsage$amino==amino,]
    #cat(amino)
    maxCount<-max(tmp$count)
    tmp$w<-tmp$count/maxCount
    
    codonUsage<-rbind(codonUsage,tmp)
  }
  tmp<-merge(deltaW,codonUsage[,c("ref","w")],
             by="ref", all.x =) 
  colnames(tmp)<-c("ref","mut","dw","wRef")
  tmp<-merge(tmp,codonUsage[,c("ref","w")],
             by.x = "mut",
             by.y = "ref",
             all.x = T)[,c("ref","mut","dw","wRef","w")]
  
  colnames(tmp)<-c("ref","mut","dw","wRef","wMut")
  tmp$dw<-apply(tmp[,c("wRef","wMut")],
                MARGIN = 1,
                FUN = function(x){
                  return(log10(x[2]/x[1]))
                  })
  return(tmp[,c("ref","mut","dw")])
  
}


plotCorrCaiOld<- function(top=86,
                    save=F,
                    workdir){
  cat("Using",workdir)
  #setwd(workdir)
  deltaW<-readDeltaW(top,workdir)
  deltaWOld<-calcDwOld(top,workdir)
  colnames(deltaWOld)<-c("ref", "mut", "dwOld" )

  tmp<-merge(deltaW,deltaWOld,by=c("ref","mut"))

  lm<-lm(tmp$dwOld~tmp$dw)
  corr<-cor.test(tmp$dwOld,tmp$dw)
  maxX<-max(tmp$dw)
  maxY<-max(tmp$dwOld)
  cat(paste("\nDelta W \n\tfrom\t",
            round(-maxX,2),"to",
            round(maxX,2),
            "\n\tand\t",round(-maxY,2),"to",
            round(maxY,2),
            "\n"))
library(ggplot2)
  g<-ggplot()+theme_bw()+
    xlab(bquote(Delta~"w New Method"))+
    ylab(bquote(Delta~"w Old Method"))+
    geom_point(data = tmp,
               aes(x=dw,y=dwOld), 
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
  print(g)
  if(save){
    ggsave(filename = "../figures/corrDwNewOld.pdf", 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 0.6,
           width = 16.610, height = 7.440, units = "in",
           dpi = 300)
  }
  
  #return(p)
}

mutPersist<- function(top=86,
                      type="High",
                  pval=0.01,
                  quant=5, 
                  normalize = F, 
                  rank=100, 
                  tail="L",
                  workdir){
  if(tail%in%c("L","H")){
    tail<- (tail == "L")
  }else{
    stop('Use values "L" (Low) or "H" (High) for tail')
  }
  setwd(workdir)
  deltaW<-readDeltaW(top,workdir)
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
  mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"|mutAllel$ref == "TAG" ),]
  #ordena por posição
  mutAllel<-mutAllel[order(mutAllel$Position),]
  
  codons<-unique(mutAllel[,c("ref","mut")])
  
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
    top<-top100$cNorm[order(top100$cNorm, decreasing = T)]
    if(top[rank] == 0){
      top100<-top100[top100$cNorm>top[rank],]
    }else{
      top100<-top100[top100$cNorm>=top[rank],]
    }
    counts100<-as.data.frame(table(top100[,c("faixa")]),stringsAsFactors = F)
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


joinHighLow<- function(workdir){
  mutAllelAll<-read.csv(paste0("./base/aHighMutAllDW.csv"),
                        header = T,
                        stringsAsFactors = F)
  mutAllelAllL<-read.csv(paste0("./base/aLowMutAllDW.csv"),
                         header = T,
                         stringsAsFactors = F)
  all<-rbind(mutAllelAll,mutAllelAllL)
  write.csv(all,
            file = "./base/aHighLowMutAllDW.csv",
            row.names = F,
            workdir=workdir)
}


totalCodons<-function(top=86,
                      workdir){
  setwd(workdir)
  deltaW<-read.csv(paste0("AuxFiles/cDeltaW",top,".csv"),
                   header = T,
                   stringsAsFactors = F)
  
  types<-c("Low","High","MutT","HighLow")
  #type<-"HighLow"
  for(type in types){
    mutAllelAll<-read.csv(paste0("base/a",type,"MutAllDW.csv"),
                          header = T,
                          stringsAsFactors = F)
    #só syn
    mutAllel<-mutAllelAll[mutAllelAll$Annotation == "synonymous",]
    #remove stop codons
    mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"| mutAllel$ref == "TAG"),]
    
    codons<-unique(mutAllel[,c("ref","mut")])
    
    somat<-mutAllel[0,c(1:2,7:128)]
    probab<-mutAllel[0,c(1:2,7:128)]
    #i=2
    for(i in 1:nrow(codons)){
      tmp<-subset(x= mutAllel,subset = (ref == codons$ref[i] & 
                                          mut == codons$mut[i]),
                  select = X0:X60500)
      #remove NA
      tmp[is.na(tmp)]<-0

      tmp2<-apply(X = tmp,MARGIN = 2, FUN = function(x){
        x<-x[x!=-1]
        x<-1-x
        return(1-prod(x))
      })

      probab[nrow(probab)+1,1:2]<-c(codons$ref[i],codons$mut[i])    
      
      probab[nrow(probab),3:124]<-as.numeric(tmp2)
      
      tmp2<-apply(X = tmp,MARGIN = 2, FUN = function(x){
        x<-x[x!=-1]
        return(sum(x))
      })
      somat[nrow(somat)+1,1:2]<-c(codons$ref[i],codons$mut[i])
      somat[nrow(somat),3:124]<-as.numeric(tmp2)
    } 
    totais<- apply(somat[,3:124],2,sum)
    MIN = min(totais[totais!=0])
    cat("Tipo:", type, " Mínimo = ",min(MIN), " Máximo = ", max(totais),"\n")
    totais<- round(totais*1/MIN,0)
    
    probab<-merge(probab,deltaW[,c(1,2,3)],by=c("ref","mut"))
    #colnames(probab)[colnames(probab) == 'LeGac'] <- 'dW'
    somat<-merge(somat,deltaW[,c(1,2,3)],by=c("ref","mut"))
    #colnames(somat)[colnames(somat) == 'LeGac'] <- 'dW'
    
    
    write.csv(x = probab,
              file = paste0("AuxFiles/bProb",type,".csv"),
              row.names = F)
    write.csv(x = somat,
              file = paste0("AuxFiles/bSum",type,".csv"),
              row.names = F)
    write.csv(x = totais,
              file = paste0("AuxFiles/bTot",type,".csv"),
              row.names = F)
    
    
  }
}

mutationsMutT<-function(workdir){
  setwd(workdir)
  deltaW<-readDeltaW(86,workdir)
  # Cria arq com mutações MutT
  codons<- unique(deltaW$ref)
  codons<-data.frame(ref=codons[grep("A|T",codons)],stringsAsFactors = F)
  codons<-merge(codons,deltaW[,1:2],by="ref")
  mutT<- codons[0,]
  for(i in 1:nrow(codons)){
    ref<-unlist(strsplit(as.character(codons$ref[i]),split = "",fixed = T))
    mut<-unlist(strsplit(as.character(codons$mut[i]),split = "",fixed = T))
    if((ref[1]=="A" & mut[1]=="C")|
       (ref[1]=="T" & mut[1]=="G")|
       (ref[2]=="A" & mut[2]=="C")|
       (ref[2]=="T" & mut[2]=="G")|
       (ref[3]=="A" & mut[3]=="C")|
       (ref[3]=="T" & mut[3]=="G")){
      mutT[nrow(mutT)+1,]<-codons[i,]
    }
  }
  write.csv(x = mutT,
            file = "AuxFiles/dMutTmutations.csv",
            row.names = F)
  
}

joinHighLow<- function(workdir){
  setwd(workdir)
  mutAllelAll<-read.csv(paste0("./base/aHighMutAllDW.csv"),
                        header = T,
                        stringsAsFactors = F)
  mutAllelAllL<-read.csv(paste0("./base/aLowMutAllDW.csv"),
                         header = T,
                         stringsAsFactors = F)
  all<-rbind(mutAllelAll,mutAllelAllL)
  write.csv(all,
            file = "./base/aHighLowMutAllDW.csv",
            row.names = F)
  
  mutAllelAll<-read.csv(paste0("./base/aHighMutSynFreqDW.csv"),
                         header = T,
                         stringsAsFactors = F)
  mutAllelAllL<-read.csv(paste0("./base/aLowMutSynFreqDW.csv"),
                         header = T,
                         stringsAsFactors = F)
  all<-rbind(mutAllelAll,mutAllelAllL)
  write.csv(all,
            file = "./base/aHighLowMutSynFreqDW.csv",
            row.names = F)
}

enrich<- function(type="High",
                  pval=0.001,
                  quant=5, 
                  normalize = T, 
                  normBy = "count",
                  rank=200, 
                  tail="L",
                  fix= T,
                  TsTv = "Ts",
                  Dw = "Cai",
                  top = 86,
                  quiet=T,
                  workdir ){
  library(dplyr)
  
  if(tail%in%c("L","H")){
    tail<- (tail == "L")
  }else{
    stop('Use values "L" (Low) or "H" (High) for tail')
    #return(0)
  }
  if(!TsTv%in%c("all","Ts","Tv")){
    stop('Use values "all","Ts" (transitions), or "Tv" (transversions) for TsTv')
    #return(0)
  }
  if(!Dw%in%c("Tai","Cai")){
    stop('Use values "all","Ts" (transitions), or "Tv" (transversions) for TsTv')
  }
  if(normalize){
    if(!normBy%in%c("CUB","count","mean")){
      stop('Normalizations can be done only by Codon Usage Bias ("CUB"), mutation count ("count"), or mean of mutations in a generation')
    }
  }
  vTsTv<-TsTv
  rankStd<-rank
  
  
  
  setwd(workdir)
  if(Dw == "Tai"){
    deltaW<-read.csv("./AuxFiles/cDeltaWTAI.csv",
                     header = T,
                     stringsAsFactors = F)
    #top is ignored when is a Tai Dw
    top<-''
  }else if(Dw == "Cai"){
    deltaW<-readDeltaW(top,workdir)
  }
  #round delta W values in 1e-10 to avoid some rounds problems
  deltaW$dw<-round(deltaW$dw,10)
  codonUsage<-read.csv("./AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  
  mutAllelAll<-read.csv(paste0("./base/a",type,"MutAllDW.csv"),
                        header = T,
                        stringsAsFactors = F)
  
  #só syn
  mutAllel<-mutAllelAll[mutAllelAll$Annotation == "synonymous",]
  #remove stop codons
  mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"| mutAllel$ref == "TAG"),]

  #count all mutations tha have frequency greater than zero
  qtCount<-mutAllel[,c(1,7:128)]
  qtCount[is.na(qtCount)]<-0
  qtCount$sum<-apply(qtCount[,-1],
                     MARGIN = 1,
                     FUN = sum)
  qtCount<-qtCount[qtCount$sum !=0,c("ref","sum")]
  
  qtCount<-qtCount%>%
    group_by(ref)%>%
    summarize(tmpCol=n())
  
  #atualize dw values 
  mutAllel<-merge(mutAllel,deltaW,
                  by=c("ref","mut"))
  mutAllel$dw.x<-NULL
  colnames(mutAllel)[135]<-"dw"
  #filtra Ts ou Tv
  if(TsTv != "all"){
    TsTvdf<-read.csv("AuxFiles/TsTv.csv",
                   header = T,
                   stringsAsFactors = F)
    deltaW<-merge(deltaW, 
                  TsTvdf[c("ref","mut","type")],
                  by=c("ref","mut"))
    deltaW<-subset(deltaW, type == vTsTv)
    
    mutAllel<-subset(mutAllel, TsTv == vTsTv)
  }
  
  if(!fix%in%c("all",T,F) ){
    cat('Argument "fix" must be TRUE, FALSE or "all".')
    return()
  }
  if(fix!="all"){
    if(fix){
      mutAllel<-mutAllel[!is.na(mutAllel$fixed),]
    }else{
      mutAllel<-mutAllel[is.na(mutAllel$fixed),]
    }
  }
  #ordena por posição
  mutAllel<-mutAllel[order(mutAllel$Position),]
  
  codons<-unique(mutAllel[,c("ref","mut")])
  
  if(rank<0){
    stop("Rank must be an integer or a number between 0 and 1.")
    #return(0)
  }else if(rank>=0 & rank <=1){
    rankStd<- round(as.numeric(nrow(mutAllel))*as.numeric(rank),0)
  }else{
    rankStd<- round(rank,0)
  }
  if(!quiet){
    cat("Using the top",rankStd,"frequencies \n")
  }
  
  
  top100Base<-mutAllel[,c(1:5,135)]#,2,127,128,3)]
  top100Base<-merge(top100Base,codonUsage[,c(1,3)],by=c("ref"))
  colnames(top100Base)[colnames(top100Base) == 'factorPerK'] <- 'factor'

  #quantilize aqui ----
  quantilizedDw<-quantilizeDw(deltaW = deltaW[,c("ref","mut","dw")],
               quant = quant,
               Dw = Dw,
               quiet = quiet)
  vdw<-quantilizedDw[[2]]
  faixasDw<-quantilizedDw[[1]]
  faixasDw$dw<-NULL
  top100Base<-merge(top100Base,faixasDw,
                    by = c("ref","mut"))

  #garante a mesma ordem me mutAllel e top100Base
  top100Base<-top100Base[order(top100Base$Position),c("Position","Gene","Allele","ref","mut","dw","factor","faixa")]

  
  if(exists("result")){rm(result)}
  #loop ----
  col="X25000"
  faixas<-seq(1:quant)
  minFreq=Inf
  for( col in colnames(mutAllel[7:128])){
    rank<-rankStd
    if(fix=="all"){
      top100<-cbind(top100Base,mutAllel[,col])
    }else if(fix){
      tmpTime<-mutAllel[,col]
      curTime<-strtoi(substr(col,2,nchar(col)))
      tmpTime[curTime<mutAllel$t0]<-NA
      top100<-cbind(top100Base,tmpTime)
    }else{
      top100<-cbind(top100Base,mutAllel[,col])
    }
    colnames(top100)<-c("Position","Gene","Allele","ref","mut","dw","factor","faixa","count")
    #remove NA and -1 frequences and zeros
    top100$count[top100$count == -1] <- NA
    top100<-na.exclude(top100)
    top100<-top100[top100$count>0,]
    #top100$count<-top100$count-1e-10
    if(nrow(top100)==0){
      next()
    }

    if(normalize){
      if(normBy =="CUB"){
      
      #normalize observations
      top100$cNorm<-apply(X = top100[,c("factor","count")],
                          MARGIN = 1,
                          FUN = function(x){
                            return(x[2]/x[1])
                          })
      }else if(normBy =="count"){
        top100<-merge(top100,
                        qtCount, 
                        by="ref",
                      all.x=T)
        top100$cNorm<-apply(X = top100[,c("tmpCol","count")],
                            MARGIN = 1,
                            FUN = function(x){
                              return(x[2]/x[1])
                            })
        top100$tmpCol<-NULL
      }else if(normBy =="mean"){
        #count all mutations that have frequency greater than zero in generation
        qtCountTop<-top100[,c("ref","count")]

        qtCountTop<-qtCountTop%>%
          group_by(ref)%>%
          summarize(tmpCol=n())

        top100<-merge(top100,
                      qtCountTop, 
                      by="ref",
                      all.x=T)
        top100$cNorm<-apply(X = top100[,c("tmpCol","count")],
                            MARGIN = 1,
                            FUN = function(x){
                              return(x[2]/x[1])
                            })
        top100$tmpCol<-NULL
        
        
        
      }
    }else{
      top100$cNorm<-top100$count
    }
    countsTot<-as.data.frame(table(top100[c("faixa")]),
                             stringsAsFactors = F)
    countsTot<-countsTot[countsTot$Freq!=0,]
    colnames(countsTot)<- c("faixa","white")
    countsTot$black<-nrow(top100)-countsTot$white
    
    #colnames(top100)<-c("Position","Gene","Allele","ref","mut","dw","factor","faixa","count")
    top<-top100$cNorm[order(top100$cNorm, decreasing = T)]
    if(rank > length(top))
      rank <- length(top)
    if(top[rank] == 0){
      top100<-top100[top100$cNorm>top[rank],]
    }else{
      top100<-top100[top100$cNorm>=top[rank],]
    }
    if(nrow(top100)==0){
      next()
    }
    minFreq<-min(minFreq,top[rank])
    counts100<-as.data.frame(table(top100[,c("faixa")]),stringsAsFactors = F)
    #counts100<-counts100[counts100$Freq!=0,]
    colnames(counts100)<-c("faixa","Freq")
    #completa a tabela de counts com zeros
    for (padrao in (faixas[!faixas%in%counts100$faixa])) {
      counts100[nrow(counts100)+1,]<-c(padrao,0)
    }
    if(sum(counts100$Freq)-nrow(top100) !=0){
      cat(col," fail!","\n")
      next
    }
    counts100<-merge(counts100,countsTot,by=c("faixa"))
    counts100$drawn<-nrow(top100)
    counts100$hyp<-phyper(counts100$Freq,
                          counts100$white,
                          counts100$black,
                          counts100$drawn,
                          lower.tail = tail)
    #exclude non significant zeros
    #occured when drawn all balls avaliable on upper tail test
    counts100$hyp[counts100$white+counts100$black == counts100$drawn &
                    counts100$hyp == 0]<- NA
    counts100<-counts100[, c("faixa","hyp")]
    #adjust p-val
    counts100$hyp<-p.adjust(counts100$hyp)
    colnames(counts100)<-c("faixa",col)
    if(exists("result")){
      result<-merge(result,counts100,by=c("faixa"), all = T)
    }else{
      result<-counts100
    }
  }
  if(!quiet){
    cat("Minimal frequence used:", minFreq,"\n")
  }

  resCorrigido<-result[,2:ncol(result)]

  resCorrigido[resCorrigido>pval]<-NA
  return(list(resCorrigido,vdw))
}

plotEnrDeplPVal<-function(type="High",
                          pval=0.001,
                          quant=5, 
                          normalize = T, 
                          normBy = "mean",
                          rank=200,
                          fix = T,
                          TsTv = "all",
                          title = T,
                          save=F,
                          Dw = "Tai",
                          top=200,
                          figName = '',
                          workdir ){
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)  
  setwd(workdir)
  if(!file.exists(paste0("./base/a",type,"MutAllDW.csv"))){
    if(type == "HighLow"){
      joinHighLow(workdir)
    }else{
      stop('Parameter "type" must be 
          \t"High" for High mutations rate populations;
          \t"Low" for Low mutations rate populations;
          \t"MutT" for MutT mutations populations; or
          \t"HighLow" for High and Low mutations together.' )
      #return(0)
    }
  }
  
  rankStd<-rank
  if(Dw == "Tai"){
    top<-""
  }else{
    if(!is.numeric(top)){
      stop('Parameter "top"must be numeric')
    }
  }

  if(title){
    titulo<-paste("Top",rankStd,"Frequencies, ",
                  ifelse(normalize,paste0("norm by ",normBy,", "), "without norm, "),
                  type, ", ",toupper(Dw),top,", ")
  }else{
    titulo=" "
  }
  if(!fix%in%c("all",T,F) ){
    cat('Argument "fix" must be TRUE, FALSE or "all".')
    return()
  }
  if(fix!="all"){
    if(fix){
      titulo<-paste(titulo,"Fixed Mutations")
      name<-"Fix"
    }else{
      titulo<-paste(titulo,"Non-fixed Mutations")
      name<-"Nonfix"
    }
  }else{
    titulo<-paste(titulo,"All Mutations")
    name<-"All"
  }
  
  if(!TsTv%in%c("all","Ts","Tv")){
    stop('Use values "all","Ts" (transitions), or "Tv" (transversions) for TsTv')
    #return(0)
  }
  
  resCorrigido<-enrich(normalize = normalize, 
                       normBy = normBy,
                       pval = pval,
                       rank = rankStd,
                       type = type, 
                       tail = "L",
                       quant = quant,
                       fix = fix,
                       TsTv = TsTv,
                       Dw = Dw,
                       top = top,
                       quiet = F,
                       workdir = workdir)

  vdw<-resCorrigido[[2]]
  resCorrigido<-resCorrigido[[1]]
  
  resPlot<-resCorrigido[,1:ncol(resCorrigido)]
  #adiciona indice de faixas
  resPlot$faixa<- c(-1:-quant)
  #resPlot$faixa<- -(quant+1)+c(1:quant)
  #elimina linhas sem resultado
  #resPlot<-resPlot[rowSums(!is.na(resPlot[,1:ncol(resCorrigido)])) > 0,]
  resCorrigido<-enrich(normalize = normalize, 
                       normBy = normBy,
                       pval = pval,
                       rank = rankStd,
                       type = type, 
                       tail = "H",
                       quant = quant,
                       fix = fix,
                       TsTv = TsTv,
                       Dw = Dw,
                       top = top,
                       workdir = workdir)
  vdw<-resCorrigido[[2]]
  resCorrigido<-resCorrigido[[1]]
  
  tmp<-resCorrigido[,1:ncol(resCorrigido)]
  #elimina 0
  #tmp[tmp==0]<-NA #mudado pq a correção agora acontece dentro do enrich()
  #adiciona indice de faixas
  tmp$faixa<-c(1:quant)
  
  #elimina linhas sem resultado -----
  #tmp     <-tmp[     rowSums(!is.na(tmp)) > 0,]
  #enriched<-enriched[rowSums(!is.na(enriched)) > 0,]
  
  #une resultados
  resPlot<-rbind(resPlot, tmp)
  
  #elimina colunas sem resultado
  #resPlot<-resPlot[colSums(!is.na(resPlot)) > 0]
  #test if some enrich/depletions was detected
  if(sum(dim(resPlot))== 0){
    warning("No enrichments or depletions detected!")
    return(invisible(0))
  }else{
    tmp<-resPlot
    tmp$faixa<-NULL
    if(sum(!is.na(tmp)) == 0 ){
      warning("No enrichments or depletions detected!")
      return(invisible(0))
    }
  }
  
  
  min<-min(resPlot[!is.na(resPlot) & resPlot >0])
  # range of p-val in powers of ten
  val<-seq(floor(log10(min)),ceiling(log10(pval)),1)
  #Colors for p-values
  mydf <- data.frame(id = rep(1, length(val)), 
                     pval = val,
                     cor=rainbow(length(val),end=0.7 ),#topo.colors
                     #cor=colorRampPalette(c("blue","green", "red"))( length(val) ),
                     stringsAsFactors = F)
  
  #legend bar code
  legenda<-
    ggplot(mydf) +
    ggtitle("p-value",subtitle = "(log10)")+
    theme(title = element_text(size = 6),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size = 6),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())+
    geom_tile(aes(x = 0, y=pval, fill = cor)) +
    scale_fill_identity()
  
  #
  linha =  seq(-0.4,(-quant*0.05)-0.4,-0.05)#c(-0.5,-0.6)
  contLin = 0
  #time line values
  x<-strtoi(sapply(X = colnames(resPlot[,1:(ncol(resPlot)-1)]),
                   substr,2,30))/1000
  i=10
  
  #change real zeros on hipergeometric in the min val to avoid inf
  resPlot[resPlot==0]<-min
  qtAll<-data.frame(x=numeric(),
                    y=numeric(),
                    pval=numeric(),
                    stringsAsFactors = F)
  countReg<-0
  for(i in 1:nrow(resPlot)){
    qt1<-data.frame(x=x,#seq(0,60,0.5),
                    y=rep(resPlot$faixa[i],length(x)),
                    pval=t(round(log10(resPlot[i,1:(ncol(resPlot)-1)]),0)),
                    stringsAsFactors = F)
    colnames(qt1)<-c("x","y","pval")
    countReg<-countReg+nrow(na.exclude(qt1))
    qtAll<-rbind(qtAll,qt1)
    
    # qt1<-merge(qt1,mydf[,2:3], by = "pval")
    # p<-p+geom_point(data = qt1,
    #                 aes(x, y, color = (cor)),
    #                 pch=15)
  }
  qtAll<-merge(qtAll,mydf[,2:3], by = "pval")
  if(countReg == nrow(qtAll)){
    cat("Count ok\n")
  }else{
    cat("Count fail\n")
  }
  
  qtAll$type<-"Enriched"
  qtAll$type[qtAll$y<0]<- "Depleted"
  qtAll$type<-factor(qtAll$type,
                       levels = c("Enriched", "Depleted"))
  qtAll$y[qtAll$y<0]<- -qtAll$y[qtAll$y<0]
  
  library(ggplot2)
  
  
    # ylim(-7,7)+
    xlab("Generations (x1000)")+
    scale_y_continuous(name =bquote(Delta~"w levels"), 
                       label=c(c(1:quant)," ",c(1:quant)),
                       #label=c(c(quant:1)," ",c(1:quant)),
                       #"7","6","5","4","3","2","1"," ",
                       #        "1","2","3","4","5","6","7"),
                       breaks= c(-quant:quant),
                       minor_breaks = NULL,
                       limits = c(-quant,quant))+
    
    geom_hline(yintercept = 0,col="white")+
    geom_hline(yintercept = 0.1,col="black")+
    geom_hline(yintercept = -0.1,col="black")+
    theme(#panel.grid.major.y = element_blank(), 
      panel.grid.minor.y = element_blank(),
      axis.ticks.y = element_blank())
  
  
  
  p<-ggplot()+
    ggtitle(titulo) +
    theme_bw()+
    scale_color_identity()+
    theme(legend.title = element_blank())+
    facet_grid(type~.)+
    scale_x_continuous(name = "Generations (x1000)",
                       limits = c(0,max(x)))+
    scale_y_continuous(name =bquote(Delta~"w levels"), 
                       minor_breaks = NULL,
                       breaks= c(1:quant),
                       limits = c(0.5,quant+0.5))+
    geom_point(data = qtAll,
                aes(x=x,y=y,
                    #group=faixa, 
                    color=as.factor(cor)),
               pch=15)
  print(p)
  
  txt<-textGrob("p-value\n (log10)",
                gp=gpar(fontsize=6))
  
  lay <- rbind(c(1,1,1,1,1,1,1,1,1,1,NA),
               c(1,1,1,1,1,1,1,1,1,1,NA),
               #c(1,1,1,1,1,1,1,1,1,1,3),
               c(1,1,1,1,1,1,1,1,1,1,2),
               c(1,1,1,1,1,1,1,1,1,1,2),
               c(1,1,1,1,1,1,1,1,1,1,NA),
               c(1,1,1,1,1,1,1,1,1,1,NA))
  g<-grid.arrange(p,legenda, layout_matrix = lay)
  
  if(save){
    if(figName == ''){
      filename <- paste0("../figures/enrich",type,rankStd,
                        name,
                        Dw,".pdf")
    }else{
      filename <-paste0("../figures/",figName,".pdf")
    }
    ggsave(filename = filename, 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 0.6, 
           width = 16.610, height = 7.440, units = "in",
           dpi = 300)
  }
  
  
  return(resPlot)
  }



plotEnrDeplAll<-function(type="High",
                          pval=0.001,
                          quant=7, 
                          normalize = F, 
                          normBy = "count",
                          rank=100,
                          fix = "all",
                          TsTv = "all",
                          title = T,
                          save=F,
                          Dw = "Cai",
                          top=200,
                          figName = '',
                          workdir ){
  library(reshape2)
  library(ggplot2)
  library(dplyr)
  
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)  
  setwd(workdir)
  if(!file.exists(paste0("./base/a",type,"MutAllDW.csv"))){
    if(type == "HighLow"){
      joinHighLow(workdir)
    }else{
      stop('Parameter "type" must be 
           \t"High" for High mutations rate populations;
           \t"Low" for Low mutations rate populations;
           \t"MutT" for MutT mutations populations; or
           \t"HighLow" for High and Low mutations together.' )
      #return(0)
    }
  }
  
  rankStd<-rank
  if(Dw == "Tai"){
    top<-""
  }else{
    if(!is.numeric(top)){
      stop('Parameter "top"must be numeric')
    }
  }
  
  if(title){
    titulo<-paste0("Hipergeometric test, Top ",rankStd," Frequencies, ",
                  ifelse(normalize,"normalized, ","without normalization, "),
                  type, ", ",toupper(Dw),top,", ")
  }else{
    titulo=" "
  }
  if(!fix%in%c("all",T,F) ){
    cat('Argument "fix" must be TRUE, FALSE or "all".')
    return()
  }
  if(fix!="all"){
    if(fix){
      titulo<-paste0(titulo,"Fixed Mutations")
      name<-"Fix"
    }else{
      titulo<-paste0(titulo,"Non-fixed Mutations")
      name<-"Nonfix"
    }
  }else{
    titulo<-paste0(titulo,"All Mutations")
    name<-"All"
  }
  
  if(!TsTv%in%c("all","Ts","Tv")){
    stop('Use values "all","Ts" (transitions), or "Tv" (transversions) for TsTv')
    #return(0)
  }
  resCorrigido<-enrich(normalize = normalize, 
                       normBy = normBy,
                       pval = 1,#pval,
                       rank = rankStd,
                       type = type, 
                       tail = "L",
                       quant = quant,
                       fix = fix,
                       TsTv = TsTv,
                       Dw = Dw,
                       top = top,
                       quiet = F,
                       workdir = workdir)
  
  vdw<-resCorrigido[[2]]
  resCorrigido<-resCorrigido[[1]]
  
  depleted<-resCorrigido[,1:ncol(resCorrigido)]
  #adiciona indice de faixas
  #resPlot$faixa<- c(-1:-quant)
  #resPlot$faixa<- -(quant+1)+c(1:quant)
  resPlot<-depleted
  #elimina linhas sem resultado
  #depleted<-depleted[rowSums(!is.na(depleted[,1:ncol(resCorrigido)])) > 0,]
  resCorrigido<-enrich(normalize = normalize, 
                       normBy = normBy,
                       pval = 1,#pval,
                       rank = rankStd,
                       type = type, 
                       tail = "H",
                       quant = quant,
                       fix = fix,
                       TsTv = TsTv,
                       Dw = Dw,
                       top = top,
                       workdir = workdir)
  vdw<-resCorrigido[[2]]
  enriched<-resCorrigido[[1]]
  
  #tmp<-resCorrigido[,1:ncol(resCorrigido)]
  #elimina 0
  #tmp[tmp==0]<-NA
  #adiciona indice de faixas
  #tmp$faixa<-c(1:quant)
  #elimina linhas sem resultado
  enriched<-enriched[rowSums(!is.na(enriched)) > 0,]
  #tmp1<-resPlot
  # tmp1[tmp1 == 0]<-1e-50
  # tmp[tmp == 0]<-1e-50
  #teste<-mapply(min, enriched,depleted)
  
  #resPlot<- log10(depleted/enriched)
  for(lin in 1:nrow(resPlot)){
    for(col in 1:ncol(resPlot)){
      dpl<-ifelse(is.na(depleted[lin,col]),1,depleted[lin,col])
      enr<-ifelse(is.na(enriched[lin,col]),1,enriched[lin,col])
      if(dpl<=enr){
        minP<- dpl
        signal<- 1
      }else{
        minP<- enr
        signal<- -1
      }
      resPlot[lin,col]<-log10(minP)*signal
    }
  }
  #une resultados
  #resPlot<-tmp2
  resPlot$faixa<-c(1:quant)
  resPlot$faixa<-paste("Level",resPlot$faixa)
  #elimina colunas sem resultado
  resPlot<-resPlot[colSums(!is.na(resPlot)) > 0]
  
  if(sum(dim(resPlot))== 0){
    warning("No enrichments or depletions detected!")
    return(invisible(0))
  }
  x<-strtoi(sapply(X = colnames(resPlot[,1:(ncol(resPlot)-1)]),
                   substr,2,30))/1000
  i=5
  qt1<-data.frame(x=numeric(),
                  pval=numeric(),
                  level=character(),
                  color=factor(),
                  stringsAsFactors = F)
  infinites<-data.frame(x=numeric(),
                        pval=numeric(),
                        level=character(),
                        color=factor(),
                        stringsAsFactors = F)
  for(i in 1:quant){
    tmp<-data.frame(x=x,#seq(0,60,0.5),
                    pval=t(resPlot[i,1:(ncol(resPlot)-1)]),
                    stringsAsFactors = F)
    colnames(tmp)<-c("x","pval")
    tmp$level<-resPlot[i,ncol(resPlot)]
    tmp$color<-as.factor(ifelse(tmp$pval >= 0,2,1))
    infinites<-rbind(infinites, tmp[is.infinite(tmp$pval),])
    tmp$pval[is.infinite(tmp$pval)&
               tmp$pval>0]<-max(tmp$pval[is.finite(tmp$pval)])
    tmp$pval[is.infinite(tmp$pval)&
               tmp$pval<0]<-min(tmp$pval[is.finite(tmp$pval)])
    infinites$pval[is.infinite(infinites$pval)&
                     infinites$pval>0]<-max(tmp$pval)
    infinites$pval[is.infinite(infinites$pval)&
                     infinites$pval<0]<-min(tmp$pval)

    qt1<-rbind(qt1,tmp)
  }
  yMax<-max(-min(qt1$pval),max(qt1$pval))
  #qt1<-merge(qt1,mydf[,2:3], by = "pval")
  yLabels<-c(paste0("1e-",floor(yMax)),"1e-10","0.001","0.01","0.05","1",
            "0.05","0.01","0.001","1e-10",paste0("1e-",floor(yMax)))
  yBreaks<-c(-(floor(yMax)),-10,-3,-2,-1.30103,0,
             1.30103,2,3,10,floor(yMax))
  colors<-qt1$color
  library(scales) 
  p<-ggplot()+
    theme_bw()+
    facet_grid(level~.)+ 
    #xlim(0,max(x))+
    # ylim(-7,7)+
    ylab("p-value")+
    xlab("Generations (x1000)")+
    ggtitle(titulo)+
    geom_col(data=qt1,aes(x=x,y = pval, fill=color))+
    scale_y_continuous(limits = c(-yMax,yMax),
                       trans = modulus_trans(1e-10),
                       labels = yLabels,
                       breaks = yBreaks,minor_breaks = F)+
    scale_fill_discrete(name="Tail",labels=c( "Lower","Upper"))+
    theme(axis.text.y = element_text(size=7)) 
    
    if(nrow(infinites)>0){
      p<-p+geom_point(data = infinites,aes(x=x,y=pval,shape = "25"),
                 col = "red", bg = "red", cex = 0.6)+
      scale_shape_manual(name="", values = c("25"=25), 
                         labels = "zeros")
        
    }
  
  print(p)
  if(save){
    if(figName == ''){
      filename <- paste0("../figures/enrichPval",type,rankStd,
                         name,
                         Dw,".pdf")
    }else{
      filename <-paste0("../figures/",figName,".pdf")
    }
    ggsave(filename = filename,
           plot = p,
           device = cairo_pdf, 
           width = 210, 
           height = 297, 
           units = "mm",
           dpi = 300)
  }
  
  return(qt1)
}




timeline<-function(population="High",
                   normalize = T,
                   normBy = "count",
                   type = "box",
                   fix = T,
                   save = F,
                   normType = "perK",
                   mutT = 0,
                   operation = "sum",
                   TsTv = "Ts",
                   workdir ){
  dirFig<-"../figures"
  if(!operation %in% c("mean","sum")){
    stop('Parameter "operation" must be "mean" or "sum"' )
  }
  if(!type %in% c("box","line","bar")){
    stop('Parameter "type" must be "box","line" or "bar"' )
  }
  
  if(!TsTv%in%c("all","Ts","Tv")){
    stop('Use values "all","Ts" (transitions), or "Tv" (transversions) for TsTv')
    #return(0)
  }
  if(normalize){
    if(!normBy%in%c("CUB","count","mean")){
      stop('Normalizations can be done only by Codon Usage Bias ("CUB"), mutation count ("count"), or mean of mutations in a generation')
    }
  }
  
  vTsTv<-TsTv
  
  setwd(workdir)
  library(reshape2)
  library(ggplot2)
  library(dplyr)
  
  deltaW<-read.csv("./AuxFiles/cDeltaWTAI.csv",
                     header = T,
                     stringsAsFactors = F)
  codonUsage<-read.csv("./AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  
  mutAllelAll<-read.csv(paste0("./base/a",population,"MutAllDW.csv"),
                        header = T,
                        stringsAsFactors = F)
  #só Non
  #mutAllel<-mutAllelAll[mutAllelAll$Annotation%in%c("missense","nonsense"),]
  
  #só syn
  mutAllel<-mutAllelAll[mutAllelAll$Annotation == "synonymous",]
  #remove stop codons
  mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"| mutAllel$ref == "TAG"),]
  
  #filtra Ts ou Tv
  if(TsTv != "all"){
    mutAllel<-subset(mutAllel, TsTv == vTsTv)
  }
  
  if(!fix%in%c("all",T,F) ){
    stop('Argument "fix" must be TRUE, FALSE or "all".')
  }
  if(fix!="all"){
    if(fix){
      mutAllel<-mutAllel[!is.na(mutAllel$fixed),]
    }else{
      mutAllel<-mutAllel[is.na(mutAllel$fixed),]
    }
  }
  
  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  codonUsage<-codonUsage[order(codonUsage$factorPerK),]
  
  codonOrder<-cbind(codonUsage[c("ref","factorPerK")], 
                    seq(1,nrow(codonUsage),1))
  colnames(codonOrder)<-c("ref","factor","ordem")
  quant<-nrow(codonOrder)
  colors<-rainbow(n=quant,start = 0.2)
  names(colors)<-1:quant
  fills<-alpha(colors,0.5)
  ticktext<-codonUsage$ref
  tickvals<-c(0:(length(ticktext)-1))

  #number of mutation
  qt<-mutAllel[,c(1,7:128)]
  if(type!="bar"){
    if(normalize){
      qt<-normalization(data = qt,
                        normBy = normBy,
                        melted = T,
                        workdir = workdir)
    }else{
      qt<-normalization(data = qt,
                        normBy = "none",
                        melted = T,
                        workdir = workdir)
    }
    if(operation == "sum"){
      qt<-qt%>%
        group_by(ref, variable)%>%
        summarize(y=sum(value))
    }else{
      qt<-qt%>%
        group_by(ref, variable)%>%
        summarize(y=mean(value))
    }
  }

  qt<-merge(qt,codonOrder,by="ref")
  qt<-na.exclude(qt)
  
  titulo<-paste0(population,
                 " (",
                 vTsTv,
                 " ",ifelse(fix,"fixed",""),
                 ") ",
                 ifelse(normalize,
                        paste0("- normalized by ",normBy),
                        "- without normalization"))
  
  p<-ggplot()+theme_bw()+theme(legend.position = "none")+
    ggtitle(titulo)+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)
  if(type=="line"){
    qt$x<-as.numeric(sub("X","",qt$variable))
    p<-p+geom_line(data=qt,
                   aes(x=x,y=y,colour=as.character(ordem)))+
      scale_x_continuous(name = "k generations")+
      scale_y_continuous(name = "Frequency")
  }else if(type=="box") {
    p<-p+geom_boxplot(data=qt,
                      aes(x=ordem,
                          y=y,
                          fill=as.character(ordem),
                          color=as.character(ordem),
                          group=ordem,
                          alpha = 0.2),
    )+
      theme(axis.text.x = element_text(angle = 90))+
      scale_x_continuous(name = "Codons",
                         breaks = 1:length(codonUsage$ref),
                         labels = codonUsage$ref)+
      scale_y_continuous(name = "Frequency")
  }else if(type=="bar"){
    if(normalize){
      warning('Normalize have no effect with "bar" type option')
    }
    #number of mutation
    qt<-mutAllel[,c(1,7:128)]
    qt[is.na(qt)]<-0
    qt$sum<-apply(qt[,-1],
                  MARGIN = 1,
                  FUN = sum)
    #just exclude lines with no frequencies
    qt<-qt[qt$sum !=0,c("ref","sum")]
    
    qt<-qt%>%
      group_by(ref)%>%
      summarize(count=n())
    
    qt<-merge(qt,codonOrder,by="ref")
    qt<-na.exclude(qt)
    qtFactors<-qt
    colnames(qt)<-c("ref","y","factor","ordem")
    colnames(qtFactors)<-c("ref","count","CUB","ordem")
    #count of mutations per age
    col="X0"
    for( col in colnames(mutAllel[7:128])){
      tmpTime<-mutAllel[mutAllel[,col]>0,c("ref",col)]
      tmpTime<-tmpTime%>%
        group_by(ref)%>%
        summarize(count=n())
      colnames(tmpTime)<-c("ref",col)
      qtFactors<-merge(qtFactors,tmpTime, 
                       by= "ref",
                       all.x = T)
    }
    
    p<-ggplot()+theme_bw()+theme(legend.position = "none")+
      #ggtitle(titulo)+
      scale_color_manual(values=colors)+
      scale_fill_manual(values=colors)
    
    p<-p+geom_col(data=qt,
                  aes(x=ordem,
                      y=y,
                      fill=as.character(ordem),
                      color=as.character(ordem),
                      group=ordem,
                      alpha = 0.2),
    )+
      theme(axis.text.x = element_text(angle = 90))+
      scale_x_continuous(name = "Codons",
                         breaks = 1:length(codonUsage$ref),
                         labels = codonUsage$ref)+
      scale_y_continuous(name = "Mutation count")
  }
  
  
  mydf <- data.frame(id = rep(1, length(quant)), 
                     usage = 1:quant,
                     cor=colors,#topo.colors
                     values=round(codonUsage$count/1000,0),
                     #cor=colorRampPalette(c("blue","green", "red"))( length(val) ),
                     stringsAsFactors = F)
  mid<-round(length(mydf$values)/2,0)
  ticks<-data.frame(pos=c(1,mid,length(mydf$values)),
                    val=c(mydf$values[1],
                          mydf$values[mid],
                          mydf$values[length(mydf$values)]))
  #legend bar code
  legenda<-
    ggplot(mydf) +
    ggtitle("",subtitle = "Codon Usage\n(x1000)")+
    theme(title = element_text(size = 6),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size = 6),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())+
    geom_tile(aes(x = id, 
                  y=usage,
                  fill=cor)) +
    scale_fill_identity()

  library(grid)
  library(gridExtra)

  lay <- rbind(c(1,1,1,1,1,1,1,1,1,1,NA),
               c(1,1,1,1,1,1,1,1,1,1,NA),
               #c(1,1,1,1,1,1,1,1,1,1,3),
               c(1,1,1,1,1,1,1,1,1,1,2),
               c(1,1,1,1,1,1,1,1,1,1,2),
               c(1,1,1,1,1,1,1,1,1,1,NA),
               c(1,1,1,1,1,1,1,1,1,1,NA))
  g<-grid.arrange(p,legenda, layout_matrix = lay)
  
  if(save){
    ggsave(filename = paste0("../figures/enrich",type,#rankStd,
                             population,
                             #Dw,
                             ".pdf"), 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 0.6, 
           width = 16.610, height = 7.440, units = "in",
           dpi = 300)
  }
  
  return(g)
}


correlation <- function(population="HighLow",
                      save = F,
                      normalize = F,
                      normBy = "count",
                      workdir ,
                      type = "sum",
                      TsTv="all"){

  if(!type %in% c("mean","sum")){
    stop('Parameter "type" must be "mean" or "sum"' )
  }
  if(!TsTv %in% c("Ts","Tv","all")){
    stop('Parameter "TsTv" must be "Ts", "Tv" or "all"' )
  }
  if(normalize){
    if(!normBy%in%c("CUB","count","mean")){
      stop('Normalizations can be done only by Codon Usage Bias ("CUB"),mean or mutation count ("count")')
    }
  }
  
  vTsTv<-TsTv
  setwd(workdir)
  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  
  TsTv<-read.csv("AuxFiles/TsTv.csv",
                 header = T,
                 stringsAsFactors = F)
  
  mutAllel<-read.csv(paste0("base/a",population,"MutAllDW.csv"),
                     header = T,
                     stringsAsFactors = F)

  #remove stop codons
  mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"|mutAllel$ref == "TAG" ),]

  #só syn
  mutAllelSyn<-mutAllel[mutAllel$Annotation == "synonymous",]

  #add fator based on codon usage
  mutAllelSyn<-merge(mutAllelSyn,codonUsage[,c(1:3)],by=c("ref"))
  colnames(mutAllelSyn)[137]<-"factor"
  
  #normalize values
  if(normalize){
    if(normBy =="count"){
      #remove CUB factor
      mutAllelSyn[,137]<-NULL
      #create mutation count
      qtCount<-mutAllelSyn[,c(1,7:128)]
      qtCount[is.na(qtCount)]<-0
      qtCount$sum<-apply(qtCount[,-1],
                         MARGIN = 1,
                         FUN = sum)
      qtCount<-qtCount[qtCount$sum !=0,c("ref","sum")]
      
      qtCount<-qtCount%>%
        group_by(ref)%>%
        summarize(factor=n())
      
      #use mutation coust as normalization factor
      mutAllelSyn<-merge(mutAllelSyn,qtCount,by=c("ref"))
      tmp<-mutAllelSyn[,7:128]/mutAllelSyn$factor
      tmp[is.na(tmp)]<- 0
      mutAllelNorm<-cbind(mutAllelSyn[,1:6],tmp,mutAllelSyn[129:137])
    }else if(normBy =="mean"){
      mutAllelNorm<-mutAllelSyn
      mutAllelNorm$order<-as.numeric(rownames(mutAllelNorm))
      mutAllelNorm[is.na(mutAllelNorm)]<-0
      mutAllelNorm<-mutAllelNorm[order(mutAllelNorm$order),]
      #mutAllelNorm<-mutTmp
      
      #colnames(qt)<-c("ref","y","factor","ordem")
      #colnames(qtFactors)<-c("ref","count","CUB","ordem")
      col="X1000"
      for( col in colnames(mutAllelSyn[7:128])){
        tmpTime<-mutAllelNorm[,c("ref","order",col)]
        tmpCount<-mutAllelNorm[mutAllelNorm[,col]>0,c("ref",col)]
        tmpCount<-tmpCount%>%
          group_by(ref)%>%
          summarize(count=n())
        colnames(tmpCount)<-c("ref","count")
        tmpTime<-merge(tmpTime,tmpCount, 
                         by= "ref",
                         all.x = T)
        tmpTime$count[is.na(tmpTime$count)]<-1
        tmpTime$norm<-tmpTime[,col]/tmpTime$count
        tmpTime<-tmpTime[order(tmpTime$order),]
        mutAllelNorm[,col]<-tmpTime$norm
      }
      tmp<-tmpTime
      rm(tmpTime,tmpCount)
      
    }else if(normBy =="CUB"){
      tmp<-mutAllelSyn[,7:128]/mutAllelSyn$factor
      tmp[is.na(tmp)]<- 0
      mutAllelNorm<-cbind(mutAllelSyn[,1:6],tmp,mutAllelSyn[129:137])
    }
  }else{
    tmp<-mutAllelSyn[,7:128]
    tmp[is.na(tmp)]<- 0
    mutAllelNorm<-cbind(mutAllelSyn[,1:6],tmp,mutAllelSyn[129:137])
  }
  rm(tmp)
  
  
  mutAllelNorm$maxFreq<-apply(X = mutAllelNorm[,7:128],
                              MARGIN = 1,
                              max)
  if(vTsTv != "all"){
    mutAllelNorm<-mutAllelNorm[mutAllelNorm$TsTv == vTsTv,]
  }
  library(dplyr)
  library(ggplot2)
  if(type == "mean"){
    sumario<-mutAllelNorm %>% 
      group_by(ref,mut,count) %>% 
      summarise(freq = mean(maxFreq))
  }else{
    sumario<-mutAllelNorm %>% 
      group_by(ref,mut,count) %>% 
      summarise(freq = sum(maxFreq))
  }
  lm<-lm(sumario$freq~sumario$count)
  corr<-cor.test(sumario$freq,sumario$count,
                 method = "pearson")
  titulo<- paste(population,vTsTv)
  g<-ggplot()+theme_bw()+
    ggtitle(titulo) +
    xlab("Codon Usage")+
    # ylab(paste("Mutation Frequency (",
    #            ifelse(normalize,"Normalized)","Absolute)")))+
    ylab("Frequency")+
    geom_point(data = sumario,
               aes(x=count,y=freq), 
               col="red", 
               pch=1)+
    geom_abline(aes(linetype = "line"),slope = lm$coefficients[2], 
                intercept = lm$coefficients[1],
                lty=2,
                color="blue", show.legend = T)+
    geom_text(aes(x=max(sumario$count)*0.9, y=max(sumario$freq)*0.9,
                  label=paste("Correlation:",
                              round(corr$estimate,3),
                              "\npvalue:",
                              format(corr$p.value,digits = 3, scientific=T))))+
    scale_colour_manual(values="blue")
  
  print(g)
  if(save){
    ggsave(filename = paste0("../figures/Corr",population,vTsTv,
                             ifelse(normalize,"Norm","NotNorm"),".pdf"), 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 0.6, 
           width = 16.610, height = 7.440, units = "in",
           dpi = 300)
  }
  return(g)
  
}

freqPerDw <- function(population="High",
                      save = F,
                      workdir ,
                      type = "sum",
                      Dw = "Tai",
                      TsTv="all",
                      top=86,
                      normBy = "count"){
  if(!type %in% c("mean","sum")){
    stop('Parameter "type" must be "mean" or "sum"' )
  }
  if(!TsTv %in% c("Ts","Tv","all")){
    stop('Parameter "TsTv" must be "Ts", "Tv" or "all"' )
  }
  library(ggplot2)
  library(grid)
  library(gridExtra)
  vTsTv<-TsTv
  setwd(workdir)
  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  
  TsTv<-read.csv("AuxFiles/TsTv.csv",
                 header = T,
                 stringsAsFactors = F)
  
  mutAllel<-read.csv(paste0("base/a",population,"MutAllDW.csv"),
                     header = T,
                     stringsAsFactors = F)
  genes<-read.csv("AuxFiles/genes.txt",
                  header = F,
                  stringsAsFactors = F)
  genes$len1<-apply(X = genes,
                    MARGIN = 1,
                    FUN = function(x){
                      return(nchar(x[5]))
                    })
  
  if(Dw == "Tai"){
    deltaW<-read.csv("./AuxFiles/cDeltaWTAI.csv",
                     header = T,
                     stringsAsFactors = F)
  }else if(Dw == "Cai"){
    deltaW<-readDeltaW(top,workdir)
  }
  
  
  # genes$len2<-genes$V4-genes$V3+1
  # # if(nrow(genes[genes$len1!=genes$len2,])!=0){
  # #   cat("Deu ruim")
  # # }
  # genes$len2<-NULL
  
  colnames(genes)<-c("Gene","strand","ini","fim","seq","len")
  #remove stop codons
  mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"|mutAllel$ref == "TAG" ),]
  duplos<-mutAllel[duplicated(mutAllel[,c(3,5)]),c(3,5)]
  
  #só syn
  mutAllelSyn<-mutAllel[mutAllel$Annotation == "synonymous",]
  teste<-mutAllel[mutAllel$Annotation == "missense" | mutAllel$Annotation == "nonsense" ,]
  
  duplos<-mutAllelSyn[duplicated(mutAllelSyn[,c(3,5)]),c(3,5)]
  
  #atualize dw values 
  mutAllelSyn<-merge(mutAllel,deltaW,
                  by=c("ref","mut"))
  mutAllelSyn$dw.x<-NULL
  colnames(mutAllelSyn)[135]<-"dw"
  
  #hist(mutAllelSyn$dw[mutAllelSyn$Position %in% duplos$Position&
  #                      mutAllelSyn$Allele %in% duplos$Allele])
  #hist(mutAllelSyn$dw)
  #acrescenta fator
  mutAllelSyn<-merge(mutAllelSyn,codonUsage[,c(1,3)],by=c("ref"))
  colnames(mutAllelSyn)[136]<-"factor"
  
  #normaliza valores
  tmp<-mutAllelSyn[,c(1,7:128)]
  tmp[is.na(tmp)]<- 0
  tmp<-normalization(data = tmp,
                     normBy = normBy,
                     melted = F,
                     workdir = workdir )
  mutAllelNorm<-cbind(mutAllelSyn[,1:6],tmp[,-1],mutAllelSyn[129:136])
  
  rm(tmp)
  
  
  mutAllelNorm$maxFreq<-apply(X = mutAllelNorm[,7:128],
                              MARGIN = 1,
                              max)
  if(vTsTv != "all"){
    mutAllelNorm<-mutAllelNorm[mutAllelNorm$TsTv == vTsTv,]
  }
  library(dplyr)
  library(ggplot2)
  library(grDevices)
  library(gridGraphics)
  if(type == "mean"){
    sumario<-mutAllelNorm %>% 
      group_by(ref,mut,dw) %>% 
      summarise(freq = mean(maxFreq))
  }else{
    sumario<-mutAllelNorm %>% 
      group_by(ref,mut,dw) %>% 
      summarise(freq = sum(maxFreq))
  }
  zeros<-deltaW[!paste(deltaW$ref,deltaW$mut)%in%paste(sumario$ref,sumario$mut),]
  #zeros<-deltaW$dw[!round(deltaW$dw,5)%in%round(sumario$dw,5)]
  #zeros<-data.frame(dw=zeros,freq=rep(0,length(zeros)))
  zeros$freq<-0
  
  #sumario<-rbind(sumario,zeros)
  regression<-lm(sumario$freq~sumario$dw)
  sumario<-merge(sumario,
                 TsTv[,c("ref","mut","type")], 
                 by = c("ref","mut"))
  zeros<-merge(zeros,
               TsTv[,c("ref","mut","type")], 
               by = c("ref","mut"))
  if(vTsTv != "all"){
    zeros<-zeros[zeros$type == vTsTv,]
  }
  if(population == "MutT"){
    mutT<-read.csv("AuxFiles/dMutTmutations.csv",
                   header = T,
                   stringsAsFactors = F)
    mutT<-merge(mutT,deltaW[,c(1,2,3)], by=c("ref","mut"))
    
    sumario$MT<-0
    sumario$MT[paste(sumario$ref,sumario$mut)%in% paste(mutT$ref,mutT$mut) ]<-1
    
    g<-ggplot()+theme_bw()+
      ggtitle(population)+
      xlab(bquote(Delta~"w"~.(toupper(Dw))))+
      ylab(paste("Mutation Frequency (",type,"Normalized)"))+
      geom_point(data = sumario[sumario$MT==0 & sumario$type == "Tv",],
                 aes(dw,freq),col="red",pch=20)+
      geom_point(data = sumario[sumario$MT==0 & sumario$type == "Ts",],
                 aes(dw,freq),col="darkmagenta",pch=18)+
      geom_point(data = sumario[sumario$MT==1,],
                 aes(dw,freq),col="darkgoldenrod2",pch=8)+
      geom_point(data = zeros,
                 aes(dw,freq),col="green",pch=2)+
      geom_abline(slope = regression$coefficients[2],
                  intercept = regression$coefficients[1],
                  lty=2,col="blue")
    plot.new()
    pWidth = 2
    pHeight = 1
    plot.window(c(0,pWidth),
                c(0,pHeight))
    
    leg <- legendGrob(c(" Tv", " Ts"," Zeros","MutT"), 
                      pch=c(20,18,2,8),
                      gp=gpar(col = c("red", "darkmagenta","green","darkgoldenrod2"),
                              fill = "gray",
                              cex=0.7))    
    
    
  }else{
    g<-ggplot()+theme_bw()+
      ggtitle(population)+
      xlab(bquote(Delta~"w"~.(toupper(Dw))))+
      ylab(paste("Mutation Frequency (",type,"Normalized)"))+
      geom_point(data = sumario[sumario$type == "Tv",],
                 aes(dw,freq),col="red",pch=20)+
      geom_point(data = sumario[sumario$type == "Ts",],
                 aes(dw,freq),col="darkmagenta",pch=18)+
      geom_point(data = zeros,
                 aes(dw,freq),col="green",pch=2)+
      geom_abline(slope = regression$coefficients[2],
                  intercept = regression$coefficients[1],
                  lty=2,col="blue")
    plot.new()
    pWidth = 2
    pHeight = 1
    plot.window(c(0,pWidth),
                c(0,pHeight))
    
    leg <- legendGrob(c(" Tv", " Ts"," Zeros"), 
                      pch=c(20,18,2),
                      gp=gpar(col = c("red", "darkmagenta","green"),
                               fill = "gray",
                              cex=0.7))    
  }
  g<-g+geom_text(aes(y=(max(sumario$freq)*2/3), 
                     x=(min(deltaW$dw)), 
                     label = paste(nrow(zeros),"zeros"),
                     hjust = "left"))
  
  lay <- rbind(c(1,1,1,1,1,1,1,1,1,1,NA),
               c(1,1,1,1,1,1,1,1,1,1,NA),
               #c(1,1,1,1,1,1,1,1,1,1,3),
               c(1,1,1,1,1,1,1,1,1,1,2),
               c(1,1,1,1,1,1,1,1,1,1,2),
               c(1,1,1,1,1,1,1,1,1,1,NA),
               c(1,1,1,1,1,1,1,1,1,1,NA))
  g2<-grid.arrange(g,leg, layout_matrix = lay)
  
  if(save){
    ggsave(filename = paste0("../figures/MutPerDw",population,".pdf"), 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 0.6, 
           width = 16.610, height = 7.440, units = "in",
           dpi = 300)
  }
  return(g2)
}

genesExpression<-function(workdir){
  
  
  setwd(workdir)
  setwd("./dataExpresion/LeGac2012/")
  
  
  #If is not avaliable, install ecoli2.db
  if(!requireNamespace("BiocManager", quietly = TRUE))
    suppressMessages(install.packages("BiocManager"))
  if(!requireNamespace("ecoli2.db",quietly = TRUE)){
    suppressMessages(BiocManager::install("ecoli2.db", version = "3.8"))
  }
  
  suppressMessages(library("org.EcK12.eg.db",quietly = T))
  suppressMessages(library("ecoli2.db",quietly = T))
  suppressMessages(library("ggplot2",quietly = T))
  
  #get genes that codify proteins
  rs<-ecoli2REFSEQ
  rms<- mappedkeys(rs)
  justProteins<-as.data.frame(rs[rms])
  #justProteins<-unique(justProteins$probe_id[substr(justProteins$accession,1,2)=="NP"])
  
  
  ## Bimap interface:
  x <- ecoli2SYMBOL
  # Get the probe identifiers that are mapped to a gene symbol
  mapped_probes <- mappedkeys(x)
  # Convert to a data frame
  symbol2Probe <- as.data.frame(x[mapped_probes])
  
  #remove non protein genes
  symbol2Probe<-symbol2Probe[symbol2Probe$probe_id%in%justProteins$probe_id,]
  
  
  expression<-read.csv("GSE30639SeriesMatrix.csv")
  colnames(expression)<-c("probe_id",colnames(expression)[2:6])
  
  expression<-merge(expression,symbol2Probe, by="probe_id")
  expression<-aggregate(. ~symbol, expression[,c(2:7)], mean)
  justProteins2<-unique(expression$symbol)
  
  #PCA test
  PCAData<-t(expression[,2:6])
  # calcular a PCA
  pca <- prcomp(PCAData,scale = T)
  
  # Plotar
  pcs<-data.frame(pca$x)
  shapes<-as.factor(rownames(pcs))
  color<-as.factor(rownames(pcs))
  names<-as.character(rownames(pcs))
  g<-ggplot(data = pcs[,1:2],
            aes(PC1,PC2,
                label=names,
                #shape=shapes,
                color = color))+
    geom_point()#+
  g
  
  samplesOk<-c("GSM759848","GSM759849","GSM759851")
  
  expression$median<-(apply(expression[,samplesOk],MARGIN = 1,median))
  expression<-expression[order(expression$median,decreasing = T),c(1,7)] 
  library(stringr)
  expression$symbol<-str_trim( expression$symbol)
  return(expression)
}

findFirstTop<-function(population="HighLow",
                       justSyn = T,
                       workdir ){
  setwd(workdir)
  source("../bin/allFunctions.R")
  mutAllelAll<-read.csv(paste0("base/a",population,"MutAllDW.csv"),
                        header = T,
                        stringsAsFactors = F)
  #só syn
  if(justSyn){
    mutAllelAll<-mutAllelAll[mutAllelAll$Annotation == "synonymous",]
  }
  expression<-genesExpression(workdir)
  mutGenes<-unique(mutAllelAll$Gene)
  i=788
  for(i in 1:nrow(expression)){
    #cat(expression$symbol[i],".")
    if(expression$symbol[i]%in%mutGenes){
      cat("Firt mutated gene is the top", i, "in expression",expression$symbol[i],"\n" )
      break()
    }
  }  
}


listZeros <- function(save = F,
                      quant=5,
                      Dw = "Cai",
                      top=86,
                      withLow = F,
                      workdir ){
  
  setwd(workdir)
  if(Dw == "Tai"){
    deltaW<-read.csv("./AuxFiles/cDeltaWTAI.csv",
                     header = T,
                     stringsAsFactors = F)
  }else if(Dw == "Cai"){
    deltaW<-readDeltaW(top = top,
                       workdir = workdir)
  }else{
    deltaW<-calcDwOld(workdir)
  }
  
  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  if(exists("zeros")){
    rm(zeros)
  }
  
  population="High"
  population="MutT"
  population="HighLow"
  population="Low"
  if(withLow){
    populations<-c("High","MutT","HighLow","Low")
  }else{
    populations<-c("High","MutT","HighLow")
  }
  for(population in populations){
    
    mutAllel<-read.csv(paste0("base/a",population,"MutAllDW.csv"),
                       header = T,
                       stringsAsFactors = F)
    #remove stop codons
    mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"|mutAllel$ref == "TAG" ),]
    duplos<-mutAllel[duplicated(mutAllel[,c(3,5)]),c(3,5)]
    
    #só syn
    mutAllelSyn<-mutAllel[mutAllel$Annotation == "synonymous",]
    teste<-mutAllel[mutAllel$Annotation == "missense" | mutAllel$Annotation == "nonsense" ,]
    #atualize dw values 
    mutAllelSyn<-merge(mutAllelSyn,deltaW,
                    by=c("ref","mut"))
    mutAllelSyn$dw.x<-NULL
    colnames(mutAllelSyn)[135]<-"dw"
    
    duplos<-mutAllelSyn[duplicated(mutAllelSyn[,c(3,5)]),c(3,5)]
    #acrescenta fator
    mutAllelSyn<-merge(mutAllelSyn,codonUsage[,c(1,3)],by=c("ref"))
    colnames(mutAllelSyn)[136]<-"factor"
    
    #normaliza valores
    tmp<-mutAllelSyn[,7:128]/mutAllelSyn$factor
    tmp[is.na(tmp)]<- 0
    mutAllelNorm<-cbind(mutAllelSyn[,1:6],tmp,mutAllelSyn[129:136])
    
    rm(tmp)
    
    
    mutAllelNorm$maxFreq<-apply(X = mutAllelNorm[,7:128],
                                MARGIN = 1,
                                max)
    library(dplyr)
    sumario<-mutAllelNorm %>% 
      group_by(ref,mut,dw) %>% 
      summarise(freq = mean(maxFreq))
    if(!exists("zeros")){
      zeros<-deltaW[!paste(deltaW$ref,deltaW$mut)%in%paste(sumario$ref,sumario$mut),]
      #zeros<-deltaW[!round(deltaW$dw,5)%in%round(sumario$dw,5),]
      #mark on zeros from HighLow
      zeros[[population]]<-"X"
    }else{
      tmp<-deltaW[!paste(deltaW$ref,deltaW$mut)%in%paste(sumario$ref,sumario$mut),
                  c("ref","mut","dw")]
      tmp[[population]]<-"X"
      zeros<-merge(zeros,tmp,
                   by=c("ref","mut","dw"),
                   all=T)
    }
  }
  MutT<-read.csv(file = "AuxFiles/dMutTmutations.csv",
                 header = T,
                 stringsAsFactors = F)
  MutT$isMutT<-"*"
  zeros<-merge(zeros,MutT,
               by=c("ref","mut"),
               all.x=T)
  zeros[is.na(zeros)]<-""
  quantil<-quantile(deltaW$dw,c(seq(from = 0,to = 1,by = 1/quant)),type = 1)
  vdw<-data.frame(min=round(quantil[seq(1,quant)],2),max=round(quantil[seq(2,quant+1)],2))
  quantil[1]<-quantil[1]-0.1
  
  for(i in 1:(length(quantil)-1)){
    if(i == 1){
      zeros$level[zeros$dw == quantil[1]]<-1
    }
    zeros$level[zeros$dw> quantil[i] & zeros$dw<=quantil[i+1]]<-i
  }
  #garante a mesma ordem me mutAllel e top100Base
  zeros<-zeros[order(zeros$dw),]
  zeros$dw<-round(zeros$dw,5)
  # colnames(zeros)<-c("ref","mut","dw", 
  #                    populations,
  #                    "isMutT","level")
  msg<-"\n\nLevels: \n"
  for (qn in 1:(length(quantil)-1)) {
    if(qn==1){
      limDown<-quantil[qn]-0.1
    }else{
      limDown<-quantil[qn]
    }
    qtd<-nrow(deltaW[deltaW$dw>limDown &
                       deltaW$dw<=quantil[qn+1],])
    msg<-paste0(msg,"\t",
                qn,"-\t",
                round(quantil[qn],2)," to ",
                round(quantil[qn+1],2),"\t- ",qtd," values\n")
  }
  print(zeros,quote = F)
  cat(msg)
  if(save){
    write.csv(x=zeros,
              file = paste0(workdir,
                            "./AuxFiles/zeros.csv"),
              row.names = F,
              quote = F)
    write.table(x=zeros,
                file = paste0(workdir,
                              "../figures/zeros.txt"),
                sep = "\t",
                row.names = F,
                col.names = T,
                quote = F)
    cat(msg,
        file = paste0(workdir,
                      "../figures/zeros.txt"),
        append = T)
    
  }
  return(zeros)
}

convertToComplement<-function(x){
  bases=c("A","C","G","T")
  xx<-unlist(strsplit(toupper(x),NULL))
  xxx<-paste(unlist(lapply(xx,function(bbb){
    if(bbb=="A") compString<-"T"
    if(bbb=="C") compString<-"G"
    if(bbb=="G") compString<-"C"
    if(bbb=="T") compString<-"A"
    if(!bbb %in% bases) compString<-"N"
    return(compString)
  })),collapse="")
  return(sapply(lapply(strsplit(xxx, NULL), rev), paste, collapse=""))
}

calcDWTAI<-function(type="Min",
                    workdir){
  setwd(workdir)
  tRNA<-read.table(paste0("TAIFiles/output_wi_file",type,".txt"),
                 header = T,
                 sep="\t",
                 stringsAsFactors = F)
  colnames(tRNA)<-c("ref","w")
  deltaW<-read.csv("AuxFiles/cDeltaW86.csv",
                   header = T,
                   stringsAsFactors = F)
  deltaW$dw<-0
  
  deltaW<-merge(deltaW,
                tRNA[,c("ref","w")],
                by ="ref")
  colnames(deltaW)<- c("ref","mut","dw","wRef")
  deltaW<-merge(deltaW,
                tRNA[,c("ref","w")],
                by.x ="mut",
                by.y="ref")
  colnames(deltaW)<- c("mut","ref","dw","wRef", "wMut")
  deltaW$dw<-log10(deltaW$wMut/deltaW$wRef)
  deltaW<-deltaW[,c("ref","mut","dw")]
  write.csv(deltaW,
            file = "AuxFiles/cDeltaWTAI.csv",
            quote = F,
            row.names = F)
}

corrTaiCai<- function(top=86,
                      save = F,
                      figName,
                      workdir){
  setwd(workdir)
  deltaWCAI<-readDeltaW(top,workdir)
  deltaWTAI<-read.csv("AuxFiles/cDeltaWTAI.csv",
                      header = T,
                      stringsAsFactors = F)
  library(ggplot2)
  deltas<-merge(deltaWCAI,deltaWTAI,
                by=c("ref","mut"))
  colnames(deltas)<-c("ref","mut","CAI","TAI")
  lm<-lm(deltas$TAI~deltas$CAI)
  corr<-cor.test(deltas$TAI,deltas$CAI)
  maxX<-max(deltas$CAI)
  maxY<-max(deltas$TAI)
  
  tmp1<-data.frame(table(deltas$TAI))
  tmp2<-data.frame(table(deltas$CAI))
  cat("Repeated dwTAI:\n")
  print(tmp1[tmp1$Freq!=1,])
  cat("Repeated dwCAI:\n")
  print(tmp2[tmp2$Freq!=1,])
  cat("Distinct dwTAI:", length(unique(deltas$TAI)),
      "\nDistinct dwCAI:", length(unique(deltas$CAI)),"\n")
  cat("Max dwTAI:", min(deltas$TAI),"to",max(deltas$TAI),
      "\nMax dwCAI:", min(deltas$CAI),"to",max(deltas$CAI),"\n")
  
  library(ggplot2)
  g<-ggplot()+theme_bw()+
    xlab(bquote(Delta~"w CAI"))+
    ylab(bquote(Delta~"w TAI"))+
    geom_point(data = deltas,
               aes(y=TAI,x=CAI), 
               col="red", 
               pch=1)+
    geom_abline(slope = lm$coefficients[2], 
                intercept = lm$coefficients[1],
                lty=2,
                col="blue")+
    geom_text(aes(x=-maxX*0.75, y=maxY*0.75,
                  label=paste("Correlation:",
                              round(corr$estimate,3),
                              "\npvalue:",
                              format(corr$p.value,digits = 3, scientific=T))))+
    scale_x_continuous(breaks = seq(-2,2,0.5),
                       minor_breaks = seq(-2,2,0.25))
  print(g)
  if(save){
    ggsave(filename = paste0("../figures/",figName,".pdf"), 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 1,
           width = 7, height = 4.75, units = "in",
           dpi = 300)
  }
}

get.ws <- function(tRNA,  # tRNA gene copy number
                   s = c(0.0, 0.0, 0.0, 
                         0.0, 0.41, 0.28, 
                         0.9999, 0.68, 0.89),     # selective constraints
                   sking) # super kingdom: 0-eukaryota, 1-prokaryota
{
  #default order of codons
  codonOrder<-c("TTT","TTC","TTA","TTG","TCT","TCC",
                "TCA","TCG","TAT","TAC","TAA","TAG",
                "TGT","TGC","TGA","TGG","CTT","CTC",
                "CTA","CTG","CCT","CCC","CCA","CCG",
                "CAT","CAC","CAA","CAG","CGT","CGC",
                "CGA","CGG","ATT","ATC","ATA","ATG",
                "ACT","ACC","ACA","ACG","AAT","AAC",
                "AAA","AAG","AGT","AGC","AGA","AGG",
                "GTT","GTC","GTA","GTG","GCT","GCC",
                "GCA","GCG","GAT","GAC","GAA","GAG",
                "GGT","GGC","GGA","GGG")
  #reorder the dataframe
  tRNA$order <- factor(tRNA$AntiCodonsList, 
                       levels = codonOrder)
  tRNA<-tRNA[order(tRNA$order),1:2]
  
  p = 1 - s
  
  # initialise w vector
  W = NULL  # don't confuse w (lowercase) and W (uppercase)
  
  # obtain absolute adaptiveness values (Ws)
  for (i in seq(1, 61, by=4))
    W = c(W,
          p[1]*tRNA$tGCN[i]   + p[5]*tRNA$tGCN[i+1],     # INN -> NNT, NNC, NNA
          p[2]*tRNA$tGCN[i+1] + p[6]*tRNA$tGCN[i],       # GNN -> NNT, NNC
          p[3]*tRNA$tGCN[i+2] + p[7]*tRNA$tGCN[i],       # TNN -> NNA, NNG
          p[4]*tRNA$tGCN[i+3] + p[8]*tRNA$tGCN[i+2])     # CNN -> NNG
  
  # check methionine
  W[36] = p[4]*tRNA$tGCN[36]
  
  # if bacteria, modify isoleucine ATA codon
  if(sking == 1) W[35] = p[9]
  
  # get ws
  w = W/max(W)
  
  if(sum(w == 0) > 0) {
    ws <- w[w != 0] # zero-less ws
    gm <- exp(sum(log(ws))/length(ws)) # geometric mean
    w[w == 0] = gm # substitute 0-ws by gm
  }
  
  tRNA$w<-w
  
  # get rid of stop codons (11, 12, 15) and methionine (36)
  tRNA <- tRNA[-c(11,12,15,36),]
  colnames(tRNA)<-c("codons", "tGCN", "w")
  return(tRNA)
}

totalize_tRNA<-function(workdir){
  setwd(workdir)
  setwd("./TAIFiles/")
  codons<-data.frame(codon=c("TTT","TTC","TTA","TTG","TCT","TCC",
                             "TCA","TCG","TAT","TAC","TAA","TAG",
                             "TGT","TGC","TGA","TGG","CTT","CTC",
                             "CTA","CTG","CCT","CCC","CCA","CCG",
                             "CAT","CAC","CAA","CAG","CGT","CGC",
                             "CGA","CGG","ATT","ATC","ATA","ATG",
                             "ACT","ACC","ACA","ACG","AAT","AAC",
                             "AAA","AAG","AGT","AGC","AGA","AGG",
                             "GTT","GTC","GTA","GTG","GCT","GCC",
                             "GCA","GCG","GAT","GAC","GAA","GAG",
                             "GGT","GGC","GGA","GGG"),
                     freq= rep(0,64) )
  
  codons[c(11, 12, 15, 36),]
  
  tRNA<-read.csv(file="data_files/TAIFiles/REL606tRNAPre.txt",
                 header =  F,
                 stringsAsFactors = F)
  tRNA<-tRNA[tRNA$V2!="NNN",]
  table(tRNA[,c("V1")])
  tmp<-as.data.frame(table(tRNA[,c("V2")]),stringsAsFactors = F)
  colnames(tmp)<-c("codon", "freq")
  for(i in 1:nrow(tmp)){
    codons$freq[codons$codon == tmp$codon[i]]<-tmp$freq[i]
  }
  #codons<-merge(codonOrder,codons, by="codon", all.x = T, sort =F)
  colnames(codons)<-c("AntiCodonsList","tGCN")
  write.table(codons,file = "REL606tRNA.txt",
              sep = "\t",
              row.names = F,
              col.names = T,
              quote = F)
}




createMonteCarloDF<-function(workdir ){
  setwd(workdir)
  
  #load MySQL Lib
  library('RMySQL')
  
  try(dbDisconnect(con),silent = T)
  #DB connection
  con <- dbConnect(MySQL(),
                   user = 'clovis',
                   password = '230209cf',
                   host = 'localhost',
                   dbname='ecoli60k')
  
  #Extract Monte Carlo mutations means
  rs = dbSendQuery(con,
                   paste('select ref, mut ,sum(freq) as freq, sum(freq)/20000 as Avg, std(freq) as Sd
                        from MCMutationsTv 
                        group by ref, mut;'))
  #Load Monte Carlo values
  MCMutations = fetch(rs, n=-1)

  rs = dbSendQuery(con,
                   paste('select ref,mut, count(*)/20000 as pEmp
                          from(select simulation, run,ref, mut,count(*)
                                from MCMutationsTv
                                group by simulation, run,ref, mut) as group1
                          group by ref,mut;'))
  #Load Monte Carlo values
  pEmp = fetch(rs, n=-1)
  
  
    
  #disconect
  dbDisconnect(con)
  
  write.csv(MCMutations, 
            file = "base/MonteCarloMutationsTv.csv",
            row.names = F )
  write.csv(pEmp, 
            file = "base/MonteCarloPEmp.csv",
            row.names = F )
}

corrMCxDw<- function(Dw = "Tai",
                     type="box",
                     quant=5,
                     save=F,
                     workdir){
  library(ggplot2)
  if(!type %in% c("box","point")){
    stop('Parameter "type" must be "point" or "box"')
    
  }
  setwd(workdir)
  if(Dw == "Tai"){
    deltaW<-read.csv("./AuxFiles/cDeltaWTAI.csv",
                     header = T,
                     stringsAsFactors = F)
  }else if(Dw == "Cai"){
    deltaW<-read.csv("./AuxFiles/cDeltaW86.csv",
                     header = T,
                     stringsAsFactors = F)
  }else{
    deltaW<-calcDwOld(workdir)
  }
  MCMut<-read.csv("base/MonteCarloMutationsTv.csv",
                       header = T,
                       stringsAsFactors = F)
  MCMut<-merge(MCMut[,c("ref","mut","freq")],
               deltaW,
               by=c("ref","mut"))
  if(type == "point"){
    lm<-lm(MCMut$freq~MCMut$dw)
    corr<-cor.test(MCMut$freq,MCMut$dw)
    maxX<-max(MCMut$dw)
    maxY<-max(MCMut$freq)
    library(ggplot2)
    g<-ggplot()+theme_bw()+
      xlab(bquote(Delta~"w"~.(Dw)))+
      ylab("Mutation frequency")+
      geom_point(data = MCMut,
                 aes(x=dw,y=freq), 
                 col="red", 
                 pch=1)+
      geom_abline(slope = lm$coefficients[2], 
                  intercept = lm$coefficients[1],
                  lty=2,
                  col="blue")+
      geom_text(aes(x=-maxX*0.75, y=maxY*0.75,
                    label=paste("Correlation:",
                                round(corr$estimate,3),
                                "\npvalue:",
                                format(corr$p.value,digits = 3, scientific=T))))
  }else{
    quantil<-quantile(deltaW$dw,c(seq(from = 0,to = 1,by = 1/quant)),type = 1)
    cat("Levels: \n")
    for (qn in 1:(length(quantil)-1)) {
      if(qn==1){
        limDown<-quantil[qn]-0.1
      }else{
        limDown<-quantil[qn]
      }
      qtd<-nrow(deltaW[deltaW$dw>limDown &
                         deltaW$dw<=quantil[qn+1],])
      msg<-paste0("\t",
                  qn,"-\t",
                  round(quantil[qn],2)," to ",
                  round(quantil[qn+1],2),"\t- ",qtd," values\n")
      cat(msg)
    }
    vdw<-data.frame(min=round(quantil[seq(1,quant)],2),max=round(quantil[seq(2,quant+1)],2))
    quantil[1]<-quantil[1]-0.1
    
    for(i in 1:(length(quantil)-1)){
      if(i == 1){
        MCMut$faixa[MCMut$dw == quantil[1]]<-1
      }
      MCMut$faixa[MCMut$dw> quantil[i] & MCMut$dw<=quantil[i+1]]<-i
    }
    MCMut$Levels<-as.factor(MCMut$faixa)
    g<-ggplot()+theme_bw()+
      xlab(bquote(Delta~"w"~.(Dw)))+
      ylab("Mutation frequency")+
      geom_boxplot(data = MCMut,
                 aes(x=dw,y=freq, 
                     group=faixa,
                     fill=Levels))
      
  }
  print(g)
  if(save){
    ggsave(filename = paste0("../figures/corrMCarloDw",
                             Dw,type,".pdf"), 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 0.6,
           width = 16.610, height = 7.440, units = "in",
           dpi = 300)
  }
}

enrichZeros <- function(save = F,
                      quant=5,
                      population = "HighLow",
                      Dw = "Tai",
                      workdir ){
  setwd(workdir)
  if(Dw == "Tai"){
    deltaW<-read.csv("./AuxFiles/cDeltaWTAI.csv",
                     header = T,
                     stringsAsFactors = F)
  }else if(Dw == "Cai"){
    deltaW<-read.csv("./AuxFiles/cDeltaW86.csv",
                     header = T,
                     stringsAsFactors = F)
  }else{
    deltaW<-calcDwOld(workdir)
  }
  
  TsTv<-read.csv("AuxFiles/TsTv.csv",
                 header = T,
                 stringsAsFactors = F)
  deltaW<-merge(deltaW,TsTv[,c("ref","mut","type")],
                by=c("ref","mut"))
  
  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  if(!population%in%c("High","MutT","HighLow")){
    stop('Parameter "population" must be "High","MutT", or "HighLow"')
  }
  mutAllel<-read.csv(paste0("base/a",population,"MutAllDW.csv"),
                     header = T,
                     stringsAsFactors = F)
  #remove stop codons
  mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"|mutAllel$ref == "TAG" ),]
  duplos<-mutAllel[duplicated(mutAllel[,c(3,5)]),c(3,5)]
  
  #só syn
  mutAllelSyn<-mutAllel[mutAllel$Annotation == "synonymous",]

  #atualize dw values 
  mutAllelSyn<-merge(mutAllelSyn,deltaW[,c("ref","mut","dw")],
                     by=c("ref","mut"))
  mutAllelSyn$dw.x<-NULL
  colnames(mutAllelSyn)[135]<-"dw"
  
  duplos<-mutAllelSyn[duplicated(mutAllelSyn[,c(3,5)]),c(3,5)]
  #add fator
  mutAllelSyn<-merge(mutAllelSyn,codonUsage[,c(1,3)],by=c("ref"))
  colnames(mutAllelSyn)[136]<-"factor"
  
  #normaliza valores
  tmp<-mutAllelSyn[,7:128]/mutAllelSyn$factor
  tmp[is.na(tmp)]<- 0
  mutAllelNorm<-cbind(mutAllelSyn[,1:6],tmp,mutAllelSyn[129:136])
  
  rm(tmp)
  
  
  mutAllelNorm$maxFreq<-apply(X = mutAllelNorm[,7:128],
                              MARGIN = 1,
                              max)
  library(dplyr)
  sumario<-mutAllelNorm %>% 
    group_by(ref,mut,TsTv) %>% 
    summarise(freq = n())
  sumario<-merge(sumario,deltaW,
                 bu=c("ref","mut"))
  sumario$TsTv<-NULL
  zeros<-deltaW[!paste(deltaW$ref,deltaW$mut)%in%paste(sumario$ref,sumario$mut),]
  zeros$freq<-0
  #just transversions
  transv<-rbind(sumario[sumario$type == "Tv",
                        c("ref","mut","type","dw","freq")],
                zeros[, c("ref","mut","type","dw","freq")])
  
  #zeros<-deltaW[!round(deltaW$dw,5)%in%round(sumario$dw,5),]
  quantil<-quantile(deltaW$dw,c(seq(from = 0,to = 1,by = 1/quant)),type = 1)
  vdw<-data.frame(min=round(quantil[seq(1,quant)],2),max=round(quantil[seq(2,quant+1)],2))
  quantil[1]<-quantil[1]-0.1
  
  #load Monte Carlo data
  MCMut<-read.csv(file = "base/MonteCarloMutationsTv.csv",
                   stringsAsFactors = F, 
                   header = T)
  MCMut<-merge(MCMut[,c("ref","mut","freq")],
               deltaW,
               by=c("ref","mut"))
  #levels classification
  for(i in 1:(length(quantil)-1)){
    if(i == 1){
      transv$level[transv$dw == quantil[1]]<-1
      MCMut$level[MCMut$dw == quantil[1]]<-1
    }
    transv$level[transv$dw> quantil[i] & transv$dw<=quantil[i+1]]<-i
    MCMut$level[MCMut$dw> quantil[i] & MCMut$dw<=quantil[i+1]]<-i
  }
  #test all zeros
  transv<-transv[order(transv$dw),]
  colnames(transv)[5]<-"obs"
  transv<-merge(transv,MCMut[c("ref","mut","freq")],
                  by=c("ref","mut"),
                  all.y=T)
  transv$try<-round(sum(MCMut$freq),0)
  transv$prob<-transv$freq/transv$try
  library(stats)
  transv$pval1<-p.adjust(apply(transv[,c("obs", "try", "prob")],
                       MARGIN = 1,
                       FUN = function(x){
                         r<-binom.test(x[1],
                                    round(x[2],0),
                                    x[3],
                                    alternative = "g",
                                    conf.level = 0.9999)
                         # return(ifelse(r$p.value<=0.001,
                         #               r$p.value,
                         #               NA))}),
                         return(r$p.value)}),
                       method = "BH")
  #testing levels
  sumario<-transv %>% 
    group_by(level) %>% 
    summarise(obs = sum(obs),
              prob = sum(prob),
              try = mean(try))
  sumario$pval1<-p.adjust(apply(sumario[,c("obs", "try", "prob")],
                               MARGIN = 1,
                               FUN = function(x){
                                 r<-binom.test(x[1],
                                               round(x[2],0),
                                               x[3],
                                               alternative = "g",
                                               conf.level = 0.9999)
                                 # return(ifelse(r$p.value<=0.001,
                                 #               r$p.value,
                                 #               NA))}),
                                 return(r$p.value)}),
                         method = "BH")
  
  
  # allZeros<-transv[,c("ref","mut","dw")]
  # allZeros<-merge(allZeros,MCMut[c("ref","mut","freq")], 
  #                 by=c("ref","mut"),
  #                 all.y=T)
  # allZeros$obs<-0
  # allZeros$try<-sum(MCMut$freq)
  # allZeros$prob<-allZeros$freq/allZeros$try
  # library(stats)
  # allZeros$pval1<-p.adjust(apply(allZeros[,c("obs", "try", "prob")],
  #                      MARGIN = 1,
  #                      FUN = function(x){
  #                        r<-binom.test(x[1],
  #                                   round(x[2],0),
  #                                   x[3],
  #                                   alternative = "l",
  #                                   conf.level = 0.9999)
  #                        # return(ifelse(r$p.value<=0.001,
  #                        #               r$p.value,
  #                        #               NA))}),
  #                        return(r$p.value)}),                       
  #                      method = "BH")
  #                        
  # allZeros$pEmp<-allZeros$freq/20000
  #garante a mesma ordem me mutAllel e top100Base
  # zeros$dw<-round(zeros$dw,5)
  # colnames(zeros)<-c("ref","mut","dw", 
  #                    populations,
  #                    "isMutT","level")
  msg<-"\n\nLevels: \n"
  for (qn in 1:(length(quantil)-1)) {
    if(qn==1){
      limDown<-quantil[qn]-0.1
    }else{
      limDown<-quantil[qn]
    }
    qtd<-nrow(deltaW[deltaW$dw>limDown &
                       deltaW$dw<=quantil[qn+1],])
    msg<-paste0(msg,"\t",
                qn,"-\t",
                round(quantil[qn],2)," to ",
                round(quantil[qn+1],2),"\t- ",qtd," values\n")
  }
  print(zeros,quote = T)
  cat(msg)
  if(save){
    write.table(x=zeros,
                file = paste0(workdir,
                              "../figures/zeros.txt"),
                sep = "\t",
                row.names = F,
                col.names = T,
                quote = F)
    cat(msg,
        file = paste0(workdir,
                      "../figures/zeros.txt"),
        append = T)
    
  }
}

analiseZeros<-function(save = F,
                       population="High",
                       Dw = "Tai",
                       top=86,
                       workdir ){
  setwd(workdir)
  # if(!file.exists("AuxFiles/zeros.csv")){
  #   stop("File TsTv.csv do not exist.
  #        Did you run listZeros(save = T)?")
  # }
  if(!population%in%c("High","MutT","HighLow","Low","all")){
    stop('Parameter "population" must be "High","MutT", "Low" or "HighLow"')
  }
  # zeros<-read.csv("AuxFiles/zeros.csv",
  #                header = T,
  #                stringsAsFactors = F)
  zeros<-listZeros(save =F,
                   quant=5,
                   Dw=Dw,
                   top=top,
                   withLow = T,
                   workdir = workdir)
  TsTv<-read.csv("AuxFiles/TsTv.csv",
                 header = T,
                 stringsAsFactors = F)
  if(population == "all"){
    populations<-c("High","MutT","HighLow","Low")
  }else{
    populations<-population
  }
  for(population in populations){
    zPop<-zeros[zeros[[population]]=="X",c("ref","mut","dw","level","isMutT")]
    zPop<-merge(TsTv,zPop, by=c("ref","mut"))
    dw<-data.frame(Positive=nrow(zPop[zPop$dw<0,]),
                   Negative=nrow(zPop[zPop$dw>=0,]),
                   Total = nrow(zPop))
    cat("\n\nZeros from Population", population,"\nDeltas",Dw,":\n")
    print(dw, row.names = F)
    level<-as.data.frame(table(zPop[,c("level")]),
                         stringsAsFactors = F)
    colnames(level)<-c("Levels  ", "Freq")
    print(level, row.names = F)
    typeT<-as.data.frame(table(zPop[,c("type")]),
                      stringsAsFactors = F)
    colnames(typeT)<-c("Mut Type", "Freq")
    print(typeT, row.names = F)
    if(population == "MutT"){
      cat("\tMutT Tv:",nrow(zPop[zPop$isMutT=="*"&
                                   zPop$type=="Tv",]),
          "\tNon MutT Tv:",nrow(zPop[zPop$isMutT!="*"&
                                       zPop$type=="Tv",]))
    }
    
  }
}


countTsTv <- function(save = F,
                      #TsTv="Tv",
                      Dw = "Tai",
                      type = "sum",
                      population = "HighLow",
                      workdir ){
  
  setwd(workdir)
  if(Dw == "Tai"){
    deltaW<-read.csv("./AuxFiles/cDeltaWTAI.csv",
                     header = T,
                     stringsAsFactors = F)
  }else if(Dw == "Cai"){
    deltaW<-read.csv("./AuxFiles/cDeltaW86.csv",
                     header = T,
                     stringsAsFactors = F)
  }else{
    deltaW<-calcDwOld(workdir)
  }
  
  
  if(!population%in%c("High","MutT","HighLow","Low","all")){
    stop('Parameter "population" must be "High","MutT", "Low" or "HighLow"')
  }
  TsTv<-read.csv("AuxFiles/TsTv.csv",
                 header = T,
                 stringsAsFactors = F)

  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)

  # population="High"
  # population="MutT"
  # population="HighLow"
  # population="Low"

  mutAllel<-read.csv(paste0("base/a",population,"MutAllDW.csv"),
                     header = T,
                     stringsAsFactors = F)
  #remove stop codons
  mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"|mutAllel$ref == "TAG" ),]

  #só syn
  mutAllelSyn<-mutAllel[mutAllel$Annotation == "synonymous",]
  #atualize dw values 
  mutAllelSyn<-merge(mutAllelSyn,deltaW,
                     by=c("ref","mut"))
  mutAllelSyn$dw.x<-NULL
  colnames(mutAllelSyn)[135]<-"dw"
  
  #acrescenta fator
  mutAllelSyn<-merge(mutAllelSyn,codonUsage[,c(1,3)],by=c("ref"))
  colnames(mutAllelSyn)[136]<-"factor"
  
  #normaliza valores
  tmp<-mutAllelSyn[,7:128]/mutAllelSyn$factor
  tmp[is.na(tmp)]<- 0
  mutAllelNorm<-cbind(mutAllelSyn[,1:6],tmp,mutAllelSyn[129:136])
  
  rm(tmp)
  
  
  mutAllelNorm$maxFreq<-apply(X = mutAllelNorm[,7:128],
                              MARGIN = 1,
                              max)
  cat("Counts of Mutation per type:\n",
      "\tTransition:",nrow(mutAllelNorm[mutAllelNorm$TsTv == "Ts",]),
      "\n\tTranversion:",nrow(mutAllelNorm[mutAllelNorm$TsTv == "Tv",]),
      "\n\tTotal:",nrow(mutAllelNorm))
  
  #unique(mutAllelNorm$Allele[mutAllelNorm$TsTv == "Tv"])
}

readDeltaW<-function(top=10,
                     workdir ){
  
    setwd(workdir)
    #treat top parameter
    top<-as.numeric(top)
    if(is.na(top) | top < 0){
      stop(paste("Parameter top must be a positive number."))
    }else if(top>0 & top<=1){
      #fraction of number of genes
      top<-top*4324
    }
    top<-floor(top)
    
    dwFile<-paste0("./AuxFiles/cDeltaW",top,".csv")
    if(!file.exists(dwFile)){
      cat("Generating cDeltaW top ",top)
      cmd<-paste0("cd ", workdir,";",
                  "cd ../bin;",
                  "python2.7 1.4CalculateDeltaWCAI.py ", top,";")
      system(cmd)
    }
    deltaW<-read.csv(dwFile,
                     header = T,
                     stringsAsFactors = F)
  return(deltaW)
}

compareGroups<-function(figName,
                             save = F,
                             workdir){
  setwd(workdir)
  population<-"High"
  
  populations<-c("High","MutT","Low")
  for(population in populations){
    
    mutAllel<-read.csv(paste0("base/a",population,"MutAllDW.csv"),
                       header = T,
                       stringsAsFactors = F)
    #remove stop codons
    mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"|mutAllel$ref == "TAG" ),]
    mutAllel<-mutAllel[,c("ref","mut","Position","Annotation","TsTv")]
    #só syn
    mutAllel$Annotation[mutAllel$Annotation == "synonymous"]<-"Syn"
    mutAllel$Annotation[mutAllel$Annotation == "missense" | mutAllel$Annotation == "nonsense" ]<-"Non"
    mutAllel$pop<-population
    if(!exists("result")){
      result<-mutAllel
    }else{
      result<-rbind(result,mutAllel)
    }
  }
  nrow(result)
  library(dplyr)
  library(ggplot2)
  result<-result %>% 
    group_by(Annotation,TsTv,pop) %>% 
    summarise(freq = n())
  result$ann<-paste(result$pop,result$Annotation)
  sum(result$freq)
  g<-ggplot() +
    theme_bw()+
    geom_bar(data = result, 
             aes( x = factor( ann ), 
                  y = freq, 
                  fill = TsTv),
             stat = 'identity' )+
    theme(legend.title = element_blank())+ 
    xlab("Populations Groups & Type of Mutation")+
    ylab("Number of Mutations (SNV)")
  print(g)
  if(save){
    ggsave(filename = "../figures/compPopulations.pdf", 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 0.6,
           width = 16.610, height = 7.440, units = "in",
           dpi = 300)
  }
  return(g)
}

comparePopulations<-function(figName,
                        save = F,
                        workdir){
  setwd(workdir)
  population<-"High"
  
  populations<-c("High","MutT","Low")
  for(population in populations){
    
    mutAllel<-read.csv(paste0("base/a",population,"MutAllDW.csv"),
                       header = T,
                       stringsAsFactors = F)
    #remove stop codons
    mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"|mutAllel$ref == "TAG" ),]
    mutAllel<-mutAllel[,c("ref","mut","Position","Annotation","TsTv","pop")]
    #só syn
    mutAllel$Annotation[mutAllel$Annotation == "synonymous"]<-"Syn"
    mutAllel<-mutAllel[mutAllel$Annotation == "Syn",]
    mutAllel$group<-population
    if(!exists("result")){
      result<-mutAllel
    }else{
      result<-rbind(result,mutAllel)
    }
  }
  nrow(result)
  library(dplyr)
  library(ggplot2)
  result<-result %>% 
    group_by(TsTv,pop,group) %>% 
    summarise(freq = n())
  result$ann<-paste0(substr(result$group,1,1),result$pop)
  result$name<-paste0(ifelse(substr(result$pop,1,1)=="m",
                             "Ara-","Ara+"),
                      substr(result$pop,2,2))
  labels<-unique(result$name[order(result$ann)])
  sum(result$freq)
  g<-ggplot() +
    theme_bw()+
    geom_bar(data = result, 
             aes( x = factor( ann ), 
                  y = freq, 
                  fill = TsTv),
             stat = 'identity' )+
    theme(legend.title = element_blank())+ 
    scale_x_discrete(name = "Populations",
                       label = labels)+
    theme(axis.text.x = element_text(angle = 90))+
    ylab("Number of Mutations (SNV)")
  
  df <- data.frame(
    x = c(2.5, 7.5,11.5),
    y = rep(1625,3),
    text = c("High", "Low", "MutT")
  )
  g<-g+
    geom_text(data = df, 
              aes(x, y, label = text),
              hjust = "center")
  print(g)
  if(save){
    ggsave(filename = paste0("../figures/",
                             figName,
                             ".pdf"), 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 0.6,
           width = 16.610, height = 7.440, units = "in",
           dpi = 300)
  }
  return(g)
}

figFunc01<-function(figName,
                    workdir){
  setwd(workdir)
  g2<-compareGroups(save = F,
                   workdir= workdir)
  ggsave(filename = paste0("../figures/",figName,".pdf"), 
         plot = g2, 
         device = "pdf", 
         path = workdir,
         scale = 0.6,
         width = 16.610, height = 7.440, units = "in",
         dpi = 300)
  
}

figFunc02<-function(figName,
                    normBy = "count",
                    separeTsTv = F,
                    workdir){
  setwd(workdir)
  g<-list()
  index<-1
  if(separeTsTv){
    groups<-c("Ts","Tv")
  }else{
    groups<-c("all")
  }
  
  for(pop in c("High", "MutT")){
    for(t in groups){
      g[[index]]<-correlation(population=pop,
                              save = F,
                              normalize = F,
                              normBy = normBy,
                              workdir = workdir,
                              type = "sum",
                              TsTv=t)
      index<-index+1
    }
    g[[index]]<-correlation(population="Low",
                            save = F,
                            normalize = F,
                            normBy = normBy,
                            workdir = workdir,
                            type = "sum",
                            TsTv="all")
  }
  library(cowplot)
  if(separeTsTv){
    g2<-plot_grid(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],
                  labels=c('A', 'B', 'C', 'D','E'),
                  nrow=3,
                  ncol = 2)
  }else{
    g2<-plot_grid(g[[1]],g[[2]],g[[3]],
                  labels=c('A', 'B', 'C'),
                  nrow=3,
                  ncol = 1)
  }
  
  ggsave(filename = paste0("../figures/",figName,".pdf"), 
         plot = g2, 
         device = "pdf", 
         path = workdir,
         #scale = 2.5,
         width = 16.610, height = 7.440, units = "in",
         dpi = 300)
  
}

figFunc03<-function(figName,
                    normBy = "count",
                    separeTsTv = F,
                    workdir){  
  setwd(workdir)
  g<-list()
  index<-1
  if(separeTsTv){
    groups<-c("Ts","Tv")
  }else{
    groups<-c("all")
  }
  
  for(pop in c("High", "MutT")){
    for(t in groups){
      g[[index]]<-correlation(population=pop,
                              save = F,
                              normalize = T,
                              normBy = normBy,
                              workdir = workdir,
                              type = "sum",
                              TsTv=t)
      index<-index+1
    }
    g[[index]]<-correlation(population="Low",
                            save = F,
                            normalize = T,
                            normBy = normBy,
                            workdir = workdir,
                            type = "sum",
                            TsTv="all")
    
  }
  library(cowplot)
  if(separeTsTv){
    g2<-plot_grid(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],
                  labels=c('A', 'B', 'C', 'D','E'),
                  nrow=3,
                  ncol = 2)
  }else{
    g2<-plot_grid(g[[1]],g[[2]],g[[3]],
                  labels=c('A', 'B', 'C'),
                  nrow=3,
                  ncol = 1)
  }
  
  
  ggsave(filename = paste0("../figures/",figName,".pdf"), 
         plot = g2, 
         device = "pdf", 
         path = workdir,
         #scale = 2.5,
         width = 16.610, height = 7.440, units = "in",
         dpi = 300)
}

figFunc04<-function(figName,
                    type="box",
                    normBy = "count",
                    workdir){
  setwd(workdir)
  g<-list()
  index<-1
  pop="High"
  for(pop in c("High","Low", "MutT")){
    #file<-paste0("AuxFiles/bSum",pop,".csv")  
    g[[index]]<-timeline(population = pop,
                         normalize = F, 
                         type = type,
                         save = F, #save option is not working here
                         normType = "perK",
                         mutT = 0,
                         fix = T,
                         TsTv = "all",
                         operation = "sum",
                         normBy = normBy,
                         workdir = workdir)
    index<-index+1
    g[[index]]<-timeline(population = pop,
                         normalize = T, 
                         type = type,
                         save = F, #save option is not working here
                         normType = "perK",
                         mutT = 0,
                         fix = T,
                         TsTv = "all",
                         operation = "sum",
                         normBy = normBy,
                         workdir = workdir)
    index<-index+1
  }
  library(cowplot)
  g2<-plot_grid(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],g[[6]],
                labels=c('A', 'B', 'C', 'D','E','F'),
                nrow=3,
                ncol = 2)
  # if(type == "box"){
  #   filename=paste0("../figures/FigS05.pdf")
  # }else if(type == "line"){
  #   filename=paste0("../figures/FigS04.pdf")
  # }
    filename=paste0("../figures/",figName,".pdf")
  ggsave(filename = filename, 
         plot = g2, 
         device = "pdf", 
         path = workdir,
         #scale = 2.5,
         width = 16.610, height = 7.440, units = "in",
         dpi = 300)
  g2
}

figFunc05<-function(figName,
                    normBy = "count",
                    type="box",
                    workdir){
  setwd(workdir)
  g<-list()
  index<-1
  pop="High"
  #file<-paste0("AuxFiles/bSum",pop,".csv")  
  
  for(type in c("line","box")){
    for (normalize in c(F,T)) {
      g[[index]]<-timeline(population = pop,
                           normalize = normalize, 
                           normBy = normBy,
                           type = type,
                           save = F, 
                           normType = "perK",
                           mutT = 0,
                           fix = T,
                           TsTv = "all",
                           operation = "sum",
                           workdir = workdir)
      index<-index+1
    }
  }
  #old.par <- par(no.readonly = TRUE)
  library(cowplot)
  g2<-plot_grid(g[[1]],g[[2]],g[[3]],g[[4]],
                labels=c('A', 'B', 'C', 'D'),
                nrow=2,
                ncol = 2)
  filename=paste0("../figures/",figName,".pdf")
  ggsave(filename = filename, 
         plot = g2, 
         device = "pdf", 
         path = workdir,
         #scale = 2.5,
         width = 16.610, height = 7.440, units = "in",
         dpi = 300)
  g2
  dev.off()
}

figFunc06<-function(figName,
                    Dw ,
                    normBy,
                    workdir){
  library(ggplot2)
  library(gridExtra)
  library(gridGraphics)
  setwd(workdir)
  g<-list()
  index<-1

  for(pop in c("HighLow","High","Low", "MutT")){  
    cat("Processing population:",pop,"\n")
  
    g[[index]]<-freqPerDw(population =  pop,
                          save =F, 
                          type = "sum",
                          TsTv = "all",
                          normBy = normBy,
                          Dw = Dw,
                          workdir = workdir) 
      index<-index+1
  }
  lay <- rbind(c(1,2),
               c(3,4))
  
  g2<-grid.arrange(grobs=g,
                  layout_matrix = lay,
                  newpage = F)
  
  filename=paste0("../figures/",figName,".pdf")
  ggsave(filename = filename, 
         plot = g2, 
         device = "pdf", 
         path = workdir,
         #scale = 2.5,
         width = 16.610, height = 7.440, units = "in",
         dpi = 300)
}



listGenes<-function(top=86,
                    workdir){
  setwd(workdir)
  genesAll<-paste0("AuxFiles/listCai",4324,"Genes.csv")
  if(!file.exists(genesAll)){
    cat("Generating cDeltaW for all genes")
    cmd<-paste0("cd ", workdir,";",
                "cd ../bin;",
                "python2.7 1.4CalculateDeltaWCAI.py ", 4324,";")
    system(cmd)
  }
  genesList<-read.csv(genesAll,
                      header = F,
                      stringsAsFactors = F)
    
  riboAll=nrow(genesList[grep("ribosomal",genesList$V2,ignore.case = T),])
  membAll=nrow(genesList[grep("membrane",genesList$V2,ignore.case = T),])
  elongAll=nrow(genesList[grep("elongation factor",genesList$V2,ignore.case = T),])
  
  genes<-paste0("AuxFiles/listCai",top,"Genes.csv")
  if(!file.exists(genes)){
    cat("Generating cDeltaW top ",top)
    cmd<-paste0("cd ", workdir,";",
                "cd ../bin;",
                "python2.7 1.4CalculateDeltaWCAI.py ", top,";")
    system(cmd)
  }
  
  topGenes<-read.csv(genes,
                  header = F,
                  stringsAsFactors = F)
  ribo=nrow(topGenes[grep("ribosomal",topGenes$V2,ignore.case = T),])
  memb=nrow(topGenes[grep("membrane",topGenes$V2,ignore.case = T),])
  elong=nrow(topGenes[grep("elongation factor",topGenes$V2,ignore.case = T),])
  
  cat("Top",top,"genes:\n\tRibossomal Protein -",
      ribo,"\n\t\t",round(ribo/top*100,2),"% of", top,
      "\n\t\t",round(ribo/riboAll*100,2),"% of ribossomal proteins",
      "\n\tMembrane Proteins -",
      memb,"\n\t\t",round(memb/top*100,2),"% of", top,
      "\n\t\t",round(memb/membAll*100,2),"% of membrane proteins",
      "\n\tElongation Factor -",
      elong,"\n\t\t",round(elong/top*100,2),"% of", top,
      "\n\t\t",round(elong/elongAll*100,2),"% of elongation factor")
  
}

unbalancedMutT<-function(top=86,
                         workdir){
  setwd(workdir)
  
  deltaWT<-read.csv("./AuxFiles/cDeltaWTAI.csv",
                    header = T,
                    stringsAsFactors = F)
  deltaWC<-readDeltaW(top,workdir)
  
  mutT<-read.csv(file = "AuxFiles/dMutTmutations.csv",
                 header = T,
                 stringsAsFactors = F)
  deltas<-merge(deltaWC, deltaWT, by=c("ref","mut"))
  colnames(deltas)<-c("ref","mut","dwC","dwT")
  
  deltas$MT<-0
  deltas$MT[paste(deltas$ref,deltas$mut)%in% paste(mutT$ref,mutT$mut) ]<-1
  deltas<-deltas[deltas$MT == 1,]
  cat("MutT\ndwCai\t Positives:",nrow(deltas[deltas$dwC >=0,]),
      "\n\tNegatives",nrow(deltas[deltas$dwC < 0,]),
      "\ndwTai\t Positives:",nrow(deltas[deltas$dwT >=0,]),
      "\n\tNegatives",nrow(deltas[deltas$dwT < 0,]))
  
}
  
countPossibleTsTv<-function(workdir){
  setwd(workdir)
  TsTv<-read.csv("AuxFiles/TsTv.csv",
                 header = T,
                 stringsAsFactors = F)
  cat("Possible Syn SNV Ts and Tv:\n\tTs -",nrow(TsTv[TsTv$type == "Ts",]),
      "\n\tTv -",nrow(TsTv[TsTv$type == "Tv",]))
  
}


plotEnrTopRange<-function(type="High",
                          pval=0.001,
                          quant=5, 
                          normalize = T, 
                          normBy = "count",
                          rank=200,
                          fix = T,
                          TsTv = "Ts",
                          title = T,
                          save=F,
                          Dw = "Cai",
                          topRange=c(50,500,10),
                          refLine=200,
                          smooth = F,
                          figName = 'teste',
                          workdir){
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)  
  library(doParallel)
  
  setwd(workdir)
  if(!file.exists(paste0("./base/a",type,"MutAllDW.csv"))){
    if(type == "HighLow"){
      joinHighLow(workdir)
    }else{
      stop('Parameter "type" must be 
          \t"High" for High mutations rate populations;
          \t"Low" for Low mutations rate populations;
          \t"MutT" for MutT mutations populations; or
          \t"HighLow" for High and Low mutations together.' )
      #return(0)
    }
  }
  
  rankStd<-rank
  if(Dw == "Tai"){
    stop('Sorry! Has no sense run this for TAI')
  }else{
    if(!is.numeric(topRange) | length(topRange) != 3){
      stop('Parameter "topRenge"must be a numeric vector with length 3 (start, end, and step)')
    }
  }
  
  if(title){
    titulo<-paste("Top",rankStd,"Frequencies, ",
                  ifelse(normalize,"normalized, ","without normalization, "),
                  type, ", ",toupper(Dw),", ")
  }else{
    titulo=" "
  }
  if(!fix%in%c("all",T,F) ){
    cat('Argument "fix" must be TRUE, FALSE or "all".')
    return()
  }
  if(fix!="all"){
    if(fix){
      titulo<-paste(titulo,"Fixed Mutations")
      name<-"Fix"
    }else{
      titulo<-paste(titulo,"Non-fixed Mutations")
      name<-"Nonfix"
    }
  }else{
    titulo<-paste(titulo,"All Mutations")
    name<-"All"
  }
  
  if(!TsTv%in%c("all","Ts","Tv")){
    stop('Use values "all","Ts" (transitions), or "Tv" (transversions) for TsTv')
    #return(0)
  }
  resPlot<-data.frame(top=numeric(),
                      faixa=numeric(),
                      value=numeric(),
                      color=factor())
  ini<-topRange[1]
  fim<-topRange[2]
  step<-topRange[3]
  cat("Please wait. This can take a while...\n")
  pb <- txtProgressBar(min = ini, max = fim, style = 3)
  cl <- makeCluster(10)
  registerDoParallel(cl)
  # foreach ----
  resPlot<-foreach(top = seq(ini,fim,step),
                   .combine=rbind,
                   .export = c("enrich",
                               "readDeltaW",
                               "pb",
                               "quantilizeDw")) %dopar% {
                     
  # top=86
  # for(top in seq(ini,fim,step)){
    #depleted
    resCorrigido<-enrich(normalize = normalize, 
                         normBy = normBy,
                         pval = pval,
                         rank = rankStd,
                         type = type, 
                         tail = "L",
                         quant = quant,
                         fix = fix,
                         TsTv = TsTv,
                         Dw = Dw,
                         top = top,
                         workdir = workdir)
    
    vdw<-resCorrigido[[2]]
    resCorrigido<-resCorrigido[[1]]
    resCorrigido[is.na(resCorrigido)]<-0
    resCorrigido[resCorrigido!=0]<-1
    values<-apply(X=resCorrigido,
                  MARGIN = 1,
                  FUN = sum,na.rm = T)
    faixas<-1:quant
    tops<-rep(top,quant)
    resPlotTmp<-data.frame(top=tops,
                           faixa=-faixas,
                           value=values,
                           color=as.factor(faixas))
    # 
    # resPlot<-rbind(resPlot,
    #                data.frame(top=tops,
    #                           faixa=-faixas,
    #                           value=values,
    #                           color=as.factor(faixas)))
    #enriched
    resCorrigido<-enrich(normalize = normalize, 
                         normBy = normBy,
                         pval = pval,
                         rank = rankStd,
                         type = type, 
                         tail = "H",
                         quant = quant,
                         fix = fix,
                         TsTv = TsTv,
                         Dw = Dw,
                         top = top,
                         workdir = workdir)
    
    vdw<-resCorrigido[[2]]
    resCorrigido<-resCorrigido[[1]]
    resCorrigido[is.na(resCorrigido)]<-0
    resCorrigido[resCorrigido!=0]<-1
    values<-apply(X=resCorrigido,
                  MARGIN = 1,
                  FUN = sum,na.rm = T)
    faixas<-1:quant
    tops<-rep(top,quant)
    resPlotTmp<-rbind(resPlotTmp,
                   data.frame(top=tops,
                              faixa=faixas,
                              value=values,
                              color=as.factor(faixas)))
    setTxtProgressBar(pb, top)
    print(resPlotTmp)
  }    
  #maxVal<-max(max(resPlot$value),max(-resPlot$value))
  #maxY<-ceiling(maxVal/10)*10
  #zeroNeg<- -maxY+2
  
  resPlot$type<-"Enriched"
  resPlot$type[resPlot$faixa<0]<- "Depleted"
  resPlot$type<-factor(resPlot$type,
                       levels = c("Enriched", "Depleted"))
  
  
#  resPlot$value[resPlot$faixa<0]<- resPlot$value[resPlot$faixa<0]+zeroNeg
  library(ggplot2)
  g<-ggplot()+theme_bw()+
    theme(legend.title = element_text(size = 9))+
    facet_grid(type~.)+
    scale_x_continuous(name = "Ranked Genes")+
    scale_y_continuous(name = "Significant generations",
                       minor_breaks = NULL)+
    scale_color_discrete(name = "Levels")
  
  if(smooth){
    g<-g+geom_smooth(data = resPlot,
                     aes(x=top,y=value,
                         group=faixa, 
                         color=as.factor(color)),
                     size =0.5,
                     span = 0.2,se = F)
  }else{
    g<-g+geom_line(data = resPlot,
                   aes(x=top,y=value,
                       group=faixa, 
                       color=as.factor(color)))
  }
  if(refLine>0){
    g<-g+geom_vline( xintercept = refLine, lty=2, col="red")
  }
  
  print(g)
  if(save){
    if(figName == ''){
      filename <- paste0("../figures/enrich",type,rank,
                         Dw,".pdf")
    }else{
      filename <-paste0("../figures/",figName,".pdf")
    }
    ggsave(filename = filename, 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 0.6, 
           width = 16.610, height = 7.440, units = "in",
           dpi = 300)
  }
}

plotEnrRankRange<-function(type="HighLow",
                          pval=0.001,
                          quant=5, 
                          normalize = T, 
                          normBy = "count",
                          rankRange=c(50,100,5),
                          fix = "all",
                          TsTv = "all",
                          title = F,
                          save=F,
                          Dw = "Cai",
                          top= 86,
                          refLine=200,
                          smooth = F,
                          figName = '',
                          workdir){
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)
  library(doParallel)
  setwd(workdir)
  if(!file.exists(paste0("./base/a",type,"MutAllDW.csv"))){
    if(type == "HighLow"){
      joinHighLow(workdir)
    }else{
      stop('Parameter "type" must be 
          \t"High" for High mutations rate populations;
          \t"Low" for Low mutations rate populations;
          \t"MutT" for MutT mutations populations; or
          \t"HighLow" for High and Low mutations together.' )
      #return(0)
    }
  }
  
  #rankStd<-rank
  if(Dw == "Tai"){
    top<-""
  }else{
    if(!is.numeric(top)){
      stop('Parameter "top" must be  numeric.')
    }
  }
  
  if(title){
    titulo<-paste("Top","Frequencies, ",
                  ifelse(normalize,"normalized, ","without normalization, "),
                  type, ", ",toupper(Dw),", ")
  }else{
    titulo=" "
  }
  if(!fix%in%c("all",T,F) ){
    cat('Argument "fix" must be TRUE, FALSE or "all".')
    return()
  }
  if(fix!="all"){
    if(fix){
      titulo<-paste(titulo,"Fixed Mutations")
      name<-"Fix"
    }else{
      titulo<-paste(titulo,"Non-fixed Mutations")
      name<-"Nonfix"
    }
  }else{
    titulo<-paste(titulo,"All Mutations")
    name<-"All"
  }
  
  if(!TsTv%in%c("all","Ts","Tv")){
    stop('Use values "all","Ts" (transitions), or "Tv" (transversions) for TsTv')
    #return(0)
  }
  resPlot<-data.frame(top=numeric(),
                      faixa=numeric(),
                      value=numeric(),
                      color=factor())
  ini<-rankRange[1]
  fim<-rankRange[2]
  step<-rankRange[3]
  cat("Please wait. This can take a while...\n")
  pb <- txtProgressBar(min = ini, max = fim, style = 3)
  
  cl <- makeCluster(10)
  registerDoParallel(cl)
  # foreach ----
  resPlot<-foreach(rankStd = seq(ini,fim,step),
          .combine=rbind,
          .export = c("enrich",
                      "readDeltaW",
                      "pb",
                      "quantilizeDw")) %dopar% {
    #depleted
    resCorrigido<-enrich(normalize = normalize, 
                         normBy =normBy,
                         pval = pval,
                         rank = rankStd,
                         type = type, 
                         tail = "L",
                         quant = quant,
                         fix = fix,
                         TsTv = TsTv,
                         Dw = Dw,
                         top = top,
                         workdir = workdir)
    
    vdw<-resCorrigido[[2]]
    resCorrigido<-resCorrigido[[1]]
    resCorrigido[is.na(resCorrigido)]<-0
    resCorrigido[resCorrigido!=0]<-1
    values<-apply(X=resCorrigido,
                  MARGIN = 1,
                  FUN = sum,na.rm = T)
    faixas<-1:quant
    tops<-rep(rankStd,quant)
    resPlotTmp<-data.frame(top=tops,
                              faixa=-faixas,
                              value=values,
                              color=as.factor(faixas))
    #enriched
    resCorrigido<-enrich(normalize = normalize, 
                         normBy=normBy,
                         pval = pval,
                         rank = rankStd,
                         type = type, 
                         tail = "H",
                         quant = quant,
                         fix = fix,
                         TsTv = TsTv,
                         Dw = Dw,
                         top = top,
                         workdir = workdir)
    
    vdw<-resCorrigido[[2]]
    resCorrigido<-resCorrigido[[1]]
    resCorrigido[is.na(resCorrigido)]<-0
    resCorrigido[resCorrigido!=0]<-1
    values<-apply(X=resCorrigido,
                  MARGIN = 1,
                  FUN = sum,na.rm = T)
    faixas<-1:quant
    tops<-rep(rankStd,quant)
    resPlotTmp<-rbind(resPlotTmp,
                   data.frame(top=tops,
                              faixa=faixas,
                              value=values,
                              color=as.factor(faixas)))
    setTxtProgressBar(pb, rankStd)
    print(resPlotTmp)
  }    
  stopCluster(cl)
  resPlot$type<-"Enriched"
  resPlot$type[resPlot$faixa<0]<- "Depleted"
  resPlot$type<-factor(resPlot$type,
                       levels = c("Enriched", "Depleted"))
  
  #maxVal<-max(max(resPlot$value),max(-resPlot$value))
  #maxY<-ceiling(maxVal/10)*10
  #zeroNeg<- -maxY+2
  #resPlot$value[resPlot$faixa<0]<- resPlot$value[resPlot$faixa<0]+zeroNeg
  library(ggplot2)
  g<-ggplot()+theme_bw()+
    theme(legend.title = element_text(size = 9))+
    facet_grid(type~.)+
    scale_x_continuous(name = "Ranked Mutations")+
    scale_y_continuous(name = "Significant generations",
                       minor_breaks = NULL)+
    scale_color_discrete(name = "Levels")
  
  if(smooth){
    g<-g+geom_smooth(data = resPlot,
                     aes(x=top,y=value,
                         group=faixa, 
                         color=as.factor(color)),
                     size =0.5,
                     span = 0.2,se = F)
  }else{
    g<-g+geom_line(data = resPlot,
                     aes(x=top,y=value,
                         group=faixa, 
                         color=as.factor(color)))
  }
  if(refLine>0){
    g<-g+geom_vline( xintercept = refLine, lty=2, col="red")
  }
    
  print(g)
  if(save){
    if(figName == ''){
      filename <- paste0("../figures/enrich",type,rank,
                         Dw,".pdf")
    }else{
      filename <-paste0("../figures/",figName,".pdf")
    }
    ggsave(filename = filename, 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 0.6, 
           width = 16.610, height = 7.440, units = "in",
           dpi = 300)
  }
}

plotExpression<-function(workdir = workdir){
  setwd(workdir)
  expression<-genesExpression(workdir)
  expression$x<-seq(1,nrow(expression),1)
  library(ggplot2)
  ggplot()+theme_bw()+
    geom_line(data = expression,
              aes(x=x, y=median),
              col="blue")
  
  
}

quantilizeDw<-function(deltaW,
                          quant,
                          Dw,
                          quiet = T){
  #deltas = dataframe with ref, mut and dw values
  #quant =  number of quantiles
  #Dw is Tai or Cai
  #Tai distribution needs a special process to quantilize
  # because it have a lot of repeated values
  deltas<-deltaW$dw
  if(Dw == "Cai"){
    quantil<-quantile(deltas,c(seq(from = 0,to = 1,by = 1/quant)),type = 7)
  }else{
    deltas<-deltas[order(deltas)]
    fIdeal<-length(deltas)/quant
    quantil<-quantile(deltas,c(seq(from = 0,to = 1,by = 1/quant)),type = 1)
    i=5
    for (i in 2:quant) {
      vtop<-max(which(deltas==quantil[i]))
      vbot<-min(which(deltas==quantil[i]))-1
      if(abs(vtop-fIdeal*(i-1)) <= abs(fIdeal*(i-1) - vbot)){
        quantil[i]<-deltas[vtop]
      }else{
        quantil[i]<-deltas[vbot]
      }
      if(round(quantil[i],10) == round(quantil[i-1],10)){
        quantil[i]<-deltas[vtop]
      }
    }
    #prepare a named vector with quantiles
    names<-c("0%")
    for (i in 2:quant) {
      names<-c(names,
               paste0(round(length(deltas[deltas<=quantil[i]])/length(deltas)*100,2),"%"))
    }
    names(quantil)<-c(names,"100%")
  }
  #create Dw table with levels
  #quantDw<-as.data.frame(deltas)
  quantDw<-deltaW
  qn = 1
  for (qn in 1:(length(quantil)-1)) {
    if(qn==1){
      limDown<-quantil[qn]-0.1
    }else{
      limDown<-quantil[qn]
    }
    #count just not 0
    quantDw$faixa[quantDw$dw > limDown &
                     quantDw$dw <= quantil[qn+1]]<- qn
  }
  #check if is a balanced distribution of Dw
  level = 1
  for(level in 1:round(quant/2,0)){
    level2<-quant-level+1#opposite side of Dw
    negDw<-quantDw$dw[quantDw$faixa == level]
    quantDw$faixa[quantDw$faixa %in% -negDw]<- level2
  }
    
  #Show counts of levels
  if(!quiet){
    cat("Quantiles: \n")
    print(quantil)
    cat("Levels: \n")
    qn=2
    for (qn in 1:(length(quantil)-1)) {
      qtd<-nrow(quantDw[quantDw == qn,])
      msg<-paste0("\t",
                  qn,"-\t",
                  round(quantil[qn],2)," to ",
                  round(quantil[qn+1],2),"\t- ",qtd," values\n")
      cat(msg)
    }
    
  }
  
  #colnames(quantDw)<-c( "dw","faixa")
  #quant<-length(quantil)-1
  return(list(quantDw,quantil))
}

normalization<-function(data,
                        normBy,
                        melted = T,
                        workdir = workdir){
  if(!normBy%in%c("CUB","count","mean","none")){
    stop('Normalizations can be done only by Codon Usage Bias ("CUB"), mutation count ("count"), "mean" for mean of mutations in a generation, or "none" for no normalization.')
  }
  library(reshape2)
  setwd(workdir)
  data[is.na(data)]<-0
  qt<-data
  qt$sum<-apply(qt[,-1],
                MARGIN = 1,
                FUN = sum)
  #just exclude lines with no frequencies
  qt<-qt[qt$sum !=0,c("ref","sum")]
  
  qt<-qt%>%
    group_by(ref)%>%
    summarize(count=n())
  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  codonUsage<-codonUsage[order(codonUsage$factorPerK),]
  
  codonOrder<-cbind(codonUsage[c("ref","factorPerK")], 
                    seq(1,nrow(codonUsage),1))
  colnames(codonOrder)<-c("ref","factor","ordem")
  
  qt<-merge(qt,codonOrder,by="ref")
  qt<-na.exclude(qt)
  qtFactors<-qt
  colnames(qt)<-c("ref","y","factor","ordem")
  colnames(qtFactors)<-c("ref","count","CUB","ordem")
  for( col in colnames(data[-1])){
    tmpTime<-data[data[,col]>0,c("ref",col)]
    tmpTime<-tmpTime%>%
      group_by(ref)%>%
      summarize(count=n())
    colnames(tmpTime)<-c("ref",col)
    qtFactors<-merge(qtFactors,tmpTime, 
                     by= "ref",
                     all.x = T)
  }
  
  qt<-data
  qt$order<-c(1:nrow(qt))
  qt<-melt(qt,id.vars = c("ref","order"))
  if(normBy =="CUB"){
    qt<-merge(qt,qtFactors[,c(1:4)],by=c("ref"))
    qt$value<-qt$value/qt$CUB
  }else if(normBy =="count"){
    qt<-merge(qt,qtFactors[,c(1:4)],by=c("ref"))
    qt$value<-qt$value/qt$count
  }else if(normBy =="mean"){
    teste<-melt(qtFactors[,c(1,5:126)],id.vars = c("ref"))
    t2<-merge(qt,teste,by=c("ref","variable"))
    t2$value<-t2$value.x/t2$value.y
    qt<-t2
  }
  if(!melted){
    qt<-dcast(qt,order+ref~variable,fill=0)
    qt<-qt[order(qt$order),]
    qt<-qt[,-1]
  }
  return(qt)
}

corrCUBtRNA<-function(workdir,
                      save=F,
                      figName="corrCUBtRNA"){
  setwd(workdir)
  # deltaW<-readDeltaW(top,workdir)
  # deltaW$dw<-NA
  deltaW<-read.csv("./AuxFiles/cDeltaWTAI.csv",
                   header = T,
                   stringsAsFactors = F)
  
  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  aminos<-read.csv("AuxFiles/dAminoCodon.csv",
                   header = T,
                   stringsAsFactors = F)
  codonUsage<-merge(codonUsage[,c("ref","count")],
                    aminos, 
                    by="ref")
  tRNAe<-read.csv("TAIFiles/REL606tRNA.csv",
                  sep = "\t",
                       header = T,
                       stringsAsFactors = F)
  tRNAe$ref<-apply(as.data.frame(tRNAe$AntiCodonsList),
                     MARGIN = 1,
                     convertToComplement)
  codonUsage<-merge(codonUsage,tRNAe,by="ref",all.y = T)
  codonUsage$ref[is.na(codonUsage$amino)]
  codonUsage<-na.exclude(codonUsage)
  tmp<-as.data.frame(unique(codonUsage$amino))
  colnames(tmp)<-"amino"
  tmp$max<-apply(tmp,MARGIN = 1,
                 FUN = function(x){
                   max(codonUsage$count[codonUsage$amino == x[1]])
                 })
  tmp$min<-apply(as.data.frame(tmp[,1]),MARGIN = 1,
                 FUN = function(x){
                   min(codonUsage$count[codonUsage$amino == x[1]])
                 })
  cor<-c("red"="red",
         "blue"="blue",
         "orange"="orange",
         "darkgreen"="darkgreen")
  linha<-1
  for(i in 1:5){
    for(j in 1:4){
      tmp$cor[linha]<-cor[j]
      tmp$shape[linha]<-i
      linha<-linha+1
    }
  }
  codonUsage<-merge(codonUsage,tmp,by="amino")
  
  codonUsage$shape<-(factor(codonUsage$shape))
  codonUsage$cores<-as.factor(codonUsage$amino)
  
  lm<-lm(codonUsage$tGCN~codonUsage$count)
  corr<-cor.test(codonUsage$tGCN,codonUsage$count)
  maxX<-max(codonUsage$count)
  maxY<-max(codonUsage$tGCN)
  library(ggplot2)
  g<-ggplot()+theme_bw()+
    ylab("tRNA Count")+
    xlab("Codon Usage")+
    geom_text(data = codonUsage,
               aes(y=tGCN,x=count,label=amino,
                   colour=cores),cex=5)+
                   # shape=shape))+ 
    scale_color_manual(guide=F, values = rainbow(20))+
               # col="red", 
               # pch=1)+
    geom_abline(slope = lm$coefficients[2], 
                intercept = lm$coefficients[1],
                lty=2,
                col="blue")+
    geom_text(aes(x=maxX*0.1, y=maxY*0.75,
                  label=paste("Correlation:",
                              round(corr$estimate,3),
                              "\npvalue:",
                              format(corr$p.value,digits = 3, scientific=T))))
  
  g
  
  tmp1<-codonUsage[,1:3]

  amino="A"
  count=0
  for(amino in unique(tmp1$amino)){
    tmp2<-tmp1[tmp1$amino==amino,]
    tmp2<-tmp2[order(tmp2$count),]
    tmp2$use<-paste0(rep("M",nrow(tmp2)),seq(1,nrow(tmp2),1))
    if(count == 0){
      tmp3<-tmp2
      count =1
    }else{
      tmp3<-rbind(tmp3,tmp2)
    }
  }
  codonUsage<-merge(codonUsage,tmp3[,c("ref","use")],by="ref")
  codonUsage$use[codonUsage$count==codonUsage$max]<-"Preferred Codon"
  codonUsage$use[codonUsage$count==codonUsage$min]<-"Less Preferred"
  codonUsage$use<-as.factor(codonUsage$use)

  pos <- position_jitter(width = 0.2, seed = 1)
  g<-ggplot(data = codonUsage,
            aes(x=amino,
                y=tGCN,
                
                color=use, 
                label=ref,
                hjust="inward"))+theme_classic()+
    ylab("tRNA Gene Count")+
    xlab("Amino Acid")+
    geom_tile(col="black",fill=alpha(colour = "white",alpha = 0))+
    geom_point(position = pos)+
    geom_text(position = pos, cex=2)+
    # scale_fill_gradient(guide=F,
    #                     low = alpha(colour = "white",alpha = 0),
    #                     high = alpha(colour = "white",alpha = 0))+
    # scale_color_gradient2(name="Usage",
    #                       low = ("red"), 
    #                       mid = "grey",
    #                       high = ("blue"), 
    #                       midpoint = max(codonUsage$count)/2)+
    scale_y_continuous(breaks=seq(0, 
                                  max(codonUsage$tGCN), 
                                  by = 1))+
    scale_color_brewer(name="Usage",guide="legend",palette="RdBu")
  # shape=shape))+ 
    # col="red", 
    # pch=1)+
  if(save){
    filename=paste0("../figures/",figName,".pdf")
    ggsave(filename = filename, 
         plot = g, 
         device = "pdf", 
         path = workdir,
         #scale = 2.5,
         width = 8, height = 4, units = "in",
         dpi = 300)
  }
  
  g
  
  cat("Zeros:\n")
  tmp2<-t(apply(codonUsage[codonUsage$tGCN ==0,c(2,7,3,6)],
        MARGIN = 1,
        FUN = function(x){
          if(x[2]==x[3]){
            return(c(x[1],"Min"))
          }else if(x[4]==x[3]){
            return(c(x[1],"Max"))
          }else{
            return(c(x[1],"Oth"))
          }
          
        }))
  print(table(tmp2[,2]))
  print(table(tmp2[,1]))
  
  cat("\nZeros in min:",tmp2[tmp2[,2]=="Min",1])
  cat("\nZeros in max:",tmp2[tmp2[,2]=="Max",1])
  cat("\nTotal:",nrow(codonUsage[codonUsage$tGCN ==0,]))
  cat("\nNon Zeros:",unique(codonUsage$amino[codonUsage$amino%in%tmp$amino[!tmp$amino%in%unique(tmp2[,1])]]))
  # 
  # teste<-deltaW
  # teste<-merge(deltaW,codonUsage[,c(1,2,5)], by="ref")
  # colnames(teste)<-c("ref","mut","dw","amino","tRef" )
  # teste<-merge(teste,codonUsage[,c(1,2,5)], by.x = "mut",by.y="ref")
  # colnames(teste)<-c("mut","ref","dw","amino","tRef","tMut" )
  # teste<-teste[,c("ref","mut","dw","amino","tRef","tMut" )]  
  # table(teste[,c(5,6)])
  # 
  # library(scales)
  # g<-ggplot()+theme_bw()+
  #   ylab("tRNA count - Mutation")+
  #   xlab("tRNA count - Reference")+
  #   geom_tile(data = teste,
  #             aes(y=tMut,x=tRef,fill=dw,color=dw))+
  #   geom_jitter(data = teste,
  #             aes(y=tMut,x=tRef,color=dw),
  #             width = 0.4)+
  #   scale_fill_gradient2(guide=F,low = alpha(colour = "white",alpha = 0),
  #                        mid = alpha(colour = "white",alpha = 0),
  #                       high = alpha(colour = "white",alpha = 0), midpoint = 0)+
  #   scale_color_gradient2(low = ("red"), mid = "lightgrey",
  #                      high = ("blue"), midpoint = 0)
  # # shape=shape))+ 
  # # col="red", 
  # # pch=1)+
  # 
  # g
  
}

