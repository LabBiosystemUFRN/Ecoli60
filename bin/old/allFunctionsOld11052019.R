# !diagnostics off
calcDwOld<-function(top=86,
                    workdir="/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
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
                    workdir="/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
  cat("Using",workdir)
  #setwd(workdir)
  deltaW<-readDeltaW(top,workdir)
  deltaWOld<-calcDwOld(top,workdir)
  colnames(deltaWOld)<-c("ref", "mut", "dwOld" )
  # codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
  #                      header = T,
  #                      stringsAsFactors = F)

  tmp<-merge(deltaW,deltaWOld,by=c("ref","mut"))

  # tmp<-merge(codonUsage[,c("ref","count")],deltaW[,c("ref","mut","dw")],by="ref")
  # tmp<-merge(tmp,codonUsage[,c("ref","count")],by.x="mut",by.y="ref")
  # colnames(tmp)<-c("mut","ref","countRef","dw","countMut")
  # tmp<-tmp[,c("ref","mut","countRef","countMut","dw")]
  # tmp$codonRatio<-log10(tmp$countMut/tmp$countRef)
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
           scale = 1.5,
           #width = 3*maxX, 
           #height = 3*maxY, 
           #units = "in",
           dpi = 300)
  }
  
  #return(p)
}

mutPersist<- function(top=86,
                      type="High",
                  pval=0.01,
                  quant=7, 
                  normalize = F, 
                  rank=100, 
                  tail="L",
                  workdir="/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
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


joinHighLow<- function(workdir="/home/clovis/Dropbox/Ecoli60/data_files/"){
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
            workdir="/home/clovis/Doutorado/Projetos/Ecoli60/data_files/")
}


totalCodons<-function(workdir="/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
  setwd(workdir)
  deltaW<-read.csv("AuxFiles/cDeltaW.csv",
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
    
    # somat<-mutAllel[0,c(127,128,5:126)]
    # probab<-mutAllel[0,c(127,128,5:126)]
    somat<-mutAllel[0,c(1:2,7:128)]
    probab<-mutAllel[0,c(1:2,7:128)]
    #i=2
    for(i in 1:nrow(codons)){
      tmp<-subset(x= mutAllel,subset = (ref == codons$ref[i] & 
                                          mut == codons$mut[i]),
                  select = X0:X60500)
      #remove NA
      tmp[is.na(tmp)]<-0
      # tmp<-subset(x= mutAllel,
      #             select = X0:X60500)
      tmp2<-apply(X = tmp,MARGIN = 2, FUN = function(x){
        #zero<-1e-15
        #one<- 1-zero
        #x<-c(.9,.9,-1,0.4)
        x<-x[x!=-1]
        x<-1-x
        #x[x==1]<-one
        #x[x==0]<-zero
        #print(x)
        return(1-prod(x))
      })
      
      #probab<-cbind(mutAllel[,c(1,2,3,127,128)],tmp2)
      
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

mutationsMutT<-function(){
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

joinHighLow<- function(workdir="/home/clovis/Dropbox/Ecoli60/data_files/"){
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
}

enrich<- function(type="HighLow",
                  pval=0.01,
                  quant=7, 
                  normalize = F, 
                  rank=0.2, 
                  tail="L",
                  fix="all",
                  TsTv ,
                  Dw = "Tai",
                  top = '',
                  workdir = "/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
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
  
  #atualize dw values 
  mutAllel<-merge(mutAllel,deltaW,
                  by=c("ref","mut"))
  mutAllel$dw.x<-NULL
  colnames(mutAllel)[135]<-"dw"
  #filtra Ts ou Tv
  if(TsTv != "all"){
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
  cat("Using the top",rankStd,"frequencies \n")
  
  
  top100Base<-mutAllel[,c(1:5,135)]#,2,127,128,3)]
  #colnames(top100Base)<-c("Position","Gene","ref","mut","Allele")
  #top100Base<-merge(top100Base,deltaW[,c(1,2,4)],by=c("ref","mut"))
  top100Base<-merge(top100Base,codonUsage[,c(1,3)],by=c("ref"))
  colnames(top100Base)[colnames(top100Base) == 'factorPerK'] <- 'factor'
  #colnames(top100Base)<-c("ref","mut","Position","Gene","Allele","dw", "factor")
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
      top100Base$faixa[top100Base$dw == quantil[1]]<-1
    }
    top100Base$faixa[top100Base$dw> quantil[i] & top100Base$dw<=quantil[i+1]]<-i
  }
  #garante a mesma ordem me mutAllel e top100Base
  top100Base<-top100Base[order(top100Base$Position),c("Position","Gene","Allele","ref","mut","dw","factor","faixa")]
  # countsTot<-as.data.frame(table(top100Base[,c("faixa")]))
  # countsTot<-countsTot[countsTot$Freq!=0,]
  # colnames(countsTot)<- c("faixa","white")
  # countsTot$black<-nrow(top100Base)-countsTot$white
  
  if(exists("result")){rm(result)}
  #loop ----
  col="X3000"
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
      #normalize observations
      top100$cNorm<-apply(X = top100[,c("factor","count")],
                          MARGIN = 1,
                          FUN = function(x){
                            return(x[2]/x[1])
                          })
    }else{
      top100$cNorm<-top100$count
    }
    #    top100<-top100[top100$cNorm>0,]
    # countsTot<-as.data.frame(table(top100[top100$cNorm>0,
    #                                       c("faixa")]),
    #                          stringsAsFactors = F)
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
  cat("Minimal frequence used:", minFreq,"\n")
  #resCorrigido<-result[0,2:ncol(result)]
  resCorrigido<-result[,2:ncol(result)]
  # linha=1
  # for( linha in faixas) {
  #   resCorrigido[nrow(resCorrigido)+1,]<-c(p.adjust(result[linha,2:ncol(result)]))
  # }
  resCorrigido[resCorrigido>pval]<-NA
  return(list(resCorrigido,vdw,rankStd))
}


plotEnrDeplPVal<-function(type="High",
                          pval=0.01,
                          quant=7, 
                          normalize = F, 
                          rank=100,
                          fix = "all",
                          TsTv = "all",
                          title = T,
                          save=F,
                          Dw = "Tai",
                          top=100,
                          workdir = "/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)  
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
  
  resPlot<-resCorrigido[,1:ncol(resCorrigido)]
  #adiciona indice de faixas
  #resPlot$faixa<- c(-1:-quant)
  resPlot$faixa<- -(quant+1)+c(1:quant)
  #elimina linhas sem resultado
  resPlot<-resPlot[rowSums(!is.na(resPlot[,1:ncol(resCorrigido)])) > 0,]
  resCorrigido<-enrich(normalize = normalize, 
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
  tmp[tmp==0]<-NA
  #adiciona indice de faixas
  tmp$faixa<-c(1:quant)
  #elimina linhas sem resultado
  tmp<-tmp[rowSums(!is.na(tmp[,1:ncol(resCorrigido)])) > 0,]
  
  #une resultados
  resPlot<-rbind(resPlot, tmp)
  
  #elimina colunas sem resultado
  resPlot<-resPlot[colSums(!is.na(resPlot)) > 0]
  
  if(sum(dim(resPlot))== 0){
    warning("No enrichments or depletions detected!")
    return(invisible(0))
  }
  
  
  min<-min(resPlot[!is.na(resPlot) & resPlot >0])
  val<-seq(round(log10(min),0),round(log10(pval),0),1)
  
  mydf <- data.frame(id = rep(1, length(val)), 
                     pval = val,
                     cor=rainbow(length(val),end=0.7 ),#topo.colors
                     #cor=colorRampPalette(c("blue","green", "red"))( length(val) ),
                     stringsAsFactors = F)
  
  #barplot(mydf$pval,col=mydf$cor)
  legenda<-
    ggplot(mydf) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          #axis.text.y=element_blank(),
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
  
  
  linha =  seq(-0.4,(-quant*0.05)-0.4,-0.05)#c(-0.5,-0.6)
  contLin = 0
  x<-strtoi(sapply(X = colnames(resPlot[,1:(ncol(resPlot)-1)]),
                   substr,2,30))/1000
  
  p<-ggplot()+
    theme_bw()+
    xlim(0,max(x))+
    # ylim(-7,7)+
    xlab("Generations (x1000)")+
    ggtitle(titulo) +
    scale_y_continuous(name =bquote(Delta~"w levels"), 
                       label=c(c(1:quant)," ",c(1:quant)),
                       #label=c(c(quant:1)," ",c(1:quant)),
                       #"7","6","5","4","3","2","1"," ",
                       #        "1","2","3","4","5","6","7"),
                       breaks= c(-quant:quant),
                       minor_breaks = NULL,
                       limits = c(-quant,quant))+
    scale_color_identity()+
    geom_hline(yintercept = 0,col="white")+
    geom_hline(yintercept = 0.1,col="black")+
    geom_hline(yintercept = -0.1,col="black")+
    theme(#panel.grid.major.y = element_blank(), 
      panel.grid.minor.y = element_blank(),
      axis.ticks.y = element_blank())
  
  p
  i=1
  
  for(i in 1:nrow(resPlot)){
    qt1<-data.frame(x=x,#seq(0,60,0.5),
                    y=rep(resPlot$faixa[i],length(x)),
                    pval=t(round(log10(resPlot[i,1:(ncol(resPlot)-1)]),0)),
                    stringsAsFactors = F)
    colnames(qt1)<-c("x","y","pval")
    qt1<-merge(qt1,mydf[,2:3], by = "pval")
    p<-p+geom_point(data = qt1,
                    aes(x, y, color = (cor)),
                    pch=15)
  }
  df <- data.frame(
    x = c(5, 5),
    y = c(-quant/2, quant/2),
    text = c("Depleted", "Enriched")
  )
  p<-p+
    geom_text(data = df, aes(x, y, label = text))
  
  txt<-textGrob("p-value\n (log10)",
                gp=gpar(fontsize=10))
  
  lay <- rbind(c(1,1,1,1,1,1,1,3),
               c(1,1,1,1,1,1,1,2),
               c(1,1,1,1,1,1,1,2),
               c(1,1,1,1,1,1,1,NA))
  g<-grid.arrange(p,legenda,txt, layout_matrix = lay)
  
  if(save){
    ggsave(filename = paste0("../figures/enrich",type,rankStd,
                             name,
                             Dw,".pdf"), 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 1.5, 
           #width = 6.72, height = 2.98, units = "in",
           dpi = 300)
  }
  
  
  #return(g)
  }

timelineOld<-function(file="bSumHigh.csv",
                   quant = 7, 
                   normalize = T,
                   type = "box",
                   save = F,
                   normType = "perK",
                   mutT = 0,
                   workdir = "/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
  dirFig<-"/home/clovis/Doutorado/Projetos/Ecoli60/figures"
  setwd(workdir)
  somat<-read.csv(file = file,
                  stringsAsFactors = F)
  somat<-na.exclude(somat)
  
  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  #usado no cálculo da frqeuencia
  # allFreq<-sum(codonUsage$count)
  # codonUsage$freqPerK<- codonUsage$count*1000/allFreq
  # colnames(codonUsage)<-c("ref","count","factorAbs","factor")
  codonUsage<-codonUsage[order(codonUsage$factorPerK),]
  # write.csv(x=codonUsage,file="dCodonUsage.csv",
  #           row.names = F)
  
  
  deltaW<-read.csv("AuxFiles/cDeltaW.csv",
                   header = T,
                   stringsAsFactors = F)
  MutT<-read.csv(file = "AuxFiles/dMutTmutations.csv",
                 header = T,
                 stringsAsFactors = F)
  
  somat<-merge(somat,codonUsage[,c(1,3,4)],by= "ref")
  if(mutT == 1){#only mutT
    somat<-merge(somat,MutT,by=c("ref","mut"))  
  }else if(mutT == 2){#exclude mutT
    MutT$del<-1
    somat<-merge(somat,MutT,by=c("ref","mut"),all=T) 
    somat<-subset(x = somat,
                  is.na(del ))
    somat$del<-NULL
  }
  titulo<-paste("Frequency of codon mutation",ifelse(normalize,"- normalized","- without normalization"))
  if(mutT == 1){#only mutT
    titulo<-paste(titulo,"- Only mutT")
  }else if(mutT == 2){#exclude mutT
    titulo<-paste(titulo,"- mutT excluded")
  }
  p<-somat
  if(normalize){
    if(normType == "perK"){
      somat[,3:124]<-somat[,3:124]/somat$factorPerK
    }else if(normType == "absolute"){
      somat[,3:124]<-somat[,3:124]/somat$factorAbs
    }else{
      cat('Invalid normType parameter. Use "perK" or "absolute".\n')
      return()
    }
  }
  somat<-somat[order(somat$factorPerK),]
  library(plotly)
  quant=nrow(codonUsage)
  
  colors<-rainbow(n=quant,start = 0.2)
  legenda<-data.frame(x=rep(0,quant),
                      usage=codonUsage$count/1000,
                      color=colors)
  ticktext<-codonUsage$ref
  tickvals<-c(0:(length(ticktext)-1))
  
  p<-plot_ly()
  if(type =="line"){
    p<-add_trace(p=p,
                 type="scatter",
                 mode="markers",
                 data = legenda,
                 x=~x,  
                 y=~usage/100,  
                 color = ~usage,
                 colors = colors,
                 visible = "legendonly")%>%
      layout(xaxis = list(title = "k Generatios"),
             yaxis = list(title = "Frequency"),
             title = titulo)
  }else if(type == "box"){
    p<-add_trace(p=p,
                 type="scatter",
                 mode="markers",
                 data = legenda,
                 x=~x,  
                 y=~usage/100,  
                 color = ~usage,
                 colors = colors,
                 visible = "legendonly")%>%
      layout(xaxis = list(title = "Codons",
                          tickvals = tickvals,
                          ticktext = ticktext,
                          tickmode = "array" ),
             yaxis = list(title = "Frequency"),
             title = titulo)
  }
  
  colorIdx=1
  i="GGC"
  for(i in codonUsage$ref){
    #cat(i,"\n")
    q1<-somat[somat$ref == i,]
    qt1<-data.frame(x=seq(0,60.5,0.5),
                    y=apply(X = q1[,3:124],MARGIN = 2, FUN = sum))
    if(type=="line"){
      p<-add_trace(p = p, 
                   y=qt1$y,                  
                   x=qt1$x , 
                   type="scatter", 
                   mode="lines",
                   line=list(color = colors[colorIdx], width = 1) )%>%
        layout(showlegend = FALSE)
    }else if(type=="box") {
      p<-add_trace(p = p, 
                   y=qt1$y,                  
                   type="box", 
                   name = i,
                   line=list(color = colors[colorIdx], width = 1) )%>%
        layout(showlegend = FALSE)
    }
    colorIdx=colorIdx+1
  }
  if(save){
    Sys.setenv("plotly_username"="cfreis")
    # Sys.setenv("plotly_api_key"="rNieSbMA7oOGV6O7IiWm")
    # api_create(p, filename = paste0(
    #                                 sub("bSum","",sub("[.]csv","",file)),
    #                                 type,
    #                                 ifelse(normalize,"Norm","")))
    # Sys.setenv('MAPBOX_TOKEN' = 'qlfts2ohpw')
    # 
    # orca(p, paste0(dirFig,
    #                sub("bSum","",sub("[.]csv","",file)),
    #                type,
    #                ifelse(normalize,"Norm",""),
    #                ".png"))
  }
  p
}

timeline<-function(file="AuxFiles/bSumHighLow.csv",
                   normalize = F,
                   type = "box",
                   save = F,
                   normType = "perK",
                   mutT = 0,
                   workdir = "/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
  dirFig<-"/home/clovis/Doutorado/Projetos/Ecoli60/figures"
  setwd(workdir)
  somat<-read.csv(file = file,
                  stringsAsFactors = F)
  somat<-na.exclude(somat)
  name<-gsub("AuxFiles/bSum","",file)
  name<-gsub(".csv","",name)
  
  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  codonUsage<-codonUsage[order(codonUsage$factorPerK),]

  deltaW<-read.csv("AuxFiles/cDeltaW.csv",
                   header = T,
                   stringsAsFactors = F)
  MutT<-read.csv(file = "AuxFiles/dMutTmutations.csv",
                 header = T,
                 stringsAsFactors = F)
  
  somat<-merge(somat,codonUsage[,c(1,3,4)],by= "ref")
  if(mutT == 1){#only mutT
    somat<-merge(somat,MutT,by=c("ref","mut"))  
  }else if(mutT == 2){#exclude mutT
    MutT$del<-1
    somat<-merge(somat,MutT,by=c("ref","mut"),all=T) 
    somat<-subset(x = somat,
                  is.na(del ))
    somat$del<-NULL
  }
  # if(mutT == 1){#only mutT
  #   titulo<-paste(titulo,"- Only mutT")
  # }else if(mutT == 2){#exclude mutT
  #   titulo<-paste(titulo,"- mutT excluded")
  # }
  p<-somat
  if(normalize){
    if(normType == "perK"){
      somat[,3:124]<-somat[,3:124]/somat$factorPerK
    }else if(normType == "absolute"){
      somat[,3:124]<-somat[,3:124]/somat$factorAbs
    }else{
      cat('Invalid normType parameter. Use "perK" or "absolute".\n')
      return()
    }
  }
  somat<-somat[order(somat$factorPerK),]
  library(plotly)
  quant=nrow(codonUsage)
  
  colors<-rainbow(n=quant,start = 0.2)
  names(colors)<-1:quant
  fills<-alpha(colors,0.5)
  ticktext<-codonUsage$ref
  tickvals<-c(0:(length(ticktext)-1))
  colorIdx=1
  posX<-1:length(colors)
  i="CTG"
  qt<-data.frame(x=numeric(),
                 y=numeric(),
                 codon=character(),
                 color=integer(),
                 ordem=integer())
    for(i in codonUsage$ref){
      #cat(i,"\n")
      q1<-somat[somat$ref == i,]
      qtTmp<-data.frame(x=seq(0,60.5,0.5),
                      y=apply(X = q1[,3:124],MARGIN = 2, FUN = sum))
      qtTmp$codon<-i
      qtTmp$color<-colors[colorIdx]
      qtTmp$ordem<-colorIdx
      qt<-rbind(qt,qtTmp)
      colorIdx=colorIdx+1
    }
  titulo<-paste(name,ifelse(normalize,"- normalized","- without normalization"))
  
  p<-ggplot()+theme_bw()+theme(legend.position = "none")+
    ggtitle(titulo)+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)
  if(type=="line"){
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
      legenda<-
        ggplot(mydf) +
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              #axis.text.y=element_blank(),
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
                      fill=cor))+
        scale_fill_identity()+
        scale_y_continuous(name = "Codon Usage",
                           breaks = ticks$pos,
                           labels = ticks$val)
        
      
      library(grid)
      library(gridExtra)
      txt<-textGrob("Codon Usage\n (x1000)",
                    gp=gpar(fontsize=10))
      
      lay <- rbind(c(1,1,1,1,1,1,1,1,3),
                   c(1,1,1,1,1,1,1,1,2),
                   c(1,1,1,1,1,1,1,1,2),
                   c(1,1,1,1,1,1,1,1,NA))
      g<-grid.arrange(p,legenda,txt, layout_matrix = lay)
      
      if(save){
        ggsave(filename = paste0("../figures/enrich",type,rankStd,
                                 name,
                                 Dw,".pdf"), 
               plot = g, 
               device = "pdf", 
               path = workdir,
               scale = 1.5, 
               #width = 6.72, height = 2.98, units = "in",
               dpi = 300)
      }
      
    return(g)
}

correlationOld<-function(file="AuxFiles/bSumHigh.csv",
                      quant = 7, 
                      normalize = T,
                      type = "box",
                      save = F,
                      normType = "perK",
                      mutT = 0,
                      workdir = "/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
  dirFig<-"/home/clovis/Doutorado/Projetos/Ecoli60/figures"
  setwd(workdir)
  somat<-read.csv(file = file,
                  stringsAsFactors = F)
  somat<-na.exclude(somat)
  
  codonUsage<-read.csv("AuxFiles/dCodonUsage.csv",
                       header = T,
                       stringsAsFactors = F)
  #usado no cálculo da frqeuencia
  # allFreq<-sum(codonUsage$count)
  # codonUsage$freqPerK<- codonUsage$count*1000/allFreq
  # colnames(codonUsage)<-c("ref","count","factorAbs","factor")
  codonUsage<-codonUsage[order(codonUsage$factorPerK),]
  # write.csv(x=codonUsage,file="dCodonUsage.csv",
  #           row.names = F)
  
  
  deltaW<-read.csv("AuxFiles/cDeltaW.csv",
                   header = T,
                   stringsAsFactors = F)
  MutT<-read.csv(file = "AuxFiles/dMutTmutations.csv",
                 header = T,
                 stringsAsFactors = F)
  
  somat<-merge(somat,codonUsage[,c(1,3,4)],by= "ref")
  if(mutT == 1){#only mutT
    somat<-merge(somat,MutT,by=c("ref","mut"))  
  }else if(mutT == 2){#exclude mutT
    MutT$del<-1
    somat<-merge(somat,MutT,by=c("ref","mut"),all=T) 
    somat<-subset(x = somat,
                  is.na(del ))
    somat$del<-NULL
  }
  titulo<-paste("Frequency of codon mutation",ifelse(normalize,"- normalized","- without normalization"))
  if(mutT == 1){#only mutT
    titulo<-paste(titulo,"- Only mutT")
  }else if(mutT == 2){#exclude mutT
    titulo<-paste(titulo,"- mutT excluded")
  }
  p<-somat
  if(normalize){
    if(normType == "perK"){
      somat[,3:124]<-somat[,3:124]/somat$factorPerK
    }else if(normType == "absolute"){
      somat[,3:124]<-somat[,3:124]/somat$factorAbs
    }else{
      cat('Invalid normType parameter. Use "perK" or "absolute".\n')
      return()
    }
  }
  somat<-somat[order(somat$factorPerK),]
  quant=nrow(codonUsage)
  
  colors<-rainbow(n=quant,start = 0.2)
  legenda<-data.frame(x=rep(0,quant),
                      usage=codonUsage$count/1000,
                      color=colors)
  ticktext<-codonUsage$ref
  tickvals<-c(0:(length(ticktext)-1))
  
  colorIdx=1
  i="GGC"
  total<-data.frame(ref=character(),
                    freq=numeric(),
                    usage=numeric(),
                    stringsAsFactors = F)
  for(i in codonUsage$ref){
    #cat(i,"\n")
    total[nrow(total)+1,1]<-i
    q1<-somat[somat$ref == i,]
    qt1<-data.frame(x=seq(0,60.5,0.5),
                    y=apply(X = q1[,3:124],MARGIN = 2, FUN = sum))
    total[nrow(total),2]<-sum(qt1$y)
    total[nrow(total),3]<-codonUsage$count[codonUsage$ref==i]
  }
  lm<-lm(total$freq~total$usage)
  corr<-cor.test(total$freq,total$usage)
  
  g<-ggplot()+theme_bw()+
    xlab("Codon Usage")+
    ylab(paste("Mutation Frequency (",
               ifelse(normalize,"Normalized)","Absolute)")))+
    geom_point(data = total,
               aes(x=usage,y=freq), 
               col="red", 
               pch=1)+
    geom_abline(aes(linetype = "line"),slope = lm$coefficients[2], 
                intercept = lm$coefficients[1],
                lty=2,
                color="blue", show.legend = T)+
    geom_text(aes(x=max(total$usage)*0.9, y=max(total$freq)*0.9,
                  label=paste("Correlation:",
                              round(corr$estimate,3),
                              "\npvalue:",
                              format(corr$p.value,digits = 3, scientific=T))))+
    scale_colour_manual(values="blue")
  
  print(g)
  if(save){
    ggsave(filename = paste0("../figures/Corr",
                             ifelse(normalize,"Norm","NotNorm"),".pdf"), 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 1.5, 
           #width = 6.72, height = 2.98, units = "in",
           dpi = 300)
  }
}

correlation <- function(population="HighLow",
                      save = F,
                      normalize = F,
                      workdir = "/home/clovis/Dropbox/Ecoli60/data_files/",
                      type = "sum",
                      TsTv="all"){
  if(!type %in% c("mean","sum")){
    stop('Parameter "type" must be "mean" or "sum"' )
  }
  if(!TsTv %in% c("Ts","Tv","all")){
    stop('Parameter "TsTv" must be "Ts", "Tv" or "all"' )
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

    #acrescenta fator
  mutAllelSyn<-merge(mutAllelSyn,codonUsage[,c(1:3)],by=c("ref"))
  colnames(mutAllelSyn)[137]<-"factor"
  
  #normaliza valores
  if(normalize){
    tmp<-mutAllelSyn[,7:128]/mutAllelSyn$factor
    tmp[is.na(tmp)]<- 0
    mutAllelNorm<-cbind(mutAllelSyn[,1:6],tmp,mutAllelSyn[129:137])
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
           scale = 1.5, 
           #width = 6.72, height = 2.98, units = "in",
           dpi = 300)
  }
  return(g)
  
}

freqPerDw <- function(population="HighLow",
                      save = F,
                      workdir = "/home/clovis/Dropbox/Ecoli60/data_files/",
                      type = "sum",
                      Dw = "Tai",
                      TsTv="all",
                      top=86){
  if(!type %in% c("mean","sum")){
    stop('Parameter "type" must be "mean" or "sum"' )
  }
  if(!TsTv %in% c("Ts","Tv","all")){
    stop('Parameter "TsTv" must be "Ts", "Tv" or "all"' )
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
  tmp<-mutAllelSyn[,7:128]/mutAllelSyn$factor
  tmp[is.na(tmp)]<- 0
  mutAllelNorm<-cbind(mutAllelSyn[,1:6],tmp,mutAllelSyn[129:136])
  
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
    sumario$MT[round(sumario$dw,5)%in%round(mutT$dw,5)]<-1
    
    g<-ggplot()+theme_bw()+
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
    l<-legend(0, 1, legend=c(" Tv", " Ts"," Zeros","MutT"),
              col=c("red", "darkmagenta","green","darkgoldenrod2"),
              pch = c(20,18,2,8), 
              cex=0.7,
              box.lwd = 0)
    leg <- recordPlot()
    
    
  }else{
    g<-ggplot()+theme_bw()+
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
    l<-legend(0, 1, legend=c(" Tv", " Ts"," Zeros"),
              col=c("red", "darkmagenta","green"),
              pch = c(20,18,2), 
              cex=0.7,
              box.lwd = 0)
    leg <- recordPlot()
    
  }
  g<-g+geom_text(aes(y=(max(sumario$freq)*2/3), 
                     x=(min(deltaW$dw)), 
                     label = paste(nrow(zeros),"zeros"),
                     hjust = "left"))
  g2<-plot_grid(g, NULL,leg,
               nrow=1,
               ncol = 3,
               rel_widths=(c(9.5,0.5,1)))
  print(g2)
  if(save){
    ggsave(filename = paste0("../figures/MutPerDw",population,".pdf"), 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 1.5, 
           #width = 6.72, height = 2.98, units = "in",
           dpi = 300)
  }
}

genesExpression<-function(baseDir="/home/clovis/Doutorado/Projetos/Ecoli60/" ){
  
  workdir = paste0(baseDir,"data_files/dataExpresion/LeGac2012/")
  
  setwd(workdir)
  
  
  
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
                       workdir = "/home/clovis/Dropbox/Ecoli60/data_files/"){
  source("/home/clovis/Doutorado/Projetos/Ecoli60/bin/allFunctions.R")
  mutAllelAll<-read.csv(paste0("base/a",population,"MutAllDW.csv"),
                        header = T,
                        stringsAsFactors = F)
  #só syn
  if(justSyn){
    mutAllelAll<-mutAllelAll[mutAllelAll$Annotation == "synonymous",]
  }
  expression<-genesExpression()
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
                      Dw = "Tai",
                      withLow = F,
                      workdir = "/home/clovis/Dropbox/Ecoli60/data_files/"){
  
  setwd(workdir)
  if(Dw == "Tai"){
    deltaW<-read.csv("./AuxFiles/cDeltaWTAI.csv",
                     header = T,
                     stringsAsFactors = F)
  }else if(Dw == "Cai"){
    deltaW<-read.csv("./AuxFiles/cDeltaW.csv",
                     header = T,
                     stringsAsFactors = F)
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

calcDWTAI<-function(tRNA,
                    workdir="/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
  setwd(workdir)
  deltaW<-read.csv("AuxFiles/cDeltaW.csv",
                   header = T,
                   stringsAsFactors = F)
  deltaW$dw<-0
  
  deltaW<-merge(deltaW,
                tRNA[,c("codons","w")],
                by.x ="ref",
                by.y="codons")
  colnames(deltaW)<- c("ref","mut","dw","wRef")
  deltaW<-merge(deltaW,
                tRNA[,c("codons","w")],
                by.x ="mut",
                by.y="codons")
  colnames(deltaW)<- c("mut","ref","dw","wRef", "wMut")
  deltaW$dw<-log10(deltaW$wMut/deltaW$wRef)
  deltaW<-deltaW[,c("ref","mut","dw")]
  write.csv(deltaW,
            file = "AuxFiles/cDeltaWTAI.csv",
            quote = F,
            row.names = F)
}

corrTaiCai<- function(save = F,
                      workdir="/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
  setwd(workdir)
  deltaWCAI<-read.csv("AuxFiles/cDeltaW.csv",
                      header = T,
                      stringsAsFactors = F)
  deltaWTAI<-read.csv("AuxFiles/cDeltaWTAI.csv",
                      header = T,
                      stringsAsFactors = F)
  library(ggplot2)
  deltas<-merge(deltaWCAI,deltaWTAI,
                by=c("ref","mut"))
  colnames(deltas)<-c("ref","mut","CAI","TAI")
  lm<-lm(deltas$CAI~deltas$TAI)
  corr<-cor.test(deltas$CAI,deltas$TAI)
  maxX<-max(deltas$TAI)
  maxY<-max(deltas$CAI)
  library(ggplot2)
  g<-ggplot()+theme_bw()+
    xlab(bquote(Delta~"w TAI"))+
    ylab(bquote(Delta~"w CAI"))+
    geom_point(data = deltas,
               aes(x=TAI,y=CAI), 
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
  print(g)
  if(save){
    ggsave(filename = "../figures/corrDwCaiTai.pdf", 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 1.5,
           #width = 3*maxX, 
           #height = 3*maxY, 
           #units = "in",
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

totalize_tRNA<-function(baseDir="/home/clovis/Doutorado/Projetos/Ecoli60/data_files/TAIFiles/"){
  setwd(baseDir)
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
  
  tRNA<-read.csv(file="REL606tRNAPre.txt",
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




createMonteCarloDF<-function(workdir = "/home/clovis/Dropbox/Ecoli60/data_files/"){
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
                     workdir = "/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
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
    deltaW<-read.csv("./AuxFiles/cDeltaW.csv",
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
           scale = 1.5,
           #width = 3*maxX, 
           #height = 3*maxY, 
           #units = "in",
           dpi = 300)
  }
}

enrichZeros <- function(save = F,
                      quant=5,
                      population = "HighLow",
                      Dw = "Tai",
                      workdir = "/home/clovis/Dropbox/Ecoli60/data_files/"){
  rm(list = ls())
  save = F
  quant=5
  population = "HighLow"
  Dw = "Tai"
  workdir = "/home/clovis/Dropbox/Ecoli60/data_files/"
  setwd(workdir)
  if(Dw == "Tai"){
    deltaW<-read.csv("./AuxFiles/cDeltaWTAI.csv",
                     header = T,
                     stringsAsFactors = F)
  }else if(Dw == "Cai"){
    deltaW<-read.csv("./AuxFiles/cDeltaW.csv",
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
                       population="Low",
                       Dw = "Tai",
                       workdir = "/home/clovis/Dropbox/Ecoli60/data_files/"){
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
                   Dw="Tai",
                   withLow = T)
  TsTv<-read.csv("AuxFiles/TsTv.csv",
                 header = T,
                 stringsAsFactors = F)
  if(population == "all"){
    populations<-c("High","MutT","HighLow","Low")
  }else{
    populations<-population
  }
  for(population in populations){
    zPop<-zeros[zeros[[population]]=="X",c("ref","mut","dw","level")]
    zPop<-merge(TsTv,zPop, by=c("ref","mut"))
    dw<-data.frame(Positive=nrow(zPop[zPop$dw<0,]),
                   Negative=nrow(zPop[zPop$dw>=0,]),
                   Total = nrow(zPop))
    cat("\n\nZeros from Population", population,"\nDeltas:\n")
    print(dw, row.names = F)
    level<-as.data.frame(table(zPop[,c("level")]),
                         stringsAsFactors = F)
    colnames(level)<-c("Levels  ", "Freq")
    print(level, row.names = F)
    typeT<-as.data.frame(table(zPop[,c("type")]),
                      stringsAsFactors = F)
    colnames(typeT)<-c("Mut Type", "Freq")
    print(typeT, row.names = F)
    
  }
}


countTsTv <- function(save = F,
                      #TsTv="Tv",
                      Dw = "Tai",
                      type = "sum",
                      population = "HighLow",
                      workdir = "/home/clovis/Dropbox/Ecoli60/data_files/"){
  
  setwd(workdir)
  if(Dw == "Tai"){
    deltaW<-read.csv("./AuxFiles/cDeltaWTAI.csv",
                     header = T,
                     stringsAsFactors = F)
  }else if(Dw == "Cai"){
    deltaW<-read.csv("./AuxFiles/cDeltaW.csv",
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
                     workdir = "/home/clovis/Dropbox/Ecoli60/data_files/"){
  
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

comparePopulations<-function(save = F,
                             workdir= "/home/clovis/Dropbox/Ecoli60/data_files/"){
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
    ylab("Number of Mutations")
  print(g)
  if(save){
    ggsave(filename = "../figures/compPopulations.pdf", 
           plot = g, 
           device = "pdf", 
           path = workdir,
           scale = 1.5,
           #width = 3*maxX, 
           #height = 3*maxY, 
           #units = "in",
           dpi = 300)
  }
  return(g)
}

figureS01<-function(workdir= "/home/clovis/Dropbox/Ecoli60/data_files/"){
  g2<-comparePopulations(save = F,
                   workdir= workdir)
  ggsave(filename = "../figures/FigS01.pdf", 
         plot = g2, 
         device = "pdf", 
         path = workdir,
         scale = 1.5,
         #width = 3*maxX, 
         #height = 3*maxY, 
         #units = "in",
         dpi = 300)
  
}

figureS02<-function(workdir= "/home/clovis/Dropbox/Ecoli60/data_files/"){
  g<-list()
  index<-1
  for(pop in c("High", "MutT")){
    for(t in c("Ts","Tv")){
      g[[index]]<-correlation(population=pop,
                              save = F,
                              normalize = F,
                              workdir = workdir,
                              type = "sum",
                              TsTv=t)
      index<-index+1
    }
    g[[index]]<-correlation(population="Low",
                            save = F,
                            normalize = F,
                            workdir = workdir,
                            type = "sum",
                            TsTv="all")
  }
  library(cowplot)
  g2<-plot_grid(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],
                labels=c('A', 'B', 'C', 'D','E'),
                nrow=3,
                ncol = 2)
  
  ggsave(filename = "../figures/FigS02.pdf", 
         plot = g2, 
         device = "pdf", 
         path = workdir,
         scale = 2.5,
         #width = 3*maxX, 
         #height = 3*maxY, 
         #units = "in",
         dpi = 300)
  
}

figureS03<-function(workdir= "/home/clovis/Dropbox/Ecoli60/data_files/"){  g<-list()
  index<-1
  for(pop in c("High", "MutT")){
    for(t in c("Ts","Tv")){
      g[[index]]<-correlation(population=pop,
                              save = F,
                              normalize = T,
                              workdir = workdir,
                              type = "sum",
                              TsTv=t)
      index<-index+1
    }
    g[[index]]<-correlation(population="Low",
                            save = F,
                            normalize = T,
                            workdir = workdir,
                            type = "sum",
                            TsTv="all")
    
  }
  library(cowplot)
  g2<-plot_grid(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],
                labels=c('A', 'B', 'C', 'D','E'),
                nrow=3,
                ncol = 2)
  
  ggsave(filename = "../figures/FigS03.pdf", 
         plot = g2, 
         device = "pdf", 
         path = workdir,
         scale = 2.5,
         #width = 3*maxX, 
         #height = 3*maxY, 
         #units = "in",
         dpi = 300)
}

figureS04e5<-function(type="box",
  workdir= "/home/clovis/Dropbox/Ecoli60/data_files/"){
  setwd(workdir)
  g<-list()
  index<-1
  pop="High"
  for(pop in c("High","Low", "MutT")){
    file<-paste0("AuxFiles/bSum",pop,".csv")  
    g[[index]]<-timeline(file=file,
                         normalize = F, 
                         type = type,
                         save = F, #save option is not working here
                         normType = "perK",
                         mutT = 0,
                         workdir = workdir)
    index<-index+1
    g[[index]]<-timeline(file=file,
                         normalize = T, 
                         type = type,
                         save = F, #save option is not working here
                         normType = "perK",
                         mutT = 0,
                         workdir = workdir)
    index<-index+1
  }
  library(cowplot)
  g2<-plot_grid(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],g[[6]],
                labels=c('A', 'B', 'C', 'D','E','F'),
                nrow=3,
                ncol = 2)
  if(type == "box"){
    filename=paste0("../figures/FigS05.pdf")
  }else if(type == "line"){
    filename=paste0("../figures/FigS04.pdf")
  }
  ggsave(filename = filename, 
         plot = g2, 
         device = "pdf", 
         path = workdir,
         scale = 2.5,
         #width = 3*maxX, 
         #height = 3*maxY, 
         #units = "in",
         dpi = 300)
  g2
}

figure01<-function(type="box",
                      workdir= "/home/clovis/Dropbox/Ecoli60/data_files/"){
  setwd(workdir)
  g<-list()
  index<-1
  pop="HighLow"
  file<-paste0("AuxFiles/bSum",pop,".csv")  
  
  for(type in c("line","box")){
    for (normalize in c(F,T)) {
      g[[index]]<-timeline(file=file,
                           normalize = normalize, 
                           type = type,
                           save = F, #save option is not working here
                           normType = "perK",
                           mutT = 0,
                           workdir = workdir)
      index<-index+1
    }
  }
  library(cowplot)
  g2<-plot_grid(g[[1]],g[[2]],g[[3]],g[[4]],
                labels=c('A', 'B', 'C', 'D'),
                nrow=2,
                ncol = 2)
  filename=paste0("../figures/Fig01.pdf")
  ggsave(filename = filename, 
         plot = g2, 
         device = "pdf", 
         path = workdir,
         scale = 2.5,
         #width = 3*maxX, 
         #height = 3*maxY, 
         #units = "in",
         dpi = 300)
  g2
}


listGenes<-function(top=86,
                    workdir){
  setwd(workdir)
  expression<-read.csv("AuxFiles/highlyEGLeGac",
                                   header = F,
                                   stringsAsFactors = F)
  expression<-expression[1:top,]
  membrane<-read.csv("AuxFiles/membraneGenes.txt",
                      header = F,
                      stringsAsFactors = F)
  ribossomal<-c("rrsA","rpsA","rpsB","rpsC","rpsD","rpsE","rpsF","rpsG","rpsH","rpsI","rpsJ","rpsK","rpsL","rpsM","rpsN","rpsO","rpsP","rpsQ","rpsR","rpsS","rpsT","rpsU","Sra","rrlA","rrfA","rplA","rplB","rplC","rplD","rplE","rplF","rplJ","rplL","rplI","rplK","rplM","rplN","rplO","rplP","rplQ","rplR","rplS","rplT","rplU","rplV","rplW","rplX","rplY","rpmA","rpmB","rpmC","rpmD","rpmE","rpmF","rpmG","rpmH","rpmI","rpmJ")
  r<-expression[expression%in%ribossomal]
  membr<-expression[expression%in%membrane$V1]
  
}