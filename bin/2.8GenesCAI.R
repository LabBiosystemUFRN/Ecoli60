rm(list = ls())

combination<-function(decNumber) {
  bsum<-vector()
  bexp<-1
  #decNumber=4
  #numero máximo de casas binárias
  maxCont<-ceiling(log2(decNumber+1))
  cont = 1
  while (decNumber > 0) {
    digit<-decNumber %% 2
    decNumber<-floor(decNumber / 2)
    if(digit == 1){
      bsum<-c(bsum,cont)
    }
    cont<-cont+1
  }
  return(bsum[order(bsum)])
}


binary<-function(p_number) {
  bsum<-0
  bexp<-1
  while (p_number > 0) {
    digit<-p_number %% 2
    p_number<-floor(p_number / 2)
    bsum<-bsum + digit * bexp
    bexp<-bexp * 10
  }
  return(bsum)
}

#codons<-codonsGene
calcCAI <- function(codons){
  #calculo CAI --
  #codons must be a dataframe with, at least, 3 columns: codon, w values and count
  #sem perda
  #Somatório dos logs de w iguais
  codons$cai<-log(codons$w)*codons$count  
  lnw<-sum(codons$cai)
  L<-sum(codons$count)
  cai<-exp(lnw/L)
  return(cai)
}

type="High"
pval=0.01
quant=7
normalize = T
plot = T
rank=100
tail="L"
TsTv = "all"
fix = "all"
graph = "box"

genesCAI<- function(type="High",
                  pval=0.01,
                  quant=7, 
                  normalize = F, 
                  rank=100, 
                  tail="L",
                  TsTv ,
                  fix = T,
                  graph = "box",
                  workdir="/home/clovis/Doutorado/Projetos/Ecoli60/data_files/"){
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
  codonW<-read.csv("cCodonsW.csv",
                       header = T,
                       stringsAsFactors = F)
  codonGenes<-read.csv("cCodonsGenes.csv",
                   header = T,
                   stringsAsFactors = F)
  codonW<-codonW[,c(1,3)]
  colnames(codonW)<-c("ref","w")
  
  mutAllel<-read.csv(paste0("a",type,"MutSynFreqDW.csv"),
                     header = T,
                     stringsAsFactors = F)
  
  mutAllelAll<-read.csv(paste0("a",type,"MutAllDW.csv"),
                     header = T,
                     stringsAsFactors = F)
  nonGenes<-unique(mutAllelAll$Gene[mutAllelAll$Annotation!="synonymous"])
  
  mutAllel<-mutAllel[!mutAllel$Gene%in%nonGenes,]
  #unique(mutAllel$Gene)
  # hiGenes<-read.csv("/home/clovis/Doutorado/Projetos/Ecoli60/additional_data/highlyEGLeGac",
  #                   header = F,
  #                   stringsAsFactors = F)
  # hiGenes<-data.frame(gene=hiGenes[-(grep(pattern = "#",hiGenes$V1)),])
  
  #remove stop codons
  mutAllel<-mutAllel[!(mutAllel$ref == "TAA"| mutAllel$ref == "TGA"),]
  mutAllelAll<-mutAllelAll[!(mutAllelAll$ref == "TAA"| mutAllelAll$ref == "TGA"),]
  #filtra Ts ou Tv
  if(TsTv != "all"){
    mutAllel<-subset(mutAllel, TsTv == vTsTv)
    mutAllelAll<-subset(mutAllelAll, TsTv == vTsTv)
  }
  #fix=T
  if(fix !="all"){
    if(fix){
      mutAllel<-mutAllel[!is.na(mutAllel$fixed),]
      mutAllelAll<-mutAllelAll[!is.na(mutAllelAll$fixed),]
    }else{
      mutAllel<-mutAllel[is.na(mutAllel$fixed),]
      mutAllelAll<-mutAllelAll[is.na(mutAllelAll$fixed),]
    }
  }
  
  mutFiltered<-merge(mutAllel[c(1:5,130:134)],codonW,by="ref")
  mutFiltered<-merge(mutFiltered,codonW,by.x="mut",by.y = "ref")
  colnames(mutFiltered)<-c("mut","ref","Position","Gene","Allele","pop","t0",
                           "fixed","transit","TsTv","wRef","wMut")
  
  genes<-unique(mutAllel$Gene)
  
  pops<-unique(mutAllel$pop)
  
  library(dplyr)
  #Teste quantidade de mutações máximas em um gene
  n<-mutAllel %>%
    group_by(Gene,pop) %>%
    summarize(count = n())
  maxCombination<-max(n$count)
  rm(n)
  #maximo de mutações são 6. Ajustar se este numero mudar
  
  
  combination(10)
  binary(62)
  
  combinations<-list()
  for(i in 1:(2^maxCombination)){
    combinations[[i]]<-combination(i)
  }

  resultPop<-list()
  
  pop="p3"
  countPop=1
  #processa pops e genes ----
  for(pop in pops){
    gene="acrB"
    
    resultGene<-data.frame(matrix(nrow = 0,ncol = 122))
    colnames(resultGene)<-c("gene",seq(0,60000,500))
    #cai = rep(0,121)
    
    # resultGene<-list()
    countGene=1
    for(gene in genes){
      mutGene<-mutFiltered[mutFiltered$pop == pop &
                         mutFiltered$Gene == gene,]
      if(nrow(mutGene)==0){
        next
      }
      codonsGene<-codonGenes[codonGenes$gene == gene,]
      codonsGene<-merge(codonsGene,
                        codonW,
                        by.x = "codon",
                        by.y = "ref")
      caiRef<-calcCAI(codonsGene)
      timeLine<-as.data.frame(matrix(data = 0,
                                     ncol = nrow(mutGene),
                                     nrow = (60000/500)+1))
      
      timeLine$time<-seq(0,60000,500)
      i=1
      #cria tabela com sobreposição das mutações
      for (i in 1:nrow(mutGene)) {
          val<-2^(i-1)
          maxTime<-ifelse(is.na(mutGene$fixed[i]),
                          (mutGene$t0[i]+mutGene$transit[i]),
                          65000)
          timeLine[timeLine$time>=mutGene$t0[i]&
                     timeLine$time<=maxTime,
                   i]<-val
      }
      if(nrow(mutGene)==1){
        timeLine$tot<-timeLine[,1]
      }else{
        timeLine$tot<-apply(timeLine[,1:nrow(mutGene)],
                          MARGIN = 1,
                          sum)
      }
      timeLine$cai<-0
      #verifica tipos de sobreposição diferentes
      sobrep<-unique(timeLine$tot[timeLine$tot != 0])
      mutations<-data.frame(ref=character(), mut = character())
      indexSob = 1
      for (indexSob in sobrep) {
        codonsTmp<-codonsGene
        antes<-sum(codonsTmp$count)
        cat(indexSob,binary(indexSob),combinations[[indexSob]],"\n")
        for (indexMut in combinations[[indexSob]]) {
          cat(indexMut,mutGene$ref[indexMut],mutGene$mut[indexMut],"\n")
          codonsTmp$count[codonsTmp$codon == mutGene$ref[indexMut]]<-codonsTmp$count[codonsTmp$codon == mutGene$ref[indexMut]]-1
          codonsTmp$count[codonsTmp$codon == mutGene$mut[indexMut]]<-codonsTmp$count[codonsTmp$codon == mutGene$mut[indexMut]]+1
        }
        depois<-sum(codonsTmp$count)
        if(antes-depois!=0){
          cat("Deu merda no gene:",gene )
        }
        caiMut<-calcCAI(codonsTmp)
        timeLine$cai[timeLine$tot == indexSob]<- log10(caiMut/caiRef)
        # timeLine$cai[timeLine$tot == indexSob]<- caiMut-caiRef
        
        #mutations<-
      }
      resultGene[countGene,1]<-gene
      resultGene[countGene,2:122]<-t(timeLine$cai)
      # resultGene[[countGene]]<-list(gene,timeLine[,c("time","cai")])
      countGene<-countGene+1
    }
    
    resultPop[[countPop]]<-list(pop,resultGene)
    countPop<-countPop+1
  }
  return(resultPop)
}

#listaCAI<-p
plotCAI<-function(listaCAI, 
                  pop = "all",
                  graph = "box"){
  if(!graph%in%c("box","line")){
    cat('Parameter graph must be "box" or "line"\n')
    return()
  }
  if(class(listaCAI) != "list"){
    cat("Parameter listaCAI must be an list\n")
    return()
  }
  library(ggplot2)
  max=0
  min=0
  result<-list()
  indexPop=1
  for(indexPop in 1:length(listaCAI)){
    pop<-listaCAI[[indexPop]][[1]]
    caiTbl<-listaCAI[[indexPop]][[2]]
    colors<-rainbow(nrow(caiTbl))
    p<-ggplot()+
      theme_bw()+
      labs(title =  paste("Population",pop))+
      xlab("Generations")+
      ylab("CAI ratio (log10)")
    if(graph == "line"){
      x<-as.numeric(colnames(caiTbl[2:ncol(caiTbl)]))
      idxLine = 1
      for (idxLine in 1:nrow(caiTbl)) {
        max<-max(max,max(caiTbl[idxLine,2:ncol(caiTbl)]))
        min<-min(min,min(caiTbl[idxLine,2:ncol(caiTbl)]))
        tmp<-data.frame(x = x, 
                        y = t(caiTbl[idxLine,2:ncol(caiTbl)]),
                        stringsAsFactors = F)
        colnames(tmp)<-c("x","y")
        p<-p+geom_line(data = tmp,
                     aes(x,y),
                     col=alpha(colors[idxLine],0.2))
        
      }  
    }else{
      plotTbl<-data.frame(x=character(),
                          y=numeric())
      idx="10000"
      for(idx in c("10000","20000","30000",
                   "40000","50000","60000")){
        max<-max(max,max(caiTbl[colnames(caiTbl)==idx]))
        min<-min(min,min(caiTbl[colnames(caiTbl)==idx]))
        
        tmp<-data.frame(x=rep(idx,nrow(caiTbl)),
                        y=caiTbl[colnames(caiTbl)==idx])
        colnames(tmp)<-c("x","y")
        plotTbl<-rbind(plotTbl,tmp)
      }
      p<-p+geom_boxplot(data = plotTbl,aes(x,y))
    }
    plot(p)
    result[[indexPop]]<-p
  }
  cat("Max value:",max,"\nMin Value:",min,"\n",10^max-1,"to",10^min-1,"\n")
  return(result)
}


listaCAI<-genesCAI(type="High",
               pval=0.001,
               quant=8, 
               normalize = T, 
               rank=100,
               graph = "line",
               fix = T,
               TsTv = "all",
               workdir = "/home/clovis/Dropbox/Ecoli60/data_files/")
p<-plotCAI(listaCAI,graph = "box")
unique(listaCAI[[1]][[2]]$gene)
  
p1<-p[[1]]
#subset(x = mutFiltered, is.na(fixed),select = c(t0 ,transit))%>%apply(.,MARGIN = 1,sum)%>%min(.)
