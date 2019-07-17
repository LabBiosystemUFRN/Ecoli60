#This script correct the error found in the original code
#where TGC e TGT are not considered as Cysteine synnonymous
#line 21 of parse_file.py
#'GAT':'D', 'GAC':'D', 'TGT':'C', 'TGC':'D'
#have no mutations GAT or GAC to TGC or vice-versa
#but exist 113 mutations TGT to TGC or vice-versa
#this mutations are reclassified as synonymous
#It also add information about delta w to file

rm(list = ls())

args = commandArgs(trailingOnly=TRUE)

baseDir=args[1]

dirIn = paste0(baseDir,"data_files/AuxFiles/")
dirOut = paste0(baseDir,"data_files/base/")
setwd(baseDir)

deltaW<-read.csv(paste0(dirIn,"cDeltaW86.csv"),
                 header = T,
                 stringsAsFactors = F)
codonUsage<-read.csv(paste0(dirIn,"dCodonUsage.csv"),
                     header = T,
                     stringsAsFactors = F)
type = "High"
#dá um erro no MutT
for(type in c("High","Low","MutT")){
  mutAllel<-read.csv(paste0(dirIn,"a",type,"MutSynFreqErroTGC.csv"),
                   header = T,
                   stringsAsFactors = F)
  mutAllelAll<-read.csv(paste0(dirIn,"a",type,"MutAllFreqErroTGC.csv"),
                      header = T,
                      stringsAsFactors = F)
  
  #Transition & transversion column
  mutAllelAll$TsTv<-"Tv"
  mutAllelAll$TsTv[mutAllelAll$Allele == "A->G"|
     mutAllelAll$Allele == "G->A"|
     mutAllelAll$Allele == "T->C"|
     mutAllelAll$Allele == "C->T"]<-"Ts"

  mutAllel2<-merge(mutAllelAll,deltaW[,c(1,2,3)], by=c("ref","mut"))
  mutAllelAll<-merge(mutAllelAll,deltaW[,c(1,2,3)], 
                     by=c("ref","mut"), all.x = T)
  nrow(mutAllel2[!mutAllel2$Position%in%mutAllel$Position,1:6])
      if(nrow(mutAllel2[mutAllel2$Annotation!="synonymous",])==
         nrow(mutAllelAll[mutAllelAll$ref=="TGC"&mutAllelAll$mut=="TGT" |
                       mutAllelAll$ref=="TGT"&mutAllelAll$mut=="TGC",])){
        mutAllel2$Annotation[mutAllel2$Annotation!="synonymous"]<-"synonymous"
        mutAllelAll$Annotation[mutAllelAll$ref=="TGC"&mutAllelAll$mut=="TGT" |
                           mutAllelAll$ref=="TGT"&mutAllelAll$mut=="TGC"]<-"synonymous"
        if(nrow(mutAllel2) == 
           nrow(mutAllelAll[mutAllelAll$Annotation=="synonymous"& (
             mutAllelAll$ref!="TAA"&mutAllelAll$ref!="TGA"),])){
          cat(type," ... Ok\n")
          
          colnames(mutAllelAll)[colnames(mutAllelAll) == 'LeGac'] <- 'dw'
          
          write.csv(mutAllelAll,file = paste0(dirOut,"a",type,"MutAllDW.csv"),
                    row.names = F)

          colnames(mutAllel2)[colnames(mutAllel2) == 'LeGac'] <- 'dw'
          write.csv(mutAllel2,file = paste0(dirOut,"a",type,"MutSynFreqDW.csv"),
                    row.names = F)
          
          }else{
            cat("Error in ", type,"\n")
            next
          }
      }else{
        cat("Error in ", type,"\n")
        next
      }
}
  
