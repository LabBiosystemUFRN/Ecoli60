rm(list=ls())

source(paste0(binDir,"allFunctions.R"))

#convertToComplement("TGA")


#require("tAI")
codonOrder<-c("TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA","TAG","TGT","TGC","TGA","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG")
antiCodons<- sapply(codonOrder, convertToComplement )

codonOrder[c(11, 12, 15, 36)]
antiCodons[c(11, 12, 15, 36)]

totalize_tRNA(baseDir="/home/clovis/Doutorado/Projetos/Ecoli60/data_files/TAIFiles/")
  

setwd("/home/clovis/Doutorado/Projetos/Ecoli60/data_files/TAIFiles/")
tRNA<-read.table("REL606tRNA.txt",
                 header = T,
                 stringsAsFactors = F)
#eco.trna <- tRNA[,2]
#teste2<-data.frame(codon= codonOrder,val=eco.trna)
#write.table(teste2,file = "trnammmn2.txt",sep="\t", row.names = F)

#tRNA<-teste2
tRNA <- get.ws(tRNA=tRNA, sking=1)

tRNA2<-read.table(file = "/home/clovis/Doutorado/Projetos/Ecoli60/data_files/TAIFiles/output_wi_file.txt",
                  header = T,
                  stringsAsFactors = F)
colnames(tRNA2)<-c("codons","w")

calcDWTAI(tRNA = tRNA2)

ecow<-data.frame(codon= codons[c(1:32,34:61)],w=eco.ws)

eco.m <- matrix(scan("ecolik12.m"), ncol=61, byrow=TRUE)
eco.m <- eco.m[,-33]
eco.tai <- get.tai(eco.m, eco.ws)

w=eco.ws
x=eco.m
w = log(w)
n = apply(x, 1, "*", w)
n = t(n)
n = apply(n, 1, sum)
L = apply(x, 1, sum)
tAI = exp(n/L)


hist(eco.tai)

df <- read.table("ecolik12.w", header=TRUE)
plot(eco.tai ~ df$Nc)

ecolik12$m