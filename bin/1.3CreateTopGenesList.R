rm(list = ls())

args = commandArgs(trailingOnly=TRUE)

baseDir=args[1]
#baseDir="/home/clovis/Dropbox/Ecoli60/"
#source("/home/clovis/Doutorado/Projetos/Ecoli60/bin/allFunctions.R")
source(paste0(baseDir,
              "bin/allFunctions.R"))

#nGenes = commandArgs(trailingOnly=TRUE)
#nGenes = 1500
expression<-genesExpression(paste0(baseDir,"data_files/"))
nGenes = nrow(expression)
topExpressed<-expression$symbol[1:nGenes]
cat("Using top",nGenes,"genes\n",
    "Ribossomal:", (length(topExpressed[grep("rpl",topExpressed)])+
                      length(topExpressed[grep("rps",topExpressed)])+
                      length(topExpressed[grep("rpm",topExpressed)])),
    "\n")
sharp<-read.csv(paste0(baseDir,
                       "additional_data/highlyEGSharp"),
                header = F)[c(1:25,27),]

# ambos<-merge(data.frame(symbol=topExpressed,ref1=rep("t",length(topExpressed))),
#              data.frame(symbol=sharp,ref2=rep("s",length(sharp))), all=F)

write.table(topExpressed,
          file = paste0(baseDir,"data_files/AuxFiles/highlyEGLeGac"),
          quote = F,
          row.names = F, 
          col.names = F,
          sep = ",")
