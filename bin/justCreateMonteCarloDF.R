#################
# main processing 
#################
rm(list = ls())

baseDir="/home/clovis/Dropbox/Ecoli60/"
#baseDir="/home/clovis/Doutorado/Projetos/Ecoli60/"

binDir=paste0(baseDir,"bin/")
workdir=paste0(baseDir,"data_files/")
setwd(workdir)
setwd(binDir)
source(paste0(binDir,"allFunctions.R"))

createMonteCarloDF(workdir = "/home/clovis/Dropbox/Ecoli60/data_files/")
  