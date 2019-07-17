###################################################
# Create the dataframe with montecarlo mutations 
###################################################
rm(list = ls())

baseDir="Your Base dir"
#baseDir="/home/clovis/Doutorado/Projetos/Ecoli60/"

binDir=paste0(baseDir,"bin/")
workdir=paste0(baseDir,"data_files/")
setwd(workdir)
setwd(binDir)
source(paste0(binDir,"allFunctions.R"))

createMonteCarloDF(workdir = workdir)
  