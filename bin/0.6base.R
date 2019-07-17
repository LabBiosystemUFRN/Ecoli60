rm(list = ls())
#Did you change it to your base location?
args = commandArgs(trailingOnly=TRUE)

baseDir=args[1]
#baseDir="/home/clovis/temp/teste/Ecoli60/"
setwd(baseDir)

#packages required:
#cowplot, doParallel, dplyr, ecoli2.db, ggplot2, grDevices, grid, gridExtra, 
#gridGraphics, lattice, org.EcK12.eg.db, reshape2, RMySQL, scales, stats, and stringr

if (!require("BiocManager")){
  install.packages("BiocManager")
}

if (!require("cowplot")){
  install.packages("cowplot")
}

if (!require("doParallel")) {
  install.packages("doParallel")
}

if (!require("dplyr")) {
  install.packages("dplyr")
}

if (!require("ecoli2.db")) {
  BiocManager::install("ecoli2.db")
}

if (!require("ggplot2")) {
  install.packages("ggplot2")
}

if (!require("grid")) {
  install.packages("grid")
}

if (!require("gridExtra")) {
  install.packages("gridExtra")
}

if (!require("ecoli2.db")) {
  BiocManager::install("ecoli2.db")
}


if (!require("org.EcK12.eg.db")){
  BiocManager::install("org.EcK12.eg.db")
}

if (!require("stats")) {
  install.packages("stats")
}

if (!require("grDevices")) {
  install.packages("grDevices")
}

if (!require("gridGraphics")) {
  install.packages("gridGraphics")
}

if (!require("lattice")) {
  install.packages("lattice")
}

if (!require("reshape2")) {
  install.packages("reshape2")
}

if (!require("RMySQL")) {
  install.packages("RMySQL")
}

if (!require("scales")) {
  install.packages("scales")
}

if (!require("stringr")) {
  install.packages("stringr")
}
