#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#           main processing 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#This script generate the figures in 
# the same order than it apears in teh text
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#inicio ----
rm(list = ls())

args = commandArgs(trailingOnly=TRUE)

baseDir=args[1]
#baseDir="/home/clovis/temp/teste/Ecoli60/"
#baseDir="/home/clovis/Doutorado/Projetos/Ecoli60/"

binDir=paste0(baseDir,"bin/")
workdir=paste0(baseDir,"data_files/")
figDir=paste0(baseDir,"figures/")
setwd(workdir)
setwd(binDir)
source(paste0(binDir,"allFunctions.R"))

#create dWTai table
calcDWTAI(type = "Max",workdir = workdir)
#create a table with all High and Low mutations
joinHighLow(workdir = workdir)
#Create a MutT list
mutationsMutT(workdir = workdir)

#Figure 01 ----
#Number of mutation per population & type
figFunc01(figName = "Fig01",
          workdir)
comparePopulations(figName = "Fig01a",
                   save = T,
                   workdir = workdir)

#Figure S02 ----
#Frequencie of mutations per codon usage A to F
figFunc02(figName = "FigS02",
          separeTsTv = T,
          workdir)

#Figure S03 ----
#Frequencie of mutations per codon usage after normalization A to F
figFunc03(figName = "FigS03Count",
          normBy="count",
          separeTsTv = T,
          workdir)
figFunc03(figName = "FigS03CUB",
          normBy="CUB",
          separeTsTv = T,
          workdir)
figFunc03(figName = "FigS03Mean",
          normBy="mean",
          separeTsTv = T,
          workdir)


#Figure S01 ----
#Distribution of Frequencies of mutations per codon usage
# time line and boxplot for population HighLow
timeline(normalize = T,
         population = "HighLow",
         normBy = "mean",
         type = "bar",
         TsTv = "all",
         fix = T,
         save = T,
         workdir = workdir)
if (file.exists("FigS01.pdf")) 
  #Delete file if it exists
  file.remove("FigS01.pdf")
file.rename(from = file.path(figDir, "enrichbarHighLow.pdf"), 
            to = file.path(figDir, "FigS01.pdf"))

#Figure 02 ---- 
# all normalizations
figFunc05(figName = "Fig02Mean",
          workdir=workdir,
          normBy = "mean")
figFunc05(figName = "Fig02Count",
          workdir=workdir,
          normBy = "count")
figFunc05(figName = "Fig02CUB",
          workdir=workdir,
          normBy = "CUB")

#Fig S04/S05 ----
#Distribution of Frequencies of mutations per codon usage
# time line and boxplot A to F
figFunc04(figName = "FigS04CUB",workdir,type = "line",normBy = "CUB")
figFunc04(figName = "FigS05CUB",workdir,type = "box",normBy = "CUB")

figFunc04(figName = "FigS04Count",workdir,type = "line",normBy = "count")
figFunc04(figName = "FigS05Count",workdir,type = "box",normBy = "count")

figFunc04(figName = "FigS04Mean",workdir,type = "line",normBy = "mean")
figFunc04(figName = "FigS05Mean",workdir,type = "box",normBy = "mean")

#Figure S06 ----
#correlation beteween CAI and TAI
corrTaiCai(top=86,
           save = T,
           figName = "FigS06",
           workdir = workdir)

#Figure S07 ----
#frequency vs dw
figFunc06(figName = "FigS07",
          Dw = "Cai",
          normBy = "mean",
          workdir = workdir)
#Figure S08 ----
figFunc06(figName = "FigS08",
          Dw = "Tai",
          normBy = "mean",
          workdir = workdir)

#Figure 03 ----
#depletion Tai
nada<-plotEnrDeplPVal(type="High",
                pval=0.001,
                quant=5, 
                normalize = T,
                normBy = "mean",
                rank=200,
                fix = T,
                TsTv = "all",
                title = T,
                Dw = "Tai",
                top = 86,
                save=T,
                figName = "Fig03",
                workdir = workdir)

#Figure 04 ----
#depletion Cai
nada<-plotEnrDeplPVal(type="High",
                pval=0.001,
                quant=5, 
                normalize = T, 
                normBy = "mean",
                rank=200,
                fix = T,
                TsTv = "all",
                title = T,
                Dw = "Cai",
                top = 86,
                save=T,
                figName = "Fig04",
                workdir = workdir)




#Figure S09 ----
#all pvalues of hypergeometric test TAI
qt1<-plotEnrDeplAll(type="High",
                      pval=0.001,
                      quant=5, 
                      normalize = T, 
                      normBy = "mean",
                      rank=200,
                      fix = T,
                      TsTv = "all",
                      title = T,
                      Dw = "Tai",
                      top = 86,
                      save=T,
                      figName = "FigS09",
                      workdir = workdir)

#Figure S10 ----
#all pvalues of hypergeometric test CAI
qt1<-plotEnrDeplAll(type="High",
                    pval=0.001,
                    quant=5, 
                    normalize = T, 
                    normBy = "mean",
                    rank=200,
                    fix = T,
                    TsTv = "all",
                    title = T,
                    Dw = "Cai",
                    top = 86,
                    save=T,
                    figName = "FigS10",
                    workdir = workdir)

#Figure S11 ----
#test all range of Top ranked mutations for CAI
plotEnrRankRange(type="High",
                 pval=0.001,
                 quant=5, 
                 normalize = T, 
                 normBy = "count",
                 rankRange=c(20,1500,10),
                 fix = T,
                 TsTv = "Ts",
                 title = T,
                 Dw = "Cai",
                 top = 86,
                 refLine=200,
                 smooth = T,
                 save=T,
                 figName = "FigS11",
                 workdir = workdir)


#Figure S12 ----
#test all range of Top ranked mutations for TAI
plotEnrRankRange(type="High",
                 pval=0.001,
                 quant=5, 
                 normalize = T, 
                 normBy = "count",
                 rankRange=c(20,1500,10),
                 fix = T,
                 TsTv = "Ts",
                 title = T,
                 Dw = "Tai",
                 top = 86,
                 refLine=200,
                 smooth = T,
                 save=T,
                 figName = "FigS12",
                 workdir = workdir)


#Figure S13 ----
#test all range of Top ranked mutations for TAI
plotEnrTopRange(type="High",
                pval=0.001,
                quant=5, 
                normalize = T, 
                normBy = "count",
                rank=200,
                fix = T,
                TsTv = "Ts",
                title = T,
                Dw = "Cai",
                topRange = c(20,1500,10),
                refLine=86,
                smooth = T,
                save=T,
                figName = "FigS13",
                workdir = workdir)

#Info ----
#Other information functions
cat("###########################\n",
    "#     Umbalanced MutT     #\n",
    "###########################\n")
unbalancedMutT(top=86,
               workdir)

cat("###########################\n",
    "#   Enrichment of zeros   #\n",
    "###########################\n")
enrichZeros(quant=5,
            population = "High",
            Dw = "Cai",
            workdir )
  

cat("###########################\n",
    "# Transition/Transversion #\n",
    "###########################\n")
countTsTv(population = "High")

cat("###########################\n")
    listZeros(save = F,
          workdir = workdir,
          top = 86)
cat("###########################\n")
countPossibleTsTv(workdir)

cat("###########################\n")
analiseZeros(population = "HighLow",
             Dw = "Cai")
cat("###########################\n")
analiseZeros(population = "High",
             Dw = "Cai",top = 86)
cat("###########################\n")
analiseZeros(population = "MutT",
             Dw = "Cai",top = 86)
cat("###########################\n")
analiseZeros(population = "Low",
             Dw = "Cai",top = 86)
cat("###########################\n")
analiseZeros(population = "High",
             Dw = "Tai")
cat("###########################\n")
analiseZeros(population = "MutT",
             Dw = "Tai")
cat("###########################\n")
analiseZeros(population = "Low",
             Dw = "Tai")

cat("###########################\n")
corrMCxDw(Dw = "Tai", save = F, type = "point",quant = 5,
          workdir = workdir)

cat("###########################\n")
listGenes(top = 86,
          workdir = workdir)


# #depletion Tai just Ts High
# plotEnrDeplPVal(type="High",
#                 pval=0.001,
#                 quant=5, 
#                 normalize = T, 
#                 normBy = "count",
#                 rank=200,
#                 fix = T,
#                 TsTv = "Ts",
#                 title = T,
#                 Dw = "Tai",
#                 top = 86,
#                 save=F,
#                 figName = "FigS08",
#                 workdir = workdir)
# 
# #depletion Cai just Ts High
# plotEnrDeplPVal(type="High",
#                 pval=0.001,
#                 quant=5, 
#                 normalize = T, 
#                 normBy = "count",
#                 rank=200,
#                 fix = T,
#                 TsTv = "Ts",
#                 title = T,
#                 Dw = "Cai",
#                 top = 86,
#                 save=T,
#                 figName = "FigS09",
#                 workdir = workdir)
# #depletion Tai MutT
# plotEnrDeplPVal(type="MutT",
#                 pval=0.001,
#                 quant=5, 
#                 normalize = T, 
#                 normBy = "count",
#                 rank=200,
#                 fix = T,
#                 TsTv = "all",
#                 title = T,
#                 Dw = "Tai",
#                 top = 86,
#                 save=T,
#                 figName = "FigS10",
#                 workdir = workdir)
# 
# #depletion Cai MutT
# plotEnrDeplPVal(type="MutT",
#                 pval=0.001,
#                 quant=5, 
#                 normalize = T, 
#                 normBy = "count",
#                 rank=200,
#                 fix = T,
#                 TsTv = "all",
#                 title = T,
#                 Dw = "Cai",
#                 top = 86,
#                 save=T,
#                 figName = "FigS11",
#                 workdir = workdir)



cat("###########################\n")
totalCodons(top=86,
            workdir = workdir)

