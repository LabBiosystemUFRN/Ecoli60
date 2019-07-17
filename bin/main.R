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
cat(sep="","###########################\n",
    "Calculating Delta w TAI \n",
    "###########################\n")
calcDWTAI(type = "Max",workdir = workdir)
#create a table with all High and Low mutations
cat(sep="","###########################\n",
    "Creating table HighLow \n",
    "###########################\n")
joinHighLow(workdir = workdir)
#Create a MutT list
cat(sep="","###########################\n",
    "Creating MutT mutation list\n",
    "###########################\n")

mutationsMutT(workdir = workdir)

#Figure 01 ----
#Number of mutation per population & type
cat(sep="","###########################\n",
    "Figure 01 \n",
    "###########################\n")

figFunc01(figName = "Fig01",
          workdir)
comparePopulations(figName = "Fig01a",
                   save = T,
                   workdir = workdir)


#Figure 02 ---- 
# all normalizations
cat(sep="","###########################\n",
    "Figure 02 \n",
    "###########################\n")

figFunc05(figName = "Fig02Mean",
          workdir=workdir,
          normBy = "mean")
figFunc05(figName = "Fig02Count",
          workdir=workdir,
          normBy = "count")
figFunc05(figName = "Fig02CUB",
          workdir=workdir,
          normBy = "CUB")

#Figure 03 ----
#depletion Tai
cat(sep="","###########################\n",
    "Figure 03  \n",
    "###########################\n")

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
cat(sep="","###########################\n",
    "Figure 04 \n",
    "###########################\n")

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

#Figure S01 ----
#Distribution of Frequencies of mutations per codon usage
# time line and boxplot for population HighLow
cat(sep="","###########################\n",
    "Figure S01 \n",
    "###########################\n")
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


#Figure S02 ----
#Frequencie of mutations per codon usage A to F
cat(sep="","###########################\n",
    "Figure S02 \n",
    "###########################\n")

figFunc02(figName = "FigS02",
          separeTsTv = T,
          workdir)

#Figure S03 ----
#Frequencie of mutations per codon usage after normalization A to F
cat(sep="","###########################\n",
    "Figure S03 \n",
    "###########################\n")
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


#Fig S04/S05 ----
#Distribution of Frequencies of mutations per codon usage
# time line and boxplot A to F
cat(sep="","###########################\n",
    "Figure S04 and S05  \n",
    "###########################\n")

figFunc04(figName = "FigS04CUB",workdir,type = "line",normBy = "CUB")
figFunc04(figName = "FigS05CUB",workdir,type = "box",normBy = "CUB")

figFunc04(figName = "FigS04Count",workdir,type = "line",normBy = "count")
figFunc04(figName = "FigS05Count",workdir,type = "box",normBy = "count")

figFunc04(figName = "FigS04Mean",workdir,type = "line",normBy = "mean")
figFunc04(figName = "FigS05Mean",workdir,type = "box",normBy = "mean")

#Figure S06 ----
#correlation beteween CAI and TAI
cat(sep="","###########################\n",
    "Figure S06 \n",
    "###########################\n")

corrTaiCai(top=86,
           save = T,
           figName = "FigS06",
           workdir = workdir)

#Figure S07 ----
#frequency vs dw
cat(sep="","###########################\n",
    "Figure S07 \n",
    "###########################\n")

figFunc06(figName = "FigS07",
          Dw = "Cai",
          normBy = "mean",
          workdir = workdir)
#Figure S08 ----
cat(sep="","###########################\n",
    "Figure S08 \n",
    "###########################\n")

figFunc06(figName = "FigS08",
          Dw = "Tai",
          normBy = "mean",
          workdir = workdir)





#Figure S09 ----
#all pvalues of hypergeometric test TAI
cat(sep="","###########################\n",
    "Figure S09 \n",
    "###########################\n")

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
cat(sep="","###########################\n",
    "Figure S10 \n",
    "###########################\n")

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
cat(sep="","###########################\n",
    "Figure S11 \n",
    "###########################\n")

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
cat(sep="","###########################\n",
    "Figure S12 \n",
    "###########################\n")

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
cat(sep="","###########################\n",
    "Figure S13 \n",
    "###########################\n")

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
cat(sep="","\n###########################\n",
    "#     Umbalanced MutT     #\n",
    "###########################\n")
unbalancedMutT(top=86,
               workdir)

cat(sep="","\n###########################\n",
    "#   Enrichment of zeros   #\n",
    "###########################\n")
enrichZeros(quant=5,
            population = "High",
            Dw = "Cai",
            workdir )


cat(sep="","\n###########################\n",
    "# Transition/Transversion #\n",
    "###########################\n")
countTsTv(population = "High")

cat(sep="","\n##########################################\n")
listZeros(save = F,
          workdir = workdir,
          top = 86)
cat(sep="","\n##########################################\n")
countPossibleTsTv(workdir)

cat(sep="","\n##########################################\n")
analiseZeros(population = "HighLow",
             Dw = "Cai")
cat(sep="","\n##########################################\n")
analiseZeros(population = "High",
             Dw = "Cai",top = 86)
cat(sep="","\n##########################################\n")
analiseZeros(population = "MutT",
             Dw = "Cai",top = 86)
cat(sep="","\n##########################################\n")
analiseZeros(population = "Low",
             Dw = "Cai",top = 86)
cat(sep="","\n##########################################\n")
analiseZeros(population = "High",
             Dw = "Tai")
cat(sep="","\n##########################################\n")
analiseZeros(population = "MutT",
             Dw = "Tai")
cat(sep="","\n##########################################\n")
analiseZeros(population = "Low",
             Dw = "Tai")

cat(sep="","\n##########################################\n")
corrMCxDw(Dw = "Tai", save = F, type = "point",quant = 5,
          workdir = workdir)

cat(sep="","\n##########################################\n")
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



cat(sep="","\n##########################################\n")
totalCodons(top=86,
            workdir = workdir)

