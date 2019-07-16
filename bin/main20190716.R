#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#           main processing 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#This script generate the figures in 
# the same order than it apears in teh text
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#inicio ----
rm(list = ls())

baseDir="/home/clovis/Dropbox/Ecoli60/"
#baseDir="/home/clovis/Doutorado/Projetos/Ecoli60/"

binDir=paste0(baseDir,"bin/")
workdir=paste0(baseDir,"data_files/")
figDir=paste0(baseDir,"figures/")
setwd(workdir)
setwd(binDir)
source(paste0(binDir,"allFunctions.R"))

calcDWTAI(type = "Max",workdir = workdir)

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
         save = T)
if (file.exists("FigS01.pdf")) 
  #Delete file if it exists
  file.remove("FigS01.pdf")
file.rename(from = file.path(figDir, "enrichbarHigh.pdf"), 
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
           figName = "figS06",
           workdir = workdir)

#Figure S07 ----
#frequency vs dw
figFunc06(figName = "figS07",
          Dw = "Cai",
          normBy = "mean",
          workdir = workdir)
#Figure S08 ----
figFunc06(figName = "figS08",
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
plotEnrDeplPVal(type="High",
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







timeline(file="AuxFiles/bSumHigh.csv",
         normalize = T, 
         type = "box",
         save = T, #save option is not working here
         normType = "perK",
         mutT = 0,
         workdir = workdir)





unbalancedMutT(top=86,
               workdir)


enrichZeros()


countTsTv(population = "High")

analiseZeros(population = "HighLow",
             Dw = "Cai")

#createMonteCarloDF()
listZeros(save = F,
          workdir = workdir,
          top = 86)

countPossibleTsTv(workdir)

analiseZeros(population = "High",
             Dw = "Cai",top = 86)
analiseZeros(population = "MutT",
             Dw = "Cai",top = 86)
analiseZeros(population = "Low",
             Dw = "Cai",top = 86)
analiseZeros(population = "High",
             Dw = "Tai")
analiseZeros(population = "MutT",
             Dw = "Tai")
analiseZeros(population = "Low",
             Dw = "Tai")


corrMCxDw(Dw = "Tai", save = F, type = "point",quant = 5,
          workdir = workdir)



calcDWTAI(type = "Min",workdir = workdir)

source(paste0(binDir,"allFunctions.R"))

plotEnrRankRange(type="High",
                pval=0.001,
                quant=5, 
                normalize = T, 
                normBy = "count",  
                rankRange=c(50,1500,100),
                fix = T,
                TsTv = "Ts",
                title = T,
                Dw = "Tai",
                top = 86,
                save=F,
                figName = "FigSXTaiTs",
                workdir = workdir)

source(paste0(binDir,"allFunctions.R"))

rankRange=c(50,1500,50)
top = 86
type="High"
pval=0.001
quant=5
normalize = T
normBy = "count"
rank=200
fix = T
TsTv = "Ts"
title = T
Dw = "Cai"
topRange = c(20,1500,100)
save=F
figName = "FigS2XCaiTs2"
workdir = workdir

source(paste0(binDir,"allFunctions.R"))


listGenes(top = 86,
          workdir = workdir)


#depletion Tai just Ts High
plotEnrDeplPVal(type="High",
                pval=0.001,
                quant=5, 
                normalize = T, 
                normBy = "count",
                rank=200,
                fix = T,
                TsTv = "Ts",
                title = T,
                Dw = "Tai",
                top = 86,
                save=F,
                figName = "FigS08",
                workdir = workdir)

#depletion Cai just Ts High
plotEnrDeplPVal(type="High",
                pval=0.001,
                quant=5, 
                normalize = T, 
                normBy = "count",
                rank=200,
                fix = T,
                TsTv = "Ts",
                title = T,
                Dw = "Cai",
                top = 86,
                save=T,
                figName = "FigS09",
                workdir = workdir)
#depletion Tai MutT
plotEnrDeplPVal(type="MutT",
                pval=0.001,
                quant=5, 
                normalize = T, 
                normBy = "count",
                rank=200,
                fix = T,
                TsTv = "all",
                title = T,
                Dw = "Tai",
                top = 86,
                save=T,
                figName = "FigS10",
                workdir = workdir)

#depletion Cai MutT
plotEnrDeplPVal(type="MutT",
                pval=0.001,
                quant=5, 
                normalize = T, 
                normBy = "count",
                rank=200,
                fix = T,
                TsTv = "all",
                title = T,
                Dw = "Cai",
                top = 86,
                save=T,
                figName = "FigS11",
                workdir = workdir)



totalCodons(workdir = workdir)

plotCorrCaiOld(workdir=workdir,
         save = F)

correlation()

corrTaiCai(save = T,figName = "corrTaiCai")


nada<-plotEnrDeplPVal(type="High",
                pval=0.001,
                quant=5, 
                normalize = T, 
                normBy = "count",
                rank=100,
                fix = T,
                TsTv = "Ts",
                title = T,
                Dw = "Tai",
                save=F,
                workdir = workdir)



tRNA<-read.table("/home/clovis/Dropbox/Ecoli60/data_files/TAIFiles/REL606tRNA.txt",
           sep="\t",
           header = T,
           stringsAsFactors = F)
TsTv<-read.csv("AuxFiles/TsTv.csv",
                     header = T,
                     stringsAsFactors = F)
codons<-TsTv[!duplicated(TsTv$ref),c('amino',"ref")]
aminos<-c("A","Ala","Alanine","R","Arg","Arginine","N","Asn","Asparagine","D","Asp","Aspartic acid","C","Cys","Cysteine","Q","Gln","Glutamine","E","Glu","Glutamic acid","G","Gly","Glycine","H","His","Histidine","I","Ile","Isoleucine","L","Leu","Leucine","K","Lys","Lysine","M","Met","Methionine","F","Phe","Phenylalanine","P","Pro","Proline","S","Ser","Serine","T","Thr","Threonine","W","Trp","Tryptophan","Y","Tyr","Tyrosine","V","Val","Valine","!","stp","Stop")
aminos<-as.data.frame(t(matrix(aminos,nrow = 3,ncol = 21)),
                      stringsAsFactors = F)
colnames(aminos)<-c("amino","abrev","name")
tRNA$codon<-apply(X=as.data.frame(tRNA$AntiCodonsList), 
                  MARGIN = 1,
                  convertToComplement)
tRNA<-merge(tRNA,codons,by.x="codon",by.y="ref", all.x = T)
tRNA$amino[tRNA$codon=="ATG"]<-"M"
tRNA$amino[tRNA$codon=="TGG"]<-"W"
tRNA$amino[tRNA$codon%in%c("TAA","TAG","TGA")]<-"!"
tRNA[is.na(tRNA$amino),]
tRNA<-merge(tRNA,aminos, by="amino")
#tRNA$AntiCodonsList<-NULL

write.csv(tRNA,file= "/home/clovis/Dropbox/Ecoli60/presentation/figuras/aux/tRNA.csv")
