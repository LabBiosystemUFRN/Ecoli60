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

qt1<-plotEnrDeplAll(type="High",
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
                      figName = "FigTeste",
                      workdir = workdir)


resPlot<-plotEnrDeplPVal(type="High",
                pval=0.001,
                quant=5, 
                normalize = T, 
                normBy = "mean",
                rank=200,
                fix = T,
                TsTv = "Ts",
                title = T,
                Dw = "Tai",
                top = 86,
                save=T,
                figName = "FigTeste2",
                workdir = workdir)



#Number of mutation per population & type
figFunc01(figName = "Fig01",
          workdir)
comparePopulations(figName = "Fig01a",
                    save = T,
                   workdir = workdir)

#Frequencie of mutations per codon usage A to F
figFunc02(figName = "FigS01",
          separeTsTv = T,
          workdir)


#Frequencie of mutations per codon usage after normalization A to F
figFunc03(figName = "FigS02Count",
          normBy="count",
          separeTsTv = T,
          workdir)

#Distribution of Frequencies of mutations per codon usage
# time line and boxplot A to F
figFunc04(figName = "FigS03",workdir,type = "line")
figFunc04(figName = "FigS04",workdir,type = "box")


timeline(file="AuxFiles/bSumHigh.csv",
         normalize = T, 
         type = "box",
         save = T, #save option is not working here
         normType = "perK",
         mutT = 0,
         workdir = workdir)

#Distribution of Frequencies of mutations per codon usage
# time line and boxplot for population HighLow
figFunc05(figName = "Fig02",workdir)


corrTaiCai(top=86,
           save = T,
           figName = "figS05",
           workdir = workdir)

figFunc06(figName = "figS06",
          Dw = "Cai",
          workdir = workdir)
figFunc06(figName = "figS07",
          Dw = "Tai",
          workdir = workdir)


unbalancedMutT(top=86,
               workdir)


enrichZeros()


countTsTv(population = "HighLow")

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
                 smooth = F,
                 save=T,
                 figName = "FigSXTaiTs2",
                 workdir = workdir)

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
                topRange = c(20,1500,100),
                refLine=86,
                smooth = T,
                save=T,
                figName = "FigS2XCaiTsS",
                workdir = workdir)

listGenes(top = 86,
          workdir = workdir)

#depletion Tai
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
                save=T,
                figName = "Fig03",
                workdir = workdir)

#depletion Cai
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
                figName = "Fig04",
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


plotEnrDeplPVal(type="High",
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
