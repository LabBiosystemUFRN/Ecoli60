#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 16:27:29 2019
Acrescenta informações do codon aos arquivos *SynFreq.csv das 
famílias pouco mutagenicas
@author: Clovis Reis
"""

import numpy
import sys
import parse_file
import timecourse_utils
import itertools
import math


dataDirIn='../data_files/AuxFiles/'
dataDirOut='../data_files/AuxFiles/'
addData='../additional_data/'

def checkSNP(permCodon):
    index = 3
    retorno = []
    for index in range(len(permCodon)):
        pair = permCodon[index]
        ref = pair[0]            
        mut = pair[1]
        cont = 0 
        for i in range(3):
            if ref[i] != mut[i]:
                cont+=1
        if cont == 1:
            retorno.append(pair)
    return(retorno)
    

def readGenes():
    #read gene_data to create a dictionary
    for index in range(len(gene_data[0])):
        #Key gene; start stop strand sequence
        if gene_data[6][index]=='forward':
            geneDict[gene_data[0][index]]=[gene_data[1][index],\
                     gene_data[2][index],\
                     gene_data[6][index],\
                     gene_data[5][index]]
        else:
            geneDict[gene_data[0][index]]=[gene_data[1][index],\
                     gene_data[2][index],\
                     gene_data[6][index],\
                     parse_file.calculate_reverse_complement_sequence(gene_data[5][index])]

def saveCodonsGenes():
    gene='sfmA'
    outFileName = dataDirOut+"cCodonsGenes.csv"
    outFile = open(outFileName,"w")        
    outFile.write("%s\n" % ('gene,codon,count'))
    #exclude stop codons
    allCodons = [k for k,v in parse_file.codon_table.items()]# if v != '!']
    for gene in geneDict.keys():
        #create a hash for codon counts. all codons start in 0
        countCodon = {}
        for k in allCodons:
            countCodon[k] = 0
        geneSequence = geneDict[gene][3]
        split_points = range(0,len(geneSequence),3)
        split_points.append(len(geneSequence))
        #split codons
        codons = [geneSequence[i: j] for i, j in zip(split_points, split_points[1:])]
        # compute each codon occurence
        index = 'CTT'
        for index in codons:
            countCodon[index] = countCodon[index] + 1
        codon = 'CTT'
        for codon in countCodon.keys():
            lineOut=gene+"," + codon + "," + str(countCodon[codon])
            #print(lineOut)
            outFile.write("%s\n" % (lineOut))

    outFile.close()
        

        
    




#gene_data Já contém os dados dos genes
gene_data = parse_file.parse_gene_list()
repeat_data = parse_file.parse_repeat_list()
mask_data = parse_file.parse_mask_list()

#dictionary of gene sequences
geneDict =	{}

readGenes()

saveCodonsGenes()

#codon couts in high expressed genes
codonCount = {k:[0,0,0,0] for k in parse_file.codon_table.keys()}
#dict.fromkeys(parse_file.codon_table)
codonW = {k:[0.0,0.0,0.0,0.0,parse_file.codon_table[k]] for k in parse_file.codon_table.keys() \
          if parse_file.codon_table[k] != '!'} #Exclude stop codons


#dict.fromkeys(parse_file.codon_table,0.0)
#codonCount = dict.copy(parse_file.codon_table)
#for k,v in codonCount.items():
#    codonCount[k] = [v,0]

line = "thrL"        
#read highly Expressed Genes
files=["highlyEGLeGac"]
for indexF in range(len(files)):
    count=0
    filename=dataDirIn+files[indexF]
    inFile = open(filename,"r")
    for line in inFile:
        if line.find("#") == -1:
            geneName = line.rstrip()
            #exclude genes unknow genes symbols
            if geneName in geneDict:
                geneSequence = geneDict[geneName][3]
                split_points = range(0,len(geneSequence),3)
                split_points.append(len(geneSequence))
                #split codons
                codons = [geneSequence[i: j] for i, j in zip(split_points, split_points[1:])]
                # compute each codon occurence
                index = 'CTT'
                for index in codons:
        #            codonCount[codons[index]] = codonCount[codons[index]]+1
                    codonCount[index][indexF] = codonCount[index][indexF]+1
                    count=count+1
            else:
                print("Symbol "+geneName+" not found. Ignoring...")
    #just checking the totalization. Must be 0
    check=0
    for k in codonCount.keys():
        check=check+codonCount[k][indexF]
    if count == check:
        print(files[indexF]+" ... Ok")
    else:
        print(files[indexF]+" ... Ok")
    
   
#print(count-sum([C for A,C in codonCount.values()]))
        
#calculate w for each syn codon
#list of synonimous per aminoacid
#create delta W
amino = "Y"
for amino in (set(parse_file.codon_table.values())):
    if amino == '!': #Exclude stop codons
        continue
    synCodons = [k for k,v in parse_file.codon_table.items() if v == amino]
    indexF=0
    #calculate w
    for indexF in range(len(files)):
        counts=[codonCount[k][indexF] for k in synCodons]
        indexS =0
        #calculate W
        for indexS in range(len(synCodons)):
            if max(counts) == 0:
                codonW[synCodons[indexS]][indexF] = 'Erro'
            else:
                codonW[synCodons[indexS]][indexF] = float(counts[indexS])/float(max(counts))



# write Ws to a file
outFileName = dataDirOut+"cCodonsW.csv"
outFile = open(outFileName,"w")        
outFile.write("%s\n" % ('amino,ref,w'))
amino = "Y"
for codon in codonW.keys():
    lineOut=codonW[codon][4]+","+codon+","+str(codonW[codon][0])#+","+str(codonW[codon][1])+","+str(codonW[codon][2])
    #print(lineOut)
    outFile.write("%s\n" % (lineOut))
outFile.close()    



#calculate delta w for each synonymous codon combination
outFileName = dataDirOut+"cDeltaW.csv"
outFile = open(outFileName,"w")
#header
outFile.write("%s\n" % ('ref,mut,dw'))

    
deltaW = []
amino = "C"
for amino in (set(parse_file.codon_table.values())):
    if amino == '!': #Exclude stop codons
        continue
    synCodons = [k for k,v in parse_file.codon_table.items() if v == amino]
    indexF=1
    #permutatios of codons
    permCodon=list(itertools.permutations(synCodons, 2))
    permCodon=checkSNP(permCodon)
#    if amino == 'C':
#        print("Erro: ", amino)
    indexP = 0
    for indexP in range(len(permCodon)):
        ref = permCodon[indexP][0]
        mut = permCodon[indexP][1]
        lineOut = ref + ',' + mut 
        for indexF in range(len(files)):
            wRef = codonW[ref][indexF]
            wMut = codonW[mut][indexF]
            if wRef == 0 or wMut == 0 or not isinstance(wRef,float) or not isinstance(wMut,float):
#                lineOut = lineOut + ',' + repr(wRef)+ ',' +repr(wMut)+ ','+'Erro'
                lineOut = lineOut + ',' + 'Erro'
                continue
            else:
                dW = math.log10(wMut/wRef)
#            #just for the first time
#            if indexF == 0:
#                deltaW.append([ref,mut,dW,0.0,0.0])
#            lineOut = lineOut + ',' + repr(wRef)+ ',' +repr(wMut)+ ',' +repr(dW) 
            lineOut = lineOut + ',' + repr(dW) 
        #print(lineOut)
        outFile.write("%s\n" % (lineOut))


outFile.close()

        
#        parse_file.codon_table['AGA']
# 'GCA': 'A',
# 'GCC': 'A',
# 'GCG': 'A',
# 'GCT': 'A',
#        
#        [codonCount['GCA'],\
#        codonCount['GCC'],\
#        codonCount['GCG'],\
#        codonCount['GCT']]
#        
#        [184.0/307.0,41.0/307.0,123.0/307.0,307.0/307.0]
#        [382.0/773.0,547.0/773.0,773.0/773.0,312.0/773.0]
#        [184.0/307.0,41.0/307.0,123.0/307.0,307.0/307.0]
#        math.log10(0.5993485342019544/0.13355048859934854)
#        math.log10(0.49417852522639066/0.7076326002587322)
#        
#        [codonW['GCA'],\
#        codonW['GCC'],\
#        codonW['GCG'],\
#        codonW['GCT']]
