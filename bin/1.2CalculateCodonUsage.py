#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 16:27:29 2019
Create a file with thecodon usage of entire genome
Field factorPerK is the normalization base in the number of codosn of entire genome per thousand,and
    factorAbs is the normalization using the lowers codon usage as reference
@author: Clovis Reis
"""

import numpy
import sys
import parse_file
import timecourse_utils
import itertools
import math

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


dataDirIn='../data_files/preProcess/'
dataDirOut='../data_files/AuxFiles/'
addData='../additional_data/'


#gene_data Já contém os dados dos genes
gene_data = parse_file.parse_gene_list()
repeat_data = parse_file.parse_repeat_list()
mask_data = parse_file.parse_mask_list()

#dictionary of gene sequences
geneDict =	{}

readGenes()

#codon couts inside entire genome
codonCount = {k:0 for k in parse_file.codon_table.keys()}
#dict.fromkeys(parse_file.codon_table)

count = 0
geneName = "thrL"        
for geneName in geneDict.keys():
    if geneName.find("#") == -1:
        geneSequence = geneDict[geneName.rstrip()][3]
        if int(len(geneSequence)/3)*3 != len(geneSequence):
            print("Numero de bases incorreto no Gene ", geneName)
            continue
            
        split_points = range(0,len(geneSequence),3)
        split_points.append(len(geneSequence))
        #split codons
        codons = [geneSequence[i: j] for i, j in zip(split_points, split_points[1:])]
        # compute each codon occurence
        index = 'CTT'
        for index in codons:
#            codonCount[codons[index]] = codonCount[codons[index]]+1
            codonCount[index] = codonCount[index]+1
            count=count+1
    #just checking the totalization. Must be 0
check=0
for k in codonCount.keys():
    check=check+codonCount[k]
if count == check:
    print("Contagem ... Ok")
else:
    print("Contagem ... Fail")
    
   
#save result to disk
outFileName = dataDirOut+"dCodonUsage.csv"
outFile = open(outFileName,"w")
#header
outFile.write("%s\n" % ('ref,count,factorPerK,factorAbs'))

#remove stop codons
codonCount.pop("TAA",None)
codonCount.pop("TGA",None)
codonCount.pop("TAG",None)
#find the codon with minimum score
vMin=min(codonCount.values())
vTot=sum(codonCount.values())

codon = "CTT"
for codon  in codonCount.keys():
    lineOut = codon + ',' + repr(codonCount[codon]) + ',' +repr(float(codonCount[codon])/float(vTot)*1000)+ ',' + repr(float(codonCount[codon])/float(vMin))
    print(lineOut)
    outFile.write("%s\n" % (lineOut))


outFile.close()

outFileName = dataDirOut+"dAminoCodon.csv"
outFile = open(outFileName,"w")
#header
outFile.write("%s\n" % ('amino,ref'))      
for codon in parse_file.codon_table.keys():
    outFile.write("%s,%s\n" % (parse_file.codon_table[codon],codon))
outFile.close()

