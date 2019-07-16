#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 16:27:29 2019
Create the files to be used by stAIcalc program
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
dataDirOut='../data_files/stAIcalcFiles/'
addData='../additional_data/'


#gene_data Já contém os dados dos genes
gene_data = parse_file.parse_gene_list()
repeat_data = parse_file.parse_repeat_list()
mask_data = parse_file.parse_mask_list()

#dictionary of gene sequences
geneDict =	{}

readGenes()


#Open fasta file to save result to disk
outFileName = dataDirOut+"REL606CdsFasta.txt"
outFile = open(outFileName,"w")


geneName = "accA"        
for geneName in geneDict.keys():
    if geneName.find("#") == -1:
        geneSequence = geneDict[geneName.rstrip()][3]
        if len(geneSequence)== 0:
            print("No bases found:", geneName,". Skipping...")
            continue
        #Write gene name 
        outFile.write(">%s\n" % (geneName))
            
        split_points = range(0,len(geneSequence),70)
        split_points.append(len(geneSequence))
        #split codons
        lines = [geneSequence[i: j] for i, j in zip(split_points, split_points[1:])]
        # compute each codon occurence
        index = 0
        for index in range(0,len(lines)):
            outFile.write("%s\n" % (lines[index]))
        outFile.write("\n")
outFile.close()

