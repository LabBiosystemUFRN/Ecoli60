#!/usr/bin/env python2.7
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



#gene_data Já contém os dados dos genes
gene_data = parse_file.parse_gene_list()
repeat_data = parse_file.parse_repeat_list()
mask_data = parse_file.parse_mask_list()

#dictionary of gene sequences
geneDict =	{}

readGenes()

print "Number of genes: "+str(len(geneDict))