#
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 16:27:29 2019
Run a Monte Carlo symulation to check if the unobserved codons are relevant 
@author: Clovis Reis
"""

import numpy
import parse_file
import timecourse_utils
import random
import sys

#file name
arq=sys.argv[1]
##this simulation number
simulation=int(sys.argv[2])
##take how many runs the simulation will take
runs=int(sys.argv[3])
##take the syn mutations number for Monte Carlo
synMut=int(sys.argv[4])
print "Simulation ",simulation, "with ", synMut, " mutations, ", runs, "runs\n\tUsing files ", (arq+"Syn.csv and "+ arq+"All.csv\n")
#exit()
#synMut=5
random.seed(a=None)
#all nucleotides
nucBase='ATCG'
transitions=['A->G','G->A','T->C','C->T']
transversions = ['C->A','T->G','A->T','T->A','A->C','G->T','G->C','C->G']
gene_data = parse_file.parse_gene_list()
#create a dictionary for strand
geneData = dict((key, value) for (key, value) in zip(gene_data[0],gene_data[6]))
#all genome sequence
genome=parse_file.parse_reference_genome()

[position_gene_map,effective_gene_lengths, substitution_specific_synonymous_fraction]=parse_file.create_annotation_map()


#position=1353839
#ini=ini-1
#fim=+6
run=0
#create output file for synonymous and all mutations
#just synonymous
outFileSyn = open(arq+"Syn.csv","w")        
outFileSyn.write("%s\n" % ('simulation, run, position, symbol, type, ref, mut, nucRef, nucMut, allele, strand'))
#all mutations
outFileAll = open(arq+"All.csv","w")        
outFileAll.write("%s\n" % ('simulation, run, position, symbol, type, ref, mut, nucRef, nucMut, allele, strand'))

while run < runs:
    cont=0
    while cont < synMut :#for i in range(len(teste2)):
        #print(cont)
        #position=teste1[i]
        position=random.randint(0,len(genome))
        #position=2588552
        #print(position)
        #randomize a genome position
        geneName=parse_file.annotate_gene(position, position_gene_map)
        if geneName == 'repeat' or geneName == 'intergenic':
    #        print("###",position)
            continue
        strand=geneData[geneName]
        #find nucleotide at position
        nucRef=genome[position-1]
        #choose a nucleotide substitute
        nucList=nucBase.replace(nucRef,'')
        nucMut=nucList[random.randint(0,2)]
        #nucMut=teste2[i]
        #invert if strand reverse
    #    if strand == 'reverse':
    #        nucRef=parse_file.base_table[nucRef]
            #nucMut=parse_file.base_table[nucMut]
        
        
        allele=nucRef+"->"+nucMut
        #transition/transversion classification
        if allele in transitions:
            TsTv="Ts"
        elif allele in transversions:
            TsTv="Tv"
        else:
            TsTv="Un"
            
        [symbol,tipo,ref,mut,strd,seq,n1,n2]=parse_file.annotate_variant2(position, allele, gene_data, position_gene_map)
        #Anotate just transversions synonymous except Stop codons
        if tipo=='synonymous' and ref not in['TAA','TGA','TAG'] and TsTv == "Tv":# and strd =='reverse'):
            #write just Syn except stop codons
            outFileSyn.write("%d,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (simulation, run, position, symbol, tipo, ref, mut, nucRef, nucMut, allele, strd,TsTv))
            cont+=1
        #write anyway
        outFileAll.write("%d,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (simulation, run, position, symbol, tipo, ref, mut, nucRef, nucMut, allele, strd,TsTv))
    run+=1
outFileAll.close
outFileSyn.close
