#
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 16:27:29 2019
Add informations of reference codon, mutated codon, mutation strand, 
    population, time when the mutation arrise, time of fixation, 
    and transit to *SynFreq.csv and MutAllFreq.csv for all populations
@author: Clovis Reis
"""

import numpy
import sys
import os
import parse_file
import timecourse_utils

import numpy
from math import log10
import pylab
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as pe
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from numpy.random import normal, choice, shuffle

#import figure_utils
import stats_utils

dataDirIn='../data_files/preProcess/'
dataDirOut='../data_files/AuxFiles/'
addData='../additional_data/'


gene_data = parse_file.parse_gene_list()
repeat_data = parse_file.parse_repeat_list()
mask_data = parse_file.parse_mask_list()

position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction = parse_file.create_annotation_map(gene_data, repeat_data, mask_data)

#dictionary of gene sequences
geneDict =	{}

def calculateFixation(population):
    #fitness_trajectories = {}
    #mutation_trajectories = {}
    #fixed_mutation_trajectories = {}
    
    transit_times = {}
    
    sys.stderr.write("\tProcessing fixation...\n")
    #population='p6'


    # calculate mutation trajectories
    # load mutations
    mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
    
    population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    
    dummy_times,fmajors,fminors,haplotype_trajectories = parse_file.parse_haplotype_timecourse(population)
    state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)
        
    times = mutations[0][10]
    #Ms = numpy.zeros_like(times)*1.0
    #fixed_Ms = numpy.zeros_like(times)*1.0
    
    mutation_idx=1
    for mutation_idx in xrange(0,len(mutations)):
 
        location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx] 
    
        #Ls = haplotype_trajectories[mutation_idx]
        state_Ls = state_trajectories[mutation_idx]
        
        good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
        
        freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
        
        masked_times = times[good_idxs]
        masked_freqs = freqs[good_idxs]
        masked_state_Ls = state_Ls[good_idxs]
        
        t0,tf,transit_time =timecourse_utils.calculate_appearance_fixation_time_from_hmm(masked_times, masked_freqs, masked_state_Ls)
        if(tf == 1000000):
            tf = 'NA'
        transit_times[(location, allele)] = (str(t0),str(tf),str(transit_time))

    sys.stderr.write("\tanalyzed %d mutations!\n" % len(mutations))
        
    return(transit_times)
        
          



def annotateCodon(inFileName, outFileName, popName, append):
    #open out file
    if append:
        outFile = open(outFileName,"a")
    else:
        outFile = open(outFileName,"w")

    #header=0
    totMut=0
    comas=""

#    print(lines)
#    filename=dataDir+("m5SynFreq.csv" % popName)
    #processing fixation
    sys.stderr.write("Processing %s...\n" % parse_file.get_pretty_name(popName))
    transit_times = calculateFixation(popName)
    
    sys.stderr.write("\tProcessing mutations...\n")

    #print("Processing population: "+ popName)
    inFile = open(inFileName,"r")
    line=inFile.readline() # skip headers
    headItens=len(line.split(","))
    if not append:
        outFile.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (line.rstrip(), "ref", "mut", "strand","pop","t0","fixed","transit"))
        append = True
    for line in inFile:
        items = line.split(",")
        geneName = items[1]
        location = int(items[0])
        allele = items[2]
        annotation=items[3]
        if annotation in ["indel","noncoding","sv"]:
            continue
        # annotate mutation
        gene_name, var_type, codon, new_codon,strand, gene_sequence, gene_start_position, gene_end_position = parse_file.annotate_variant2(location, allele, gene_data, position_gene_map)
        geneDict[gene_name]=[strand,gene_start_position,gene_end_position,gene_sequence]

        if gene_name != geneName:
            print("Deu merda!")
        #print(gene_name+" "+var_type+" "+ codon+" "+ allele+" "+ new_codon+" "+ strand+" "+popName)
        if headItens-len(items) >= 0:
            comas=""
        
            #adjust to header size
            for cComa in range(headItens-len(items)):
                comas=comas+","
        else:
            #limit extra fields
            tmp=items[0]
            for contItens in range(1,126):
                tmp=tmp+","+items[contItens]
            #teste= tmp.split(",")
            line=tmp
        #print(comas)
        t0,tf,transit=transit_times[(location,allele)]
        outFile.write("%s,%s%s,%s,%s,%s,%s,%s,%s\n" % (line.rstrip(),comas, codon, new_codon, strand,popName,t0,tf,transit))
        totMut += 1
            
    outFile.close()
    inFile.close()
    return(totMut)


#############################################################


i=0

mutator_lines = [["High",['m2','m3','m4','p3']],["Low",['m5','m6','p1','p2','p4','p5']],["MutT",['m1','p6']]]

#all lines Synonymous
for i in range(3):
    lines = mutator_lines[i][1]
    #to remove file if the first time in lineage
    totMut = 0
    append = False
    #delFile = 1
    #popName="p6"
    for popName in lines:
        inFileName = dataDirIn+("%sSynFreq.csv" % popName)        
        outFileName = dataDirOut+"a"+mutator_lines[i][0]+"MutSynFreqErroTGC.csv"
        #remove preexistent file
#        if os.path.exists(outFileName) and delFile:
#            os.remove(outFileName)
#            delFile = 0
            
        totMut += annotateCodon(inFileName, outFileName, popName, append)
        append = True
    sys.stderr.write('Mutations processed for '+ mutator_lines[i][0]+": "+ str(totMut)+"\n")

#all lines all mutations
for i in range(3):
    lines = mutator_lines[i][1]
    #to remove file if the first time in lineage
    totMut = 0
    append = False
    #delFile = 1
    #popName="p6"
    for popName in lines:
        inFileName = dataDirIn+("%sFreq.csv" % popName)        
        outFileName = dataDirOut+"a"+mutator_lines[i][0]+"MutAllFreqErroTGC.csv"
        #remove preexistent file
#        if os.path.exists(outFileName) and delFile:
#            os.remove(outFileName)
#            delFile = 0
            
        totMut += annotateCodon(inFileName, outFileName, popName, append)
        append = True
    sys.stderr.write('Mutations processed for All '+mutator_lines[i][0]+": "+ str(totMut)+"\n")
