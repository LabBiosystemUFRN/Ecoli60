#usa python 2.7
###############################
#
# Rest of script begins here
#
################################
import os
#path="/home/clovis/Doutorado/Projetos/Ecoli60/bin/"
#os.chdir(path)

import pylab
import numpy
import sys
from math import log10
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from numpy.random import binomial
import bz2
import parse_file
import matplotlib
import matplotlib.pyplot as plt
import timecourse_utils

debug=False
remove_data_for_copyediting = False

focal_populations = parse_file.complete_nonmutator_lines + parse_file.mutator_lines
all_populations = parse_file.complete_nonmutator_lines + parse_file.mutator_lines
all_colors = parse_file.nonmutator_line_colors+parse_file.mutator_line_colors

remaining_populations = ['m5','p1','p4','p5','m2','m3','m4','p3','p6']


########################################
#
# Calculating
#
########################################



if debug==True:
    starting_idx = -5
else:
    starting_idx = 0
   
theory_times = numpy.arange(0,122)*500

gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = parse_file.parse_gene_list()

gene_name_position_map = {gene_names[i]: (start_positions[i], end_positions[i]) for i in xrange(0,len(gene_names))}

allIdx = numpy.arange(0,122)
  
for i in xrange(0,len(focal_populations)):

    sys.stderr.write('Processing focal population: %s...\t' % focal_populations[i])
    filename = parse_file.data_directory+("%sFreq.csv" % focal_populations[i])
    outFile = open(filename, "w")
    outFile.write("Position,Gene,Allele,Annotation,")
    outFile.write("%s\n" % ", ".join([str(t) for t in theory_times]))

    filename = parse_file.data_directory+("%sSynFreq.csv" % focal_populations[i])
    synFile = open(filename, "w")
    synFile.write("Position,Gene,Allele,Annotation,")
    synFile.write("%s\n" % ", ".join([str(t) for t in theory_times]))

    

    # load mutations
    mutations, depth_tuple = parse_file.parse_annotated_timecourse(focal_populations[i])
    
    population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    
    dummy_times,fmajors,fminors,haplotype_trajectories = parse_file.parse_haplotype_timecourse(focal_populations[i])
    
    for mutation_idx in range(0,len(mutations))[starting_idx:]:
 
        location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx] 
    
        Ls = haplotype_trajectories[mutation_idx]
    
        good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
        
        freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
        data={}
        for cont in range(0,len(times[good_idxs])):
            data[times[cont]] = freqs[cont]
        masked_freqs2 = []    
        for cont in range(0,len(theory_times)):
            if theory_times[cont] in data.keys():
                masked_freqs2.append(  data[theory_times[cont]]  )
            else:
                masked_freqs2.append(-1)
            
#        badIdx = list(set(allIdx).difference(good_idxs))
        masked_times = times
        masked_freqs = freqs
        
#        for bad in badIdx:
#            masked_freqs[bad] = -1
        
        
        interpolation_function = timecourse_utils.create_interpolation_function(masked_times, masked_freqs, tmax=100000)
        haplotype_interpolation_function = timecourse_utils.create_interpolation_function(times, Ls, tmax=100000,kind='nearest')
    
        clone_freqs = clone_alts*1.0/(clone_depths+(clone_depths==0))
        clone_depth_fold_changes = timecourse_utils.estimate_depth_fold_changes(clone_avg_depths, clone_depths)
        
#        if not remove_data_for_copyediting:
#            focal_axes[i].plot(masked_times, masked_freqs, '-o', alpha=0.5, markersize=1, markeredgecolor='none', zorder=4, linewidth=0.5) #,rasterized=True)
    
        outFile.write("%s,%s,%s,%s," % (location, gene_name, allele, var_type))
        outFile.write("%s\n" % ", ".join([str(l) for l in masked_freqs]))
            
        if var_type=='synonymous':
            synFile.write("%s,%s,%s,%s," % (location, gene_name, allele, var_type))
            synFile.write("%s\n" % ", ".join([str(l) for l in masked_freqs2]))

    sys.stderr.write("Done!\n")
    outFile.close
    synFile.close
    
    
    
sys.stderr.write("Done!\n")
#sys.exit(1)
