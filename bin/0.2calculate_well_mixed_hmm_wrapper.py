import timecourse_utils
import parse_file
from math import log,exp
import numpy
from scipy.special import gammaln as loggamma
import sys
from scipy.optimize import fmin   
from scipy.interpolate import interp1d
from scipy.stats import nbinom     

import uniform_prior_well_mixed_hmm as well_mixed_model

populations = parse_file.all_lines

min_coverage = 5

for population in populations:
 
    mutations = []
    good_mutations = []
 
    sys.stderr.write("Processing %s\n" % population)
 
    mutation_data, depth_tuple = parse_file.parse_annotated_timecourse(population)
    
    times = mutation_data[0][10]
    
    for mutation_idx in xrange(0,len(mutation_data)):
 
        location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutation_data[mutation_idx] 
        
        good_idxs, masked_alts, masked_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue, min_coverage)
 
        max_freq = (masked_alts*1.0/(masked_depths+(masked_depths<1))).max()
        num_good_timepoints = (masked_depths>0.5).sum()
        
        if num_good_timepoints > 5:
            mutations.append((masked_alts*1.0, masked_depths*1.0))
            if var_type!='sv' and var_type!='indel' and masked_depths[-1]>1:
            #if masked_depths[-1]>1:
                good_mutations.append((masked_alts*1.0, masked_depths*1.0))
        else:
            mutations.append((numpy.zeros_like(times)*1.0, numpy.zeros_like(times)))
    
    
    A = []
    D = []
    for i in xrange(0,len(mutations)):
        Ai,Di = mutations[i]
        A.append(Ai)
        D.append(Di)
        
    A = numpy.array(A)
    D = numpy.array(D)
    
    Pstate, Ls, p0, p1 = well_mixed_model.infer_hmm(A,D,num_iterations=10)
    
    # decode ML states
    #Ls = numpy.zeros_like(A) 
    #for i in xrange(0,A.shape[0]):
    #    Ls[i] = (Pstate[i,:,:]).argmax(axis=1)
    
    # make thing that says if filtered:
    
    filtered_mutations = numpy.ones_like(A)
    for i in xrange(0,A.shape[0]):
        if D[i,:].sum() < 0.5:
            filtered_mutations[i,:] = numpy.zeros_like(times)
    
    ns = numpy.zeros((Pstate.shape[2],Pstate.shape[1]))
    
    for l in xrange(0,ns.shape[0]):
        ns[l,:] = (Pstate[:,:,l]*filtered_mutations).sum(axis=0)
      
    #Ls = numpy.zeros_like(A) 
    #for i in xrange(0,A.shape[0]):
    #    if filtered_mutations[i,0] < 0.5:
    #        Ls[i] = numpy.zeros_like(times)
    #    else:
    #        Ls[i] = (Pstate[i,:,:]).argmax(axis=1)
    
    hard_ns = numpy.zeros_like(ns) 
    for l in xrange(0,ns.shape[0]):
        hard_ns[l,:] = (Ls==l).sum(axis=0)
         
    # Write stuff
    haplotype_filename = parse_file.data_directory+("%s_well_mixed_state_timecourse.txt" % population)
    haplotype_file = open(haplotype_filename, "w")
    haplotype_file.write("%s\n" % ", ".join([str(t) for t in times]))
    haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[0,:]]))
    haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[1,:]]))
    haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[2,:]]))
    haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[3,:]]))
    for i in xrange(0,Ls.shape[0]):
        haplotype_file.write("%s\n" % ", ".join([str(l) for l in Ls[i,:]]))
        
    haplotype_file.close()