# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 00:50:44 2013

@author: Talmo
@editor: David  

This script does a permutation analysis on the scores of the LexA binding motif
in the MetaHit database in order to support the findings in Cornish et al.
(2013).

Dependencies:
    - NumPy (http://www.numpy.org/)
    - matplotlib (http://matplotlib.org/)
    - gpu_pssm (https://github.com/ErillLab/gpu_pssm)
    - MetaHit database (Qin et al., 2010; see comments below)
"""
import sys
import os.path
import glob
import numpy as np
import matplotlib.pyplot as plt
import time
import metagenomics as mg
# Import gpu_pssm from the parent directory
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from gpu_pssm import gpu_pssm

def main():
    
    ##### Parameters #####
    #Want to run with promoter regions
    promoter = True
    
    # Path to the MetaHit database.
    # The patient files can be downloaded from:
    #   http://www.bork.embl.de/~arumugam/Qin_et_al_2010/
    # Make sure these are extracted from their packages! (i.e.: "gunzip *.gz")
    # Note: After extraction, the 84 patient files occupy a total of 6.56 GB on
    # disk!
    # The original paper can be found at:
    #   http://www.nature.com/nature/journal/v464/n7285/full/nature08821.html
    if promoter:
        metahit_path = "./MetaHit/Pruned"
    else:
        metahit_path = "./MetaHit/Data"
    
    # The collection of binding sites to generate the PSSM from.
    # LexA.seq.fa is a collection of 115 experimentally-determined binding
    # sites reported in literature. See Table S1 in Cornish et al. (2013).
    # binding_sites_path = "./LexA.seq.fa"
    #binding_sites_path= "./LexA_Gamma_collection.fas"
    binding_sites_path = "./LexA_Grampos_collection.fas"
    
    # Number of permutations to run. Note that the COLUMNS of the PSSM
    permutations = 50
    
    # Scores below this number of bits will not be reported.
    # Lower values will give more (false-positive) results and also slow down
    # the execution of the program since more memory needs to be allocated to
    # store the score values.
    score_threshold = -50.0
    
    ### Parameters below this line should *probably* not be changed. ###
    
    # The background frequency of the bases. An equiprobable frequency
    # distribution assumes that each base has an equal probability of occuring,
    # that is: P(A) = P(C) = P(G) = P(T) = 0.25 => GC-content = 0.5
    # If set to False, the background frequency will be calculated based on the
    # nucleotide composition of each patient.
    # For better comparison across patients, this should be set to True.
    equiprobable_nuc_freqs = True
    
    # The ranges that the scores should be binned into.
    # Since no scores will be saved below the score_threshold, it serves as a
    # lower bound to the bins. A good upper bound is around ~30 bits since the
    # maximum theoretical score a sequence can have from a PSSM is 32.
    # Under equiprobable frequencies, the LexA consensus sequence has a score
    # of ~24 bits, so no sequence should score higher than that.
    bins = range(int(score_threshold), 32, 1)
    
    # If you'd like to compare the output of this program by scoring the
    # E. coli genome, set this to True.
    # Make sure NC_000913.fna is in the parent directory.
    # The genome sequence can be download from:
    #  ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/
    # Use this for comparison/debugging.
    score_ecoli_instead = False			
    
    #Calculate the total number of sites, scaffolds scanned
    total_num_sites = 0
    total_size = 0
    total_scaffold = 0
    ##########

    # Find patient files on disk
    if promoter:
        metahit_db = glob.glob(metahit_path + "/Pruned_MH[0-9]*.seq.fa")
    else:
        metahit_db = glob.glob(metahit_path + "/MH[0-9]*.seq.fa")
       #  metahit_db = ["Eco_300_1_50_P.txt"]
         
    if score_ecoli_instead:
        metahit_db = ["../NC_000913.fna"] # E. coli genome (for debugging)
    
    # gather total time
    start = time.time()
				
    # For debugging, truncate to just first patient
    # An alert just incase I am clumsy and forget these lines are uncommented.
    #print "USING ONLY ONE PATIENT!!"
    #metahit_db = ["Eco_300_1_50_P.txt"]
    #metahit_db = [metahit_db[0]]
    # Assume equiprobable mononucleotide frequencies
    mg_frequencies = [0.25] * 4
        
    if not equiprobable_nuc_freqs:
        # Calculate the background nucleotide frequency for the metagenome
        mg_frequencies = np.bincount(metagenome).astype(np.float) / metagenome.size
        
    # Calculate the original PSSM from binding sites
    original_pssm = gpu_pssm.create_pssm(binding_sites_path, genome_frequencies = mg_frequencies)
        
    # Print the unpermuted PSSM        
    print "Unpermuted PSSM:"
    mg.print_pssm(original_pssm)

    #preallocate the array to speed up process
    permute_pssm = np.empty(shape=(permutations+1,len(original_pssm)), dtype=object)
    patient_scores = np.zeros(shape=(permutations+1,len(bins)-1))
	
    #calculate predetermined pssm
    for permutation in range(permutations):	
       
        # Permute the PSSM
        permute_pssm[permutation] = mg.permute_pssm(original_pssm)

    #Cycle through every patient file
    for patient_file in metahit_db:

        #Status Update
        print "File: ", patient_file

        # Load the sequence into memory
        metagenome, scaffolds = mg.load_scaffolds(patient_file)
        total_size += metagenome.size
        total_scaffold += scaffolds.size
        print "Genome size:", metagenome.size, "| Scaffolds:", scaffolds.size
           
        # Score the metagenome using the original PSSM
        print "Scoring without permuting..."
        original_scores,partial_num_sites = score_patient(metagenome, scaffolds, original_pssm, score_threshold, bins)

        # Keep the distributions of the scores
        patient_scores[0]+= original_scores
        total_num_sites += partial_num_sites

        #For each permutated pssm
        for permutation in range(permutations):

            #which permutation it is on
            print "Permutation %d/%d..." % (permutation + 1, permutations)            
            
            # Re-score using permuted PSSM
            perm_scores,partial_num_sites = score_patient(metagenome, scaffolds, permute_pssm[permutation], score_threshold, bins)            
            
            # Save the distribution of the scores
            patient_scores[permutation+1] += perm_scores

    # Plot results
    plt.figure()
    for score in patient_scores[1:]:
        cdf = np.cumsum(score)
        plt.plot(bins[1:], cdf, "D-r", alpha=0.5, label="Permutation")
    
    cdf = np.cumsum(patient_scores[0])
    plt.plot(bins[1:], cdf, "D-b", lw=3, label="Original")
    plt.xlabel("Site score (bits)")
    plt.ylabel("# of Sites Found")
    plt.title("Cumulative Density Function")
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[-2:], labels[-2:], loc="best")
    plt.grid()
    
    plt.figure()
    
    for score in patient_scores[1:]:
        plt.plot(bins[1:], score, "D-r", alpha=0.5, label="Permutation")
    plt.plot(bins[1:], patient_scores[0], 'D-b', lw=3, label="Original")
    plt.xlabel("Site score (bits)")
    plt.ylabel("# of Sites Found ")
    plt.title("Probability Density Function")
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[-2:], labels[-2:], loc="best")
    plt.grid()

    #Calculate p-values
    total_fake_patient_scores = patient_scores[1]
    for score in patient_scores[2:]:
	   total_fake_patient_scores += score
    total_fake_patient_scores = np.vectorize(lambda x: x if x > 0 else .000001)(np.float64(total_fake_patient_scores))

    print
    real_prob = np.float64(patient_scores[0])/total_num_sites
    fake_prob = total_fake_patient_scores/(total_num_sites * permutations)

    #print out the p-values
    print "Probability (True Matrix | Score):"
    matrix_prob = (real_prob/(fake_prob+real_prob)) * 100
    print "\n".join("%d:%.2f" % (s, p) for s, p in zip(bins,matrix_prob))

    #Plot the Probability Values
    plt.figure()
    plt.bar(bins[1:],matrix_prob)
    plt.xlabel("Site score (bits)")
    plt.ylabel("Confidence Level (%)")
    plt.grid()
    
    #final diagnostics
    end = time.time()
    print "Total time: %.2f seconds" % (end-start)	
    print "Total Metagenome Size: %d bp" % (total_size)
    print "Total Scaffolds Scanned: %d" % (total_scaffold)	

def score_patient(metagenome, scaffolds, pssm, score_threshold, bins):
    # Score the metagenome with the new PSSM
    scores = gpu_pssm.score_long_sequence(metagenome, pssm, keep_strands=False)
    #scores, __ = gpu_pssm.score_sequence_with_cpu(metagenome, pssm, benchmark=False)

    # Invalidate cross-scaffold scores
    for pos in scaffolds:
        if pos - (pssm.size / 4) + 1 > 0:
            scores[pos - (pssm.size / 4) + 1:pos] = -np.inf

    partial_num_sites = len(np.where(scores > -np.inf)[0])
    # Eliminate scores below threshold
    high_scores = np.where(scores >= score_threshold)
    #print np.where(scores >= score_threshold)l
    #print scores.size
    #print len(scaffolds)
    print "Scores >= 20 bits:", len(np.where(scores >= 20)[0])
    #print len(high_scores[0])
    scores = scores[high_scores]
    
    # Bin the scores
    score_hist, __ = np.histogram(scores, bins, density=False)
    
    return (score_hist, partial_num_sites)

if __name__ == "__main__":
    main()