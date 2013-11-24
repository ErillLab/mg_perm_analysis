# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 00:50:44 2013

@author: Talmo
"""

from gpu import gpu_pssm
import metagenomics as mg
import glob
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

def main():
    # Parameters
    metahit_path = "./Metahit"
    binding_sites_path = "./lexA.seq.fa"
    score_threshold = -100
    permutations = 10
    bins = range(-100, 24) # See figure S4   
    ##########
    
    # Load patient files into memory
    metahit_db = glob.glob(metahit_path + "/MH[0-9]*.seq.fa")
    
    # For debugging, truncate to just first patient
    metahit_db = [metahit_db[0]]
    
    for patient_file in metahit_db:
        # Load the sequence into memory
        metagenome, scaffolds = mg.load_scaffolds(patient_file)
        
        # Calculate the background nucleotide frequency for the metagenome
        mg_frequencies = np.bincount(metagenome).astype(np.float) / metagenome.size
        
        # Calculate the original PSSM from binding sites
        original_pssm = gpu_pssm.create_pssm(binding_sites_path, mg_frequencies)
        
        # Score the metagenome using the original PSSM
        original_scores = score_patient(metagenome, scaffolds, original_pssm, score_threshold, bins)
        
        # Keep the distributions of the scores
        patient_scores = [original_scores]
                
        for permutation in range(permutations):
            print "Permutation %d/%d..." % (permutation + 1, permutations)
            # Break the PSSM into columns
            pssm_cols = original_pssm.reshape((original_pssm.size / 4, 4))
            
            # Shuffle the columns and re-flatten the PSSM
            pssm = np.random.permutation(pssm_cols).flatten()
            
            # Re-score using permuted PSSM
            perm_scores = score_patient(metagenome, scaffolds, pssm, score_threshold, bins)            
            
            # Save the distribution of the scores
            patient_scores.append(perm_scores)
    
    # Plot results
    #print patient_scores
    for score in patient_scores[1:]:
        cdf = np.cumsum(score)
        plt.plot(bins[1:], cdf, "D-r", alpha=0.5, label="Permutation")
        
    cdf = np.cumsum(patient_scores[0])
    plt.plot(bins[1:], cdf, "D-b", lw=3, label="Original")
    plt.xlabel("Site score (bits)")
    plt.ylabel("Cumulative distribution")
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[-2:], labels[-2:], loc="best")
    plt.grid()
    
    plt.figure()
    for score in patient_scores[1:]:
        plt.plot(bins[1:], score, "D-r", alpha=0.5, label="Permutation")
    plt.plot(bins[1:], patient_scores[0], 'D-b', lw=3, label="Original")
    plt.xlabel("Site score (bits)")
    plt.ylabel("Probability istribution function")
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[-2:], labels[-2:], loc="best")
    plt.grid()
    
def score_patient(metagenome, scaffolds, pssm, score_threshold, bins):
    # Score the metagenome with the new PSSM
    scores = gpu_pssm.score_long_sequence(metagenome, pssm, keep_strands=False)

    # Invalidate cross-scaffold scores
    for pos in scaffolds:
        if pos - (pssm.size / 4) + 1 > 0:
            scores[pos - (pssm.size / 4) + 1:pos] = -np.inf

    # Eliminate scores below threshold        
    high_scores = np.where(scores >= score_threshold)
    scores = scores[high_scores]
    
    # Bin the scores
    score_hist, __ = np.histogram(scores, bins, density=True)
    
    return score_hist

if __name__ == "__main__":
    main()