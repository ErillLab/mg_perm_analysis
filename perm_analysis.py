# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 00:50:44 2013

@author: Talmo

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
from os import mkdir
import glob
import datetime
import time
import numpy as np
import matplotlib.pyplot as plt
import metagenomics as mg
# Import gpu_pssm from the parent directory
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from gpu_pssm import gpu_pssm


def main():
    ######################
    ##### Parameters #####
    ######################
    # Path to the MetaHit database.
    # The patient files can be downloaded from:
    #   http://www.bork.embl.de/~arumugam/Qin_et_al_2010/
    # Make sure these are extracted from their packages! (i.e.: "gunzip *.gz")
    # Note: After extraction, the 84 patient files occupy a total of 6.56 GB on
    # disk!
    # The original paper can be found at:
    #   http://www.nature.com/nature/journal/v464/n7285/full/nature08821.html
    metahit_path = "../MetaHit"
    
    # The collection of binding sites to generate the PSSM from.
    # LexA.seq.fa is a collection of 115 experimentally-determined binding
    # sites reported in literature. See Table S1 in Cornish et al. (2013).
    binding_sites_path = "./LexA.seq.fa"
    
    # The path to where the results of the permutation analysis will be saved.
    # A directory will be created with the date that the data was collected on.
    # If data already exists in this directory, a new folder will be created.
    # The result of the permutation analysis will be saved in CSV format for
    # each patient to its own file.
    # The output is a count of scores in each bin (a histogram).
    output_path = "../Permutations"
    save_output = True
    
    # Number of permutations to run. Note that the COLUMNS of the PSSM are the
    # ones permuted.
    num_permutations = 1
    
    # The probability density function and empirical cumulative distribution
    # function can be plotted by binning the results of scoring the entire
    # metagenome.
    # Note: See TODO below for why this is not currently working.
    plot_figures = False
    
    # Scores below this number of bits will not be reported.
    # Lower values will give more (false-positive) results and also slow down
    # the execution of the program since more memory needs to be allocated to
    # store the score values.
    score_threshold = -50.0
    
    # The ranges that the scores should be binned into.
    # Since no scores will be saved below the score_threshold, it serves as a
    # lower bound to the bins. A good upper bound is around ~30 bits since the
    # maximum theoretical score a sequence can have from a PSSM is 32.
    # Under equiprobable frequencies, the LexA consensus sequence has a score
    # of ~24 bits, so no sequence should score higher than that.
    bin_low = score_threshold
    bin_high = 24.5
    bin_step = 0.5
    
    
    ####################################################################
    ### Parameters below this line should *probably* not be changed. ###
    ####################################################################
    
    # The background frequency of the bases. An equiprobable frequency
    # distribution assumes that each base has an equal probability of occuring,
    # that is: P(A) = P(C) = P(G) = P(T) = 0.25 => GC-content = 0.5
    # If set to False, the background frequency will be calculated based on the
    # nucleotide composition of each patient.
    # For better comparison across patients, this should be set to True.
    equiprobable_nuc_freqs = True
    
    # For debugging, if you'd like to just score just the first patient, set
    # this parameter to True. This is useful wehn you want to just run a quick
    # set of permutations.
    only_score_first_patient = True
    
    # If you'd like to compare the output of this program by scoring the
    # E. coli genome, set this to True.
    # Make sure NC_000913.fna is in the parent directory.
    # The genome sequence can be download from:
    #  ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/
    # Use this for comparison/debugging.
    score_ecoli_instead = False
    path_to_ecoli_genome = "../NC_000913.fna"
    
    # Pre-load patient metagenomes into memory.
    # WARNING: This may  reduce runtime by avoiding the need to read the
    # sequences from the disk when every patient is scored, but it requires a
    # large amount of memory (RAM)!
    # To load all the patients in the MetaHit database (~6.9bi bp), a minimum
    # of 32 GB of memory is recommended.
    preload_patients = False
    
    # Set this to your preferred timing function from the time module.
    # On Windows, you should use time.clock, whereas on Unix you should use
    # time.time. See the comments in the gpu_pssm module for more info.
    timer = time.clock if sys.platform == 'win32' else time.time
    ####################################################################
    start_time = timer()
    
    # Find patient files on disk
    metahit_db = glob.glob(metahit_path + "/MH[0-9]*.seq.fa")
    
    # For debugging, truncate to just first patient
    if only_score_first_patient:
        metahit_db = [metahit_db[0]]
    
    # For debugging, score the E. coli genome instead of the MetaHit DB
    if score_ecoli_instead:
        metahit_db = [path_to_ecoli_genome]
        
    # Calculate bin array for score binning based on parameters
    bins = np.arange(bin_low, bin_high + bin_step, bin_step)
    
    # Load all patient sequences to memory before permuting
    if preload_patients:
        preload_start_time = timer()
        print "Preloading patient files (this may take a while)..."
        patients = []
        for patient_file in metahit_db:
            patients.append(mg.load_scaffolds(patient_file))
        preload_time = timer() - preload_start_time
        print "Done. All patient sequences in memory. Time elapsed: %.2fs" % preload_time
    
    # Assume equiprobable mononucleotide frequencies for PSSM calculation
    mg_frequencies = [0.25] * 4
    
    if not equiprobable_nuc_freqs:
        # Calculate the background nucleotide frequency for the metagenome
        #mg_frequencies = np.bincount(metagenome).astype(np.float) / metagenome.size
        # TODO: If we did want to compensate for background nucleotide
        # composition, what would we use? Frequency of entire metagenome? Of
        # each individual patient?
        pass
    
    # Calculate the original PSSM from binding sites
    original_pssm = gpu_pssm.create_pssm(binding_sites_path, genome_frequencies = mg_frequencies)
    
    # Calculate the randomly permuted PSSMs
    perm_pssms = [mg.permute_pssm(original_pssm) for _ in range(num_permutations)]

    # Initialize container for score distributions (the output of the analysis)
    patient_distributions = [[] for _ in range(len(metahit_db))]
    
    # Score each patient    
    bases_scored = np.int64(0) # for analytics
    for i in range(len(metahit_db)):
        print "Scoring patient %d/%d" % (i + 1, len(metahit_db)),
        
        # Load the patient sequence data
        if preload_patients:  # Get sequence from memory
            sequence, scaffolds = patients[i]
        else:  # Load the sequence from file
            patient_file = metahit_db[i]
            sequence, scaffolds = mg.load_scaffolds(patient_file)
        print "| Size:", sequence.size, "bp | Scaffolds:", scaffolds.size
        
        # Score the patient using each of the permuted PSSMs
        for permutation, pssm in enumerate([original_pssm] + perm_pssms):
            bases_scored += sequence.size
            score_start_time = timer()
            
            # Score the sequence
            scores = score_patient(sequence, scaffolds, pssm, score_threshold, bins)
            
            # Save the distribution of scores (the histogram)
            patient_distributions[i].append(scores)
            
            score_time = timer() - score_start_time
            print "=> Patient: %d/%d | Permutation: %d/%d => Elapsed: %.2fs (%.2g bp/s)" % (i + 1, len(metahit_db), permutation, num_permutations, score_time, sequence.size / score_time)
    
    total_time = timer() - start_time
    print "\nFinished permutation analysis. Total time elapsed: %.2fs (%.2g bp/s)" % (total_time, bases_scored / total_time)
        
    # Plot figures based on the patient scores
    if plot_figures:
        pass
    # TODO:
    # - Have this plot just using the first patient?
    # - Save figures to data structure that can be saved as PNG in the
    # save_output block below.
    # - Generate histogram, ECDF, PDF
    # - Use all patients?
#        plt.figure()
#        for score in patient_scores[1:]:
#            cdf = np.cumsum(score)
#            plt.plot(bins[1:], cdf, "D-r", alpha=0.5, label="Permutation")
#        
#        cdf = np.cumsum(patient_scores[0])
#        plt.plot(bins[1:], cdf, "D-b", lw=3, label="Original")
#        plt.xlabel("Site score (bits)")
#        plt.ylabel("Cumulative distribution")
#        handles, labels = plt.gca().get_legend_handles_labels()
#        plt.legend(handles[-2:], labels[-2:], loc="best")
#        plt.grid()
#        
#        plt.figure()
#        
#        for score in patient_scores[1:]:
#            plt.plot(bins[1:], score, "D-r", alpha=0.5, label="Permutation")
#        plt.plot(bins[1:], patient_scores[0], 'D-b', lw=3, label="Original")
#        plt.xlabel("Site score (bits)")
#        plt.ylabel("Probability distribution function")
#        handles, labels = plt.gca().get_legend_handles_labels()
#        plt.legend(handles[-2:], labels[-2:], loc="best")
#        plt.grid()
    
    
    # Save the patient scores to files
    if save_output:
        # Create folder to save output of the analysis
        path = create_output_folder(output_path)
        
        # Save parameters used in the analysis
        with open(path + "/" + "info.txt", "w") as f:
            f.write("Analysis date: %s\n" % datetime.datetime.now())
            f.write("\n")
            f.write("Patients scored: %d\n" % len(patient_distributions))
            f.write("Total bases: %d bases\n" % (bases_scored / num_permutations))
            f.write("Total bases scored: %d bases\n" % bases_scored)
            f.write("Time elapsed in analysis: %f sec\n" % total_time)
            f.write("Average speed: %f bp/sec\n" % (bases_scored / total_time))
            f.write("GPU: %s\n" % gpu_pssm.cuda.get_current_device().name)
            f.write("\n")
            f.write("Patients preloaded: %s\n" % preload_patients)
            if preload_patients:
                f.write("Preloading time elapsed: %f sec\n" % preload_time)
            f.write("\n")
            f.write("Binding sites: %s\n" % os.path.basename(binding_sites_path))
            f.write("Number of permutations: %d\n" % num_permutations)
            f.write("Score threshold: %f\n" % score_threshold)
            f.write("Bins: %f to %f in steps of %f\n" % (bin_low, bin_high, bin_step))
            f.write("Background frequencies: %s\n" % ", ".join([base + " = " + str(mg_frequencies[i]) for i, base in enumerate(mg.bases)]))
            f.write("Score threshold: %f\n" % score_threshold)
            f.write("\n")
            f.write("Original PSSM:\n%s\n" % mg.pssm2str(original_pssm, max_width=1000))
        
        # Save the original and permuted PSSMs used to score
        with open(path + "/" + "PSSMs.csv", "w") as f:
            for i, pssm in enumerate([original_pssm] + perm_pssms):
                header = "Permutation:,%d%s\n" % (i, "," * (pssm.size / 4 - 1))
                f.write(header + mg.pssm2csv(pssm))
        
        # Save the scores for each patient to separate CSV        
        for patient_file, scores in zip(metahit_db, patient_distributions):
            # Get the "MH_00##" name of the patient
            patient_name = os.path.basename(patient_file).split(".")[0]
            
            # Header for columns (number corresponds to permutation number)
            header_row = "bins," + ",".join([str(p) for p in range(len(scores))])
            
            # Column of bin labels (left-most column)
            bins_col = np.array(bins[:-1]).reshape(-1, 1)
            
            # Set up the data table array
            scores_table = np.hstack([score.reshape(-1, 1) for score in scores])

            # Save patient scores to CSV
            np.savetxt(path + "/data/" + patient_name + ".csv", np.hstack((bins_col, scores_table)), delimiter=",", header=header_row, comments="")
        print "Output of permutation analysis saved to:"
        print os.path.abspath(path)
    

def score_patient(metagenome, scaffolds, pssm, score_threshold, bins):
    # Score the metagenome with the new PSSM
    scores = gpu_pssm.score_long_sequence(metagenome, pssm, keep_strands=False)

    # Invalidate cross-scaffold scores by setting them below threshold
    for pos in scaffolds:
        if pos - (pssm.size / 4) + 1 > 0:
            scores[pos - (pssm.size / 4) + 1:pos] = -np.inf

    # Eliminate scores below threshold
    high_scores = np.where(scores >= score_threshold)
    scores = scores[high_scores]
    
    # Bin the scores
    score_hist, __ = np.histogram(scores, bins, density=False)
    
    return score_hist


def create_output_folder(output_path):
    # The base path is a folder with today's date in the output folder
    base_path = output_path + "/" + str(datetime.date.today())
    
    # Try the base path without a number
    path = base_path
    i = 0
    
    # Keep adding to the number at the end of the filename until folder doesn't
    # already exist
    while os.path.exists(path):
        i += 1
        path = base_path + " (%d)" % i
    
    # Create the folders for the output files
    mkdir(path)
    mkdir(path + "/data")
    mkdir(path + "/figures")

    return path


if __name__ == "__main__":
    main()