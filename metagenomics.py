# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 00:48:00 2013

@author: Talmo
"""
import numpy as np
import os

# Constants
bases = "ACGT"


def load_scaffolds(filename):
    """
    Loads a file into memory from a FASTA-formatted sequence file.
    
    Typically, genomic sequence data is stored in "scaffolds", i.e.:
        >scaffold1_1_MH0012
        TACTCTGGAAGGAGATATT...
    The scaffold is a reconstruction of one or more contigs (contiguous reads).
    
    Args:
        filename: The path to the file containing the sequence data.
    
    Returns:
        sequence: The entire patient file loaded into one array. This is a
            uint8 NumPy array where A = 0, C = 1, G = 2, T = 3.
        scaffolds: Since all the scaffolds are joined into one array, this is
            a list of the indices where each scaffold begins. This can be used
            to keep track of which segments of the sequence are contiguous.
        
    """
    scaffolds = []
    
    with open(filename, "r+") as f:
        # Pre-allocate sequence to size of file
        f.seek(0, os.SEEK_END)
        size = f.tell()
        f.seek(0, os.SEEK_SET)
        sequence = np.tile(-1, size)
        
        # Process the lines in the file
        pos = 0
        for line in f:
            if ">" not in line:
                # Keep sequence as the charcode integers
                scaffold_seq = np.array(line.strip(), "c").view(np.uint8)
                
                # Replace the sequence chunk in memory
                sequence[pos:pos + len(scaffold_seq)] = scaffold_seq
                             
                # Keep track of the beginning of scaffolds
                scaffolds.append(pos)
                
                # Update position in the sequence array
                pos += len(scaffold_seq)
                
    # Truncate the sequence array to fit the data
    sequence = sequence[0:np.where(sequence == -1)[0][0]]
    
    # Replace character codes with 0-3 int representation of the nucleotides
    for value, base in enumerate(bases):
        sequence[np.where(sequence == ord(base))] = value
        
    return sequence, scaffolds


def permute_pssm(pssm):
    # Break PSSM into columns
    pssm_cols = pssm.reshape((pssm.size / 4, 4))
    
    # Shuffle columns and return flattened array
    return np.random.permutation(pssm_cols).flatten()
