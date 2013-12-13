# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 00:48:00 2013

@author: Talmo
"""
import numpy as np
import os

##### Constants #####
# The nucleotide dictionary used in all the functions.
bases = "ACGT"


##### Methods #####

def load_scaffolds(filename):
    """
    Loads a file into memory from a FASTA-formatted sequence file.
    
    Typically, genomic sequence data is stored in "scaffolds", i.e.:
        >scaffold1_1_MH0012
        TACTCTGGAAGGAGATATT...
    The scaffold is a reconstruction of one or more contigs (contiguous reads).
    
    This function also supports files that have multiple lines of sequence data
    per scaffold, i.e.:
        >gi|49175990|ref|NC_000913.2| Escherichia coli str. K-12 substr. MG1655 chromosome, complete genome
        AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
        TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA
        ...
    
    Args:
        filename: The path to the file containing the sequence data.
    
    Returns:
        sequence: The entire patient file loaded into one array. This is an
            integer numpy array where A = 0, C = 1, G = 2, T = 3.
        scaffolds: Since all the scaffolds are joined into one array, this is
            a list of the indices where each scaffold begins. This can be used
            to keep track of which segments of the sequence are contiguous.
            
    Example:
        >>> sequence, scaffolds = load_scaffolds("../NC_000913.fna") # The E. coli genome
        >>> sequence
        array([0, 2, 1, ..., 3, 3, 1])
        >>> sequence.dtype
        dtype('int32')
        >>> scaffolds
        array([0])
        >>> scaffolds.dtype
        dtype('int32')
        >>> int2nt(sequence[0:70]) # Display as a string
        'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
        
    Confirm that the first line is being read correctly:
        >>> ground_truth = open("../NC_000913.fna").readlines(1)[1].strip() # Read the first line
        >>> sequence, scaffolds = load_scaffolds("../NC_000913.fna")
        >>> hypothesis = int2nt(sequence[0:70])
        >>> hypothesis
        'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
        >>> hypothesis == ground_truth
        True
        
    """
    scaffolds = np.array([], dtype=np.int32)
    
    with open(filename, "r+") as f:
        # Pre-allocate sequence to size of file
        f.seek(0, os.SEEK_END)
        size = f.tell()
        f.seek(0, os.SEEK_SET)
        sequence = np.tile(np.array([-1], dtype=np.int32), size)
        
        # Process the lines in the file
        pos = 0
        next_line_starts_scaffold = False
        for line in f:
            if ">" in line:
                next_line_starts_scaffold = True
            else:
                # Keep sequence as the charcode integers
                scaffold_seq = np.array(line.strip(), "c").view(np.uint8)
                
                # Replace the sequence chunk in memory
                sequence[pos:pos + len(scaffold_seq)] = scaffold_seq
                
                if next_line_starts_scaffold:
                    # Keep track of the beginning of scaffolds
                    scaffolds = np.append(scaffolds, pos)
                    next_line_starts_scaffold = False
                
                # Update position in the sequence array
                pos += len(scaffold_seq)
                
    # Truncate the sequence array to fit the data
    sequence = sequence[0:np.where(sequence == -1)[0][0]]
    
    # Replace character codes with 0-3 int representation of the nucleotides
    for value, base in enumerate(bases):
        sequence[np.where(sequence == ord(base))] = value
        
    return sequence, scaffolds


def pssm2cols(pssm):
    """
    Returns the PSSM as a 2-dimensional array such that each 4 element array
    is a column vector of the PSSM, that is, each element in the first
    dimension represents a column of the PSSM, and each element in the second
    dimension represents the scores for each nucleotide.
    
    This can be used for any matrix structure, not just PSSMs, however.

    The reverse of this function is simply pssm2cols(pssm).flatten().
    
    To get an array of ROW vectors, transpose the output of this function.
    
    Args:
        pssm: An array of size 4 * w, where w is the width of the PSSM, or the
            number of columns it has. This array must be 1-dimensional.
    
    Returns:
        pssm_cols: A two dimensional array containing the values of pssm.
    
    Example:
        >>> cols = 4
        >>> pssm = np.array(range(0, 4 * cols))
        >>> pssm
        array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15])
        >>> pssm2cols(pssm)
        array([[ 0,  1,  2,  3],
               [ 4,  5,  6,  7],
               [ 8,  9, 10, 11],
               [12, 13, 14, 15]])
       >>> pssm2cols(pssm).transpose()
       array([[ 0,  4,  8, 12],
              [ 1,  5,  9, 13],
              [ 2,  6, 10, 14],
              [ 3,  7, 11, 15]])
       
       For a more readable representation, using print_pssm() which breaks the
       pssm into columns automatically if necessary:
       >>> print_pssm(pssm)
            |      0 |      1 |      2 |      3 |
        -----------------------------------------
        | A |  0.000 |  4.000 |  8.000 | 12.000 |
        | C |  1.000 |  5.000 |  9.000 | 13.000 |
        | G |  2.000 |  6.000 | 10.000 | 14.000 |
        | T |  3.000 |  7.000 | 11.000 | 15.000 |
        -----------------------------------------
    """
    pssm_cols = pssm.reshape((pssm.size / 4, 4))
    return pssm_cols


def permute_pssm(pssm):
    """
    Randomly shuffles the *columns* of the PSSM. This means that the value of
    each base is the same, but the positions are different.
    
    See the example for a clearer picture of how this works.
    
    Args:
        pssm: A 1D or 2D numpy array representing the PSSM.
    
    Returns:
        perm_pssm: The column-permuted pssm. This will maintain the dimensions
            of the PSSM, so if a 1D array was inputted, a 1D array is returned.
        
    Example:
        >>> cols = 4
        >>> pssm = np.array(range(0, 4 * cols))  # Generate a naive PSSM
        >>> print_pssm(pssm)
            |      0 |      1 |      2 |      3 |
        -----------------------------------------
        | A |  0.000 |  4.000 |  8.000 | 12.000 |
        | C |  1.000 |  5.000 |  9.000 | 13.000 |
        | G |  2.000 |  6.000 | 10.000 | 14.000 |
        | T |  3.000 |  7.000 | 11.000 | 15.000 |
        -----------------------------------------
        >>> pssm = permute_pssm(pssm)
        >>> print_pssm(pssm)
            |      0 |      1 |      2 |      3 |
        -----------------------------------------
        | A | 12.000 |  4.000 |  8.000 |  0.000 |
        | C | 13.000 |  5.000 |  9.000 |  1.000 |
        | G | 14.000 |  6.000 | 10.000 |  2.000 |
        | T | 15.000 |  7.000 | 11.000 |  3.000 |
        -----------------------------------------
    """
    _pssm = pssm
    if len(pssm.shape) == 1:
        # Break PSSM into columns if it is flat (i.e., 1D)
        _pssm = pssm2cols(pssm)
        
    # Shuffle columns of the PSSM
    perm_pssm = np.random.permutation(_pssm)
    
    if len(pssm.shape) == 1:
        # Re-flatten if it was originally a 1D array
        perm_pssm = perm_pssm.flatten()
    
    return perm_pssm


def pssm2str(pssm, max_width=130):
    """
    Renders the PSSM in a string in tabular form that is easily readable.
    
    Args:
        pssm: A 1- or 2-dimensional numpy array of PSSM values.
            Note: If the array is flat, or 1-dimensional, it will be converted
            internally to a 2d col array, so it does not need to be broken into
            columns before being passed to the function.
        max_width: The width of the console in characters. Columns of the PSSM
            will be wrapped past this number of characters.
    
    Returns:
        table: A string containing the table.
    
    See also: print_pssm() which is a wrapper for this function.
    """
    # Convert to 2D array if needed
    _pssm = pssm
    if len(pssm.shape) == 1:
        _pssm = pssm2cols(pssm)
    cols = int(_pssm.shape[0])
    rows = int(_pssm.shape[1])    
    
    # Calculate in how many rows we'll display the output
    cols_per_output = int(max_width - 4) / 9
    
    # Initialize output string
    table = ""
    
    # Construct the table
    for metacol in range(0, cols, cols_per_output):
        # Table header
        table += " " * 4 # Nucleotide label column
        for col in range(metacol, min(metacol + cols_per_output, cols)):
            table += "| %6d " % col
        table += "|\n"
        table += "-" * (9 * min(cols_per_output, min(metacol + cols_per_output, cols) - metacol) + 5)
        table += "\n"
        
        # Table body
        for row in range(rows):
            table += "| " + bases[row] + " "
            for col in range(metacol, min(metacol + cols_per_output, cols)):
                table += "| %6.3f " % (_pssm[col][row])
            table += "|\n"
            
        table += "-" * (9 * min(cols_per_output, min(metacol + cols_per_output, cols) - metacol) + 5)
        table += "\n"

    return table


def pssm2csv(pssm, delim=","):
    """
    Returns the PSSM as a comma-separated value string.
    
    Args:
        pssm: A 1- or 2-dimensional numpy array of PSSM values.
            Note: If the array is flat, or 1-dimensional, it will be converted
            internally to a 2d col array, so it does not need to be broken into
            columns before being passed to the function.
        delim: Separator character for cells. Alternatives to a comma might be
            a space (" ") or a tab character ("\t").
        
    Returns:
        csv: A CSV-formatted string representation of the PSSM.
    
    Example:
        >>> cols = 4
        >>> pssm = np.array(range(0, 4 * cols))  # Generate a naive PSSM
        >>> print pssm2csv(pssm)
        A,0,4,8,12
        C,1,5,9,13
        G,2,6,10,14
        T,3,7,11,15
        >>> print_pssm(pssm) # For comparison
            |      0 |      1 |      2 |      3 |
        -----------------------------------------
        | A |  0.000 |  4.000 |  8.000 | 12.000 |
        | C |  1.000 |  5.000 |  9.000 | 13.000 |
        | G |  2.000 |  6.000 | 10.000 | 14.000 |
        | T |  3.000 |  7.000 | 11.000 | 15.000 |
        -----------------------------------------
    """
    # Convert to 2D array if needed
    _pssm = pssm
    if len(pssm.shape) == 1:
        _pssm = pssm2cols(pssm)
    
    # Transpose the column vectors into row vectors
    _pssm = _pssm.transpose()

    # Initialize CSV string
    csv = ""
    
    # Build each row of the PSSM string
    for row, base in enumerate(bases):
        csv += base + delim
        csv += delim.join([str(cell) for cell in _pssm[row].tolist()])
        csv += "\n"
    
    return csv


def print_pssm(pssm, max_width=130):
    """
    Prints the PSSM in tabular form that is easily readable. This function is
    just a wrapper for pssm2str().
    
    Args:
        pssm: A 1- or 2-dimensional numpy array of PSSM values.
            Note: If the array is flat, or 1-dimensional, it will be converted
            internally to a 2d col array, so it does not need to be broken into
            columns before being passed to the function.
        max_width: The width of the console in characters. Columns of the PSSM
            will be wrapped past this number of characters.
        
    Example:
        >>> cols = 4
        >>> pssm = np.random.uniform(-7.5, 2.0, 4 * cols) # 1D array
        >>> print_pssm(pssm)
            |      0 |      1 |      2 |      3 |
        -----------------------------------------
        | A | -1.517 | -4.536 | -5.642 | -1.505 |
        | C | -7.157 |  0.742 | -3.313 | -3.493 |
        | G | -3.226 | -5.244 |  1.085 | -1.207 |
        | T | -2.886 | -4.886 |  1.919 | -4.437 |
        -----------------------------------------
    """
    print pssm2str(pssm, max_width=max_width)


def int2nt(int_seq):
    """
    Returns the string representation of a sequence from an integer array of
    nucleotides.
    
    The dictionary for this is at the top of this file and defaults to:
        bases = "ACGT"
    So A = 0, C = 1, G = 2, T = 3.

    The inverse of the function is nt2int().
    
    Args:
        int_seq: A 1D array of integer representations of nucleotides.
    
    Returns:
        seq: A string representation of the nucleotide sequence.
        
    Example:
        >>> sequence = np.random.randint(0, 3, 20) # Generate random 20 nucleotide sequence
        >>> sequence
        array([2, 1, 0, 2, 2, 2, 2, 0, 0, 1, 0, 2, 0, 1, 0, 1, 2, 0, 1, 2])
        >>> int2nt(sequence)
        'GCAGGGGAACAGACACGACG'
        
    Confirm that the function can convert the sequence reversibly:
        >>> np.array_equal(nt2int(int2nt(sequence)), sequence)
        True
    """
    seq = "".join([bases[b] for b in int_seq])
    return seq


def nt2int(nt_seq):
    """
    Returns an integer representation of a string sequence of nucleotides.
    
    This function is the inverse of int2nt().
    
    Args:
        nt_seq: A string representation of the nucleotide sequence.
    
    Returns:
        int_seq: A 1D array of integer representations of nucleotides.
        
    Example:
        >>> from random import choice
        >>> sequence = "".join([choice(bases) for _ in range(20)]) # Generate random 20 nucleotide sequence
        >>> sequence
        'CTTTCCCAAGTCTAGGGTGT'
        >>> nt2int(sequence)
        array([1, 3, 3, 3, 1, 1, 1, 0, 0, 2, 3, 1, 3, 0, 2, 2, 2, 3, 2, 3])
        >>> nt2int(sequence).dtype
        dtype('int32')
        
    Confirm that the function converts the sequence reversibly:
        >>> int2nt(nt2int(sequence)) == sequence
        True
    """
    # Convert the string into an array of characters and then as integers
    int_seq = np.array(nt_seq, "c").view(np.uint8).astype(np.int32)
    
    # Rescale the array of character codes to [0, 3]
    for value, base in enumerate(bases):
        int_seq[np.where(int_seq == ord(base))] = value
    
    return int_seq
    
    
def consensus(pssm):
    """
    Returns the consensus sequence of a PSSM.
    
    This is the sequence that has the best score for a given PSSM.
    
    Args:
        pssm: A 1D or 2D numpy array representing the PSSM.
    
    Returns:
        seq: An integer sequence representing the consensus of the PSSM.
    
    Example:
        >>> cols = 4
        >>> pssm = np.random.uniform(-7.5, 2.0, 4 * cols)
        >>> print_pssm(pssm)
            |      0 |      1 |      2 |      3 |
        -----------------------------------------
        | A | -7.021 | -5.289 |  1.621 | -3.932 |
        | C | -3.224 | -6.341 | -4.706 |  1.788 |
        | G | -3.844 | -3.741 | -1.344 | -4.757 |
        | T | -3.157 | -2.437 | -7.475 | -4.794 |
        -----------------------------------------
        >>> consensus(pssm)
        array([3, 3, 0, 1])
        >>> int2nt(consensus(pssm))
        'TTAC'
    """
    # Break PSSM into columns
    _pssm = pssm    
    if len(pssm.shape) == 1:
        _pssm = pssm2cols(pssm)
        
    seq = np.array([np.argmax(col) for col in _pssm], dtype=np.int32)
    return seq


def rev_pssm(pssm):
    """
    Returns the reverse-complement of a PSSM based on the Watson-Crick model of
    nucleotide pairing.
    
    Reversed PSSMs can be to score a sequence of nucleotides as if it were on
    the reverse strand.
    
    Args:
        pssm: A 1D or 2D numpy array representing the PSSM.
    
    Returns:
        rev_pssm: A 1D or 2D numpy array representing the reverse-complement of
            the PSSM.
    
    Example:
        >>> cols = 4
        >>> pssm = np.array(range(0, 4 * cols)) # Generate a naive PSSM
        >>> print_pssm(pssm)
            |      0 |      1 |      2 |      3 |
        -----------------------------------------
        | A |  0.000 |  4.000 |  8.000 | 12.000 |
        | C |  1.000 |  5.000 |  9.000 | 13.000 |
        | G |  2.000 |  6.000 | 10.000 | 14.000 |
        | T |  3.000 |  7.000 | 11.000 | 15.000 |
        -----------------------------------------
        >>> r_pssm = rev_pssm(pssm)
        >>> print_pssm(r_pssm)
            |      0 |      1 |      2 |      3 |
        -----------------------------------------
        | A |  3.000 |  2.000 |  1.000 |  0.000 |
        | C |  4.000 |  3.000 |  2.000 |  1.000 |
        | G |  5.000 |  4.000 |  3.000 |  2.000 |
        | T |  6.000 |  5.000 |  4.000 |  3.000 |
        -----------------------------------------
    """
    # Convert to 1D array if needed
    _pssm = pssm    
    if len(pssm.shape) == 2:
        _pssm = _pssm.flatten()
    
    # Get the reverse-complement
    _rev_pssm = np.array([_pssm[i / 4 + (3 - (i % 4))] for i in range(_pssm.size)][::-1])
    return _rev_pssm


def score_site(sequence, pssm):
    """
    Scores a sequence against a PSSM.
    
    This function is not optimized for large sequences, nor does it score the
    reverse strand of the sequence. If used on a large sequence, it will
    consume a very large amount of memory!
    
    Use it as a convenience function for scoring individual sites on a genome.
    
    Note that the width of the PSSM (the number of columns) must be equal to
    the length of the sequence!
    
    Args:
        sequence: A 1D array of integer representations of nucleotides.
        pssm: A 1D or 2D numpy array representing the PSSM.
    
    Returns:
        score: The number calculated by summing up the scores for each base at
            each position in the sequence.

    Example:
        >>> cols = 4
        >>> seq = np.random.randint(0, 3, cols) # Generate random sequence
        >>> pssm = np.random.uniform(-7.5, 2.0, 4 * cols) # Generate random PSSM
        >>> score_site(seq, pssm)
        -14.142570531696526
        >>> score_site(consensus(pssm), pssm) # The best possible score for this PSSM
        -2.5526269119209908
    """
    score = sum([pssm[i * 4 + j] for i, j in enumerate(sequence)])
    return score

