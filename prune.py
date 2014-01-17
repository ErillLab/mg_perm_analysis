"""
@author: David 

This script parses the BGI_GeneSet20090523_annotation file and grabs the intergenic regions
from patient files. 

Dependencies:
	- Numpy

Assumes:
	 1. Files to be parsed are in lexicographic order.
	 2. Ignores original strand regions that start at the very beginning of a sequence
	 3. Ignores compliment strand regions that end of the very end of a sequence
"""

import re, os, sys, time
import numpy as np

### Global Settings ###
#Name of the file that holds the gene annotations
gene_prediction = "../MetaGenome/BGI_GeneSet20090523_annotation"
#Regex pattern that matches the patient files
filepattern = "MH[\d]+"
#Regex pattern that matches the scaffold/contig lines in the annotation and patient files
gene_region_pattern = "((?:scaffold|c)[\d]+[^\d][\d]+):([\d]+):([\d]+):([+-])"
#Regex pattern that matches the scaffold/contig label 
scaffold_pattern  = "((?:scaffold|c)[\d]+[^\d][\d]+)"
#Regex pattern that matches the start and stop regions
region_pos_pattern = ":([\d]+):([\d]+):([+-])"
#Grab the accession number from each line
accession_number_pattern = "GL[\d]+"
accession_number_pattern_object = re.compile(accession_number_pattern, re.IGNORECASE)
#A set that holds the file listings to keep track of unique files
file_list = set()

### Methods ###

def preprocess_file(filename):
	"""

	Loads the gene annotation data into memory.

	Typical line in the annotation file:
	GL0017522_MH0001_[Lack_5'-end]_[mRNA]_locus=scaffold14_9:2:700:+	COG0850	K03610

	The patterns specified in the global variable section will match the MH#### part of the line:
	GL0017522_(MH0001)_[Lack_5'-end]_[mRNA]_locus=scaffold14_9:2:700:+	COG0850	K03610
	Then it will add the entire line into the array to be parsed by the prune_files method below.
	
	Args:
		filename: The file path towards the gene annotation file.

	Returns:
		Returns an array that holds all the relevant lines from the annotation file. 
		It only holds lines that match with the pattern mentioned above. 

	"""
	
	#Compile the pattern to decrease the time it takes to search
	file_pattern_object = re.compile(filepattern)
	gene_region_pattern_object = re.compile(gene_region_pattern, re.IGNORECASE)

	with open(filename) as g:
		#Trick borrowed from Talmo
		#Preallocates array based on the number of lines in the file
		g.seek(0,os.SEEK_END)
		pos = g.tell()
		gene_locations = np.empty(pos, dtype="S128")
		g.seek(0,os.SEEK_SET)
		index = 0

		#process the lines in the file
		for line in g:
			#Attempts to match each line with the regex pattern MH[\d]+
			match = file_pattern_object.search(line)

			#if match is found insert entire line into list
			if match:
				file_list.add(match.group(0))
				line_match = gene_region_pattern_object.search(line)
				accession_number_match = accession_number_pattern_object.search(line)
				gene_locations[index] = accession_number_match.group(0) + ":" + line_match.group(0)
				index += 1
	
	return gene_locations

def prune_files(gene_locations):
	"""
	Finds the intergenic regions and writes a new file containing the regions
	Each file will hold as output:
		> <scaffold name>:<the start of the intergenic region>:<the end of the intergenic region>

	Example:
		>scaffold1_2_MH0001:1:692
		CACGTCGTCAGGAAGCTGACCCAGCTCGTGGTTGGCCTCGGCAGCAGCGAGCTTCACGTATGCAA
		TGGCCTTGACATACTCGGGATAGTCGCACATGTGCTTGCCCGAAATCTTGTAATTGTTGATG
		GCGCGCTGTGTCTGCACACCATAGTAAGCCTCGGCGGGTACCTGGAGCTCACCCAAGAGGTCGCT...
		
	Args: 
		gene_locations: An array that holds the lines from the annotation file

	Return:
		This method doesn't return anything; however, this method will write out files
		that are titled: Pruned_<Patient File>.seq.fa.

	"""

	#Convert set into list and sort it to match with global array
	sorted_file_list = sorted(list(file_list))
	annotation_index = 0

	#Scaffold pattern to match 
	scaffold_pattern_object = re.compile(scaffold_pattern, re.IGNORECASE)
	region_pos_pattern_object = re.compile(region_pos_pattern)

	for file_name in sorted_file_list:
		#Open both the patient file and the output file
		f = open("../MetaGenome/Data/%s.seq.fa" % file_name, "r")
		pruned_file = open("Pruned_%s.seq.fa" % file_name, "w")

		#Keep user updated towards execution status
		print "Starting: ", file_name
		file_start = time.time()

		#Python's ghetto do-while 
		while True:
			# Grab the first two lines of the patient file
			header = f.readline().strip()
			sequence = f.readline().strip()
			sequence_len = len(sequence)
			#If at the end of the file break out of while loop
			if header == "" or sequence == "":
				break

			total_lines = list()
			scaffold_match = scaffold_pattern_object.search(header)

			#If the annotation line has scaffold###_####
			if re.search(scaffold_match.group(0),gene_locations[annotation_index]):

				#Annotation file can hold locations on the same scaffold
				#Grab all relevant lines
				while True:
					if re.search(scaffold_match.group(0), gene_locations[annotation_index]):
						region_pos_match = region_pos_pattern_object.search(gene_locations[annotation_index])
						accession_number_match = accession_number_pattern_object.search(gene_locations[annotation_index])
						total_lines.append([accession_number_match.group(0), region_pos_match.group(1), region_pos_match.group(2), region_pos_match.group(3)])
						annotation_index += 1
					else:
						break

				#Cycle through each match
				for indicies in total_lines:
					number = indicies[0]
					start = int(indicies[1])
					end = int(indicies[2])
					strand = indicies[3]

					#If the original strand
					if strand == "+":
						#If the gene region doesn't start at the beginning
						if start != 1:
							if (start - 300) < 0:
								pruned_file.write("\n%s:%d:%d:%s:%s\n" % (header,1,start+50, strand, number))
								pruned_file.write("%s" % sequence[0:start+49])
							else:
								pruned_file.write("\n%s:%d:%d:%s:%s\n" % (header,start-300,start+50, strand, number))
								pruned_file.write("%s" % sequence[start-301:start+50])

					#The compliment strand
					#start -> end | end-> start
					else:
						#If the beginning of the region isn't at the end of the sequence
						if end != sequence_len:
							if (end + 300) > sequence_len:
								pruned_file.write("\n%s:%d:%d:%s:%s\n" % (header,end-50,sequence_len, strand, number))
								pruned_file.write("%s" % sequence[end-50:sequence_len])
							else:
								pruned_file.write("\n%s:%d:%d:%s:%s\n" % (header,end-50,end+300, strand, number))
								pruned_file.write("%s" % sequence[end-51:end+300])

		#Close both files
		f.close()
		pruned_file.close()
		file_end = time.time()
		print "Patient File Parsing Time: ", (file_end-file_start)

#main 
start = time.time()
locations = preprocess_file(gene_prediction)
end = time.time()
print "Annotation Parsing Time: ", (end-start)

start = time.time()
prune_files(locations)
end = time.time()
print "Total Process Time: ", (end-start)