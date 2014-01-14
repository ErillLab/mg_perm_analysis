
"""
@author: David 

This script parses the BGI_GeneSet20090523_annotation file and grabs the intergenic regions
from patient files. 

Dependencies:
	- Numpy

Assumes that files to be parsed are in lexicographic order.
"""

import re, sys, os, time
import numpy as np

### Global Settings ###
#Name of the file that holds the gene annotations
gene_prediction = "Metahit/BGI_GeneSet20090523_annotation"
#Regex pattern that matches the patient files
filepattern = "MH[\d]+"
#Regex pattern that matches the scaffold/contig lines in the annotation and patient files
scaffoldpattern = "(?:scaffold|c)[\d]+[^\d][\d]+"
#Regex pattern that matches the start and end positions in the annotation file
genepattern = ":[\d]+|:[\d]+"
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
				gene_locations[index] = line
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

	#Compile patterns to speed up search process
	gene_pattern_object = re.compile(genepattern, re.IGNORECASE)
	scaffold_pattern_object = re.compile(scaffoldpattern, re.IGNORECASE)

	#A checker to verify that the program is working with one file at a time.
	file_checker = ""

	#Convert set into list and sort it to match with global array
	sorted_file_list = sorted(list(file_list))
	annotation_index = 0


	for file_name in sorted_file_list:
		#Open both the patient file and the output file
		f = open("MetaHit/Data/%s.seq.fa" % file_name, "r")
		pruned_file = open("MetaHit/Pruned/Pruned_%s.seq.fa" % file_name, "w")

		#Keep user updated towards execution status
		print "Starting: ", file_name
		file_start = time.time()

		while True:
			# Grab the first two lines of the patient file
			header = f.readline().strip()
			sequence = f.readline().strip()

			#If at the end of the file break out of while loop
			if header == "" or sequence == "":
				break

			total_lines = list()
			scaffold_match = scaffold_pattern_object.search(header)

			#If the annotation line has scaffold###_####
			if re.search(scaffold_match.group(0),gene_locations[annotation_index]):

				#Annotation file can hold locations on the same scaffold
				#Grab all relevant lines
				while re.search(scaffold_match.group(0),gene_locations[annotation_index]):
					total_lines.append(gene_locations[annotation_index])
					annotation_index += 1

				#Grab the start and stop positions 
				#Creates a 2d array: [[":2", ":100"], [":1", ":50"] ...]
				match = map(gene_pattern_object.findall, total_lines)

				#If only one line for a particular sequence
				if len(match) == 1:
					start = int(match[0][0][1:])
					end = int(match[0][1][1:])

					# Write out the intergenic region
					# If region is not at start of sequence then write from beginning to start
					if start != 1:
						pruned_file.write("%s:%d:%d\n" % (header,1,start))
						pruned_file.write("%s\n" % sequence[0:start-1])

						#If the reqion does finish at the end of the sequence, then write from end of region to end of sequence
						if end != len(sequence):
							pruned_file.write("%s:%d:%d\n" % (header,end,len(sequence)))
							pruned_file.write("%s\n" % sequence[end:len(sequence)])

					#If the gene region starts of at the beginning of the sequence just write from end of region to end of sequence
					else:
						pruned_file.write("%s:%d:%d\n" % (header,end,len(sequence)))
						pruned_file.write("%s\n" % sequence[end:len(sequence)])
				else:
					#If there are more then one line
					begin = 1
					for indicies in match:
						start = int(indicies[0][1:])
						end = int(indicies[1][1:])

						#Grab the intergenic regions using similar method specified above in the if len(match) section
						if begin != start:
							pruned_file.write("%s:%d:%d\n" % (header,begin,start))
							pruned_file.write("%s\n" % sequence[begin-1:start])
							begin = end
						else:
							begin = end
					if end != len(sequence):
						pruned_file.write("%s:%d:%d\n" % (header,end,len(sequence)))
						pruned_file.write("%s\n" % sequence[end:len(sequence)])

			#If the scaffold isn't mentioned in the annotation file just write it out into the file
			else:
				pruned_file.write("%s:%d:%d\n" % (header,1,len(sequence)))
				pruned_file.write("%s\n" % sequence)

		#Close both files
		f.close()
		pruned_file.close()
		file_end = time.time()
		print "Patient File ParsingTime: ", (file_end-file_start)
	

#main 
start = time.time()
locations = preprocess_file(gene_prediction)
end = time.time()
print "File Parsing Time: ", (end-start)

start = time.time()
prune_files(locations)
end = time.time()
print "Total Process Time: ", (end-start)