"""
Created February 19,2014
@author David Nicholson
This file is the backend of the filter process. It relies on parameters being passed in from extract.py
"""
import re, os, sys, time
import numpy as np

def load_annotations(filename, pattern, filetype= "FASTA"):
	"""
	This function will parse a gene annotation file.
	The parse data will be in the format: [scaffold,source,start, end, strand, accession#]

	The function will require a passed in regular expression.
	A tutorial for beginners to regular expressions: http://www.regular-expressions.info/tutorial.html
	To get a more pratical feel for regular-expressions (requires java): http://regex.powertoy.org/
	Python's Regular Expression: http://docs.python.org/2/howto/regex.html
	
	Each passed in pattern needs to use the symbolic grouping method. (?P<name>...)
	This allows for groupdict function to be called and guarentees that the annotation array will contain
	all necessary information.


	Arguments:
		filename: name of the gene annotation file
		pattern: a regex pattern to search files
		filetype: specifies how the annotation file is formatted

	Return:
		annotations: a 2d array that contains essential data
			Ex. [[scaffold,source, start, end, strand, accession#,cog1, cog2 ...],[scaffold, source, start, end, strand, accession#, cog1, cog2 ...]]

	Example:
		>>> annotations = load_annotations("GENE FILE.fa", "FASTA", "(?P<scaffold>scaffold):(?P<source>MH[\d]+):(?<start>[\d]+):(?<end>[\d]+):(?<strand>[+|-]):(?<accession>accession)(?<cog>cog)")
		>>> annotations
			[["scaffold01", "MH00001", "1", "100", "+", "accession", "cog"], ["scaffold02", "MH00001", "300", "900", "-","accession", "cog"]] ...]
	"""
	#compile the regular expression to speed up search process
	pattern_object = re.compile(pattern)

	#For now one type of format: FASTA. Expecting to add code to parse other file formats
	if filetype == "FASTA":
		#open file
		with open(filename, "r") as g:
			#pre allocate array to hold [scaffold, start, end, strand, accession #]
			g.seek(0,os.SEEK_END)
			pos = g.tell()
			gene_locations = np.empty(shape=(pos,1),dtype=object)
			g.seek(0,os.SEEK_SET)
			index = 0

			#run through line by line of the annotation file
			for line in g:

				#Regex match to gather essential data 
				match = pattern_object.search(line)

				#ignore lines other than gene annotation lines. Ignore MH0006 because it causes problems
				if match and "MH0006" not in line:
					#insert essential data into the array
					terms = match.groupdict()
					gene_locations[index] = [checkStrand(terms)]
					index += 1

		return gene_locations

def checkStrand(entries):
	"""
	This function checks to see if it can find all the necessary information within an annotation line.
	If it can find information it will attempt to add information where appropiate.
	
	Hopefully an entry will contain a dictionary in the following format:
		{"scaffold": "scaffold14_9, "source":"MH0006", "start":"300", "end":"500", "strand":"+", "accession":"GL50000","cog":"COG09123 K0129398"}

	Arguments:
		a line from a gene annotation file in dictionary format
	
	Return:
		an array of all the essential data in a specified order

	Example:
		>>>> entry = checkStrand({"scaffold": "scaffold14_9", "source":"MH0006", "start":"300", "end":"500", "strand":"+", "accession":"GL50000", "cog":"COG09123 K0129398"})
		>>>> entry
			["scaffold14_9", "MH0006", "300", "500" , "+", "GL50000", "COG09123", "K0129398"]
	"""
	#pre-define arrays to hold elements
	missing_terms = []
	important_missing_terms = []
	return_list = []

	#gather necessary key information
	keys = ["scaffold", "source", "start", "end", "strand", "accession", "cog"]
	for key in keys:
		if key not in entries:
			#if the program doesn't have information: scaffold, start, end it needs to break
			if key != "source" and key != "accession" and key != "strand":
				important_missing_terms.append(key)
			missing_terms.append(key)

	#if missing important information break it
	if "scaffold" in  missing_terms or "start" in missing_terms or "end" in missing_terms:
		print "Please make sure to have the necessary information:"
		print important_missing_terms
		sys.exit(1)

	#guarenteed to have a scaffold
	return_list.append(entries["scaffold"])
	
	#if source is not missing then put it into the return list
	if "source" not in missing_terms:
		return_list.append(entries["source"])

	#if missing strand (+-) ususally it is in the terms of end < start for compliment strand
	if "strand" in missing_terms:
		start = int(entries['start'])
		end = int(entries['end'])
		if start < end:
			return_list.append(entries["start"])
			return_list.append(entries["end"])
			return_list.append("+")
		else:
			return_list.append(entries["end"])
			return_list.append(entries["start"])
			return_list.append("-")
	else:
		return_list.append(entries["start"])
		return_list.append(entries["end"])
		return_list.append(entries["strand"])

	#if accession is in the list
	if "accession" not in missing_terms:
		return_list.append(entries["accession"])

	if "cog" not in missing_terms:
		#add listing of all the cog information
		cogs = re.split("\s|;", entries["cog"])
		for cog in cogs:
			if cog != "":
				return_list.append(cog)

	return return_list
	

def parse_files(annotations,filenames, upstream = 300, downstream = 50):
	"""
	This function parses out the non-coding regions of the metagenome files. 
	The ouput is a pruned file that contains a fasta file with the regions.

	The file will contain:
		>scaffold|source of the metagenome|start position|end position|strand|accession number
		ATGCACACACAGTCTTTTTGGGGGGGGGGCTTTAAACAGTACGTAGGAGGGAGATTAGTCA

	Arguments:
		annotations: a parsed array for the gene annotations
		filenames: the names of the files to be read in
		upstream: an optional variable to specify the number of desired base pairs upstream from the start codon
		downstream: an optional variable to specify the number of desired base pairs down stream from the start codon

	Returns:
		A newly created file with the non-coding regions. 

	Ex.
		>>>parse_files([["scaffold01", "MH001", "1", "100", "+", "accession", "cog"], ...], ["MH0001.seq.fa",...])
		>>> 
			*No output above^. A file should be generated in the same directory as this file.*
	"""
	index = 0

	#run through the avaiable metagenome files
	for metagenome_file in filenames:

		#grab name of file
		metagenome_file_header = os.path.splitext(os.path.splitext(os.path.basename(metagenome_file))[0])[0]

		#Update user about the progress
		print "Working on: ", metagenome_file_header
		start_time = time.time()

		#shortcut to get list of annotations that pertain to the specified file
		#not really a short cut need to figure out a better way to do this...
		#sub_annotations = [item for item in [annotation for annotation in annotations] if metagenome_file_header in item]
		
		#open file for writing
		g = open("MetaHit/Pruned/Pruned_%s.fa" % metagenome_file_header, "w")

		with open(metagenome_file, "r") as f:

			while True:
				#gather header and sequence
				line = f.readline().strip()
				seq = f.readline().strip()

				#if at end of the file break out of loop
				if line == "" or seq == "":
					break

				annotation_line = []
				#search for the header in the array of annotations
				#annotation_line = filter(lambda x: x if "%s_%s" % (x[0], x[1]) in line else None, sub_annotations)
				#annotation_line = [x for x in sub_annotations if "%s_%s" % (x[0], x[1]) in line]
				while(annotations[index][0] != None and "%s_%s" % (annotations[index][0][0],annotations[index][0][1]) in line):
					annotation_line.append(annotations[index][0])
					index = index + 1

				#if found:
				if len(annotation_line) > 0:
					
					#go through the annotated array
					for extract in annotation_line:
						
						#grab the start and the stop positions of the annotated entries
						start = int(extract[2])
						end = int(extract[3])

						#if the original strand
						if extract[4] == "+":
							if start - upstream >= 0:
								
								#verify that the written sequence is equal to the specified bounds (up,down)stream
								assert (upstream+downstream) == len(seq[start-upstream:start+downstream]), "%s\n%s" % (line,seq)

								#write out custom format for output
								extract[2] = str(start-upstream)
								extract[3] = str(start+downstream)
								g.write(">%s\n" % ("|".join(extract)))
								g.write("%s\n" % (seq[start-upstream:start+downstream]))

						#The reverse compliment strand
						else:
							
							if end + upstream <= len(seq):
								
								#verify that the written sequence is equal to the specified bounds (up,down)stream
								assert (upstream+downstream) == len(seq[end-downstream:upstream+end]), "%s\n%s" % (line,seq)

								#write out custom format for output
								extract[2] = str(end-downstream)
								extract[3] = str(end+upstream)
								g.write(">%s\n" %("|". join(extract)))
								g.write("%s\n" % seq[end-downstream:upstream+end])

		#Diagnositics for the User
		end_time = time.time()
		g.close()
		print "Time Taken: ", (end_time-start_time)
