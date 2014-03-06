"""
Created February 19, 2014
@author David Nicholson

This script uses prune.py to extract non-coding regions from metagenomics data.
Upstream and downstream positions, relative to the start of the called gene, can be specified via command line arguments.
These positions are optional the program can run with a default of 300bp upstream and 50bp downstream. 

Command is: python extracy.py <upstream start> <downstream end>
Ex. python extract.py 300 50
"""
import sys, time, glob
import prune as pr

#File regex pattern to match the metahit gene annotation format
metahit_pattern = "(?P<accession>GL[\d]+)[^\d]+(?P<source>MH[\d]+)[\S]+(?P<scaffold>(?:scaffold|c)[\d]+_[\d]+):(?P<start>[\d]+):(?P<end>[\d]+):(?P<strand>[+-])"
#Name of the gene annotation file
gene_annotation = "MetaGenome/BGI_GeneSet20090523_annotation"

#The directory that holds all the patient files
directory = "MetaGenome/Data/"

#Grab the names of files that end in .seq.fa
files = glob.glob(directory + "*.seq.fa")

#main
def main():
	
	start = time.time()
	
	#parse the annotation file
	if len(sys.argv) > 3:
		annotations = pr.load_annotations(gene_annotation,metahit_pattern, int(sys.argv[1]), int(sys.argv[2]))
	else:
		annotations = pr.load_annotations(gene_annotation,metahit_pattern)

	an_end = time.time()
	print "Annotation File: ", (an_end-start)

	#create the non-coding region files
	pr.parse_files(annotations, files)

	end = time.time()
	print "Parsing Time: ", (end-an_end)
	print "Total Time: ", (end-start)
main()