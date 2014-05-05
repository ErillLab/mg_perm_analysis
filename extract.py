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

at_lab = False
#File regex pattern to match the metahit gene annotation format
metahit_pattern = "(?P<accession>GL[\d]+)[^\d]+(?P<source>MH0086)[\S]+(?P<scaffold>(?:scaffold|C)[\d]+_[\d]+):(?P<start>[\d]+):(?P<end>[\d]+):(?P<strand>[+-])(?P<cog>(((eu|me)?[KNCOG]+[\d]+[;]?)|[NA]+|\s)+)"
#Name of the gene annotation file
if at_lab:
	gene_annotation = "MetaHit/BGI_GeneSet20090523_annotation"
else:
	gene_annotation = "/Users/Dave/Documents/Erill Lab/Metagenome/BGI_GeneSet20090523_annotation"

#The directory that holds all the patient files
if at_lab:
	directory = "MetaHit/Data/"
else:
	directory = "/Users/Dave/Documents/Erill Lab/Metagenome/Data/"

#Grab the names of files that end in .seq.fa
files = glob.glob(directory + "*.seq.fa")

#main
def main():
	
	start = time.time()
	
	#parse the annotation file
	annotations = pr.load_annotations(gene_annotation,metahit_pattern)

	an_end = time.time()
	print "Annotation File: ", (an_end-start)

	if len(sys.argv) > 3:
		pr.parse_files(annotations, files, int(sys.argv[1]), int(sys.argv[2]))
	else:
		#create the non-coding region files
		pr.parse_files(annotations, [directory+"MH0086.seq.fa"])

	end = time.time()
	print "Parsing Time: ", (end-an_end)
	print "Total Time: ", (end-start)
main()