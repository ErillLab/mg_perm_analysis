# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 15:52:05 2014

@author: David Nicholson
"""
import glob
import mysql.connector as db

#When writing the sequence to the file take the reverse compliment
#work on implementing our found data with the metahit databse.
#the length of motif is 16 
# the cog information is 7th pos
def parseStrandScores(strand, scores,patient_file, db_cursor, file_write):
	metahit_path = "./MetaHit/Pruned"
	header_attr = []
	index = 0
	sequence = ""
	#open the metahit patient file
	with open(metahit_path + patient_file[8:] + ".seq.fa", "r") as g:
		for line in g:
			if ">" in line:
				header_attr = line[1:].split("|")
			else:
				#gather the necssary information from the database
				#TODO gather the taxnomic field headers
				db_cursor.execute("SELECT * FROM joe_gene_ortho_mhcogs WHERE GENE_COG_NAME = %s LIMIT = 1", (header_attr[6]))
				all_rows = db_cursor.fetchall()
				if len(all_rows) != 0:
					#cycle through the long array of scores and get sequences pertaining to the first scaffold
					while((index%300) < 284):
						#if higher score is on the original direction
						if strand[index]:
							header_attr[4] = "+"
						else:
							header_attr[4] = "-"
						#grab the sequence from the patneit data
						sequence = line[index:index+16]
						#write the header and the sequence
						file_write.write("|".join(header_attr) + "|". join(all_rows) + "\n")
						file_write.write(sequence +"\n")
						index = index + 1
					#make the index point to the beginning of the new sequence
					index = index + 16
				else:
					#make the index point to the beginning of the new sequence
					index = index + 300
		file_write.close()
		db_cursor.close()

def main():
	strand = []
	is_scores = True
	scores = []
	#connect to the database
	con = db.connect(user = 'root', password = 'root', host = 'localhost', database = 'metagenome')
	results_metahit_db = glob.glob("results_MH[0-9]+")

	#for each patient file result convert the numbers into two seperate arrays and then write out each strand.
	for patient_file in results_metahit_db:
		with open(patient_file, "r") as f:
			for line in f:
				if ">" in line:
					if "Strand" in line:
						is_scores = False
				else:
					if is_scores:
						scores = map(lambda x: float(x), line.split(" "))
					else:
						strand = map(lambda x: True if x == 0 else False, line.split(" "))
			parseStrandScores(strand, scores, patient_file,con.cursor(), open("Regulatory_" + patient_file[8:] + ".seq.fa", "w"))
		is_scores = True
					