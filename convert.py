# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 15:52:05 2014

@author: David Nicholson
"""
import glob
import mysql.connector as db

#When writing the sequence to the file take the reverse compliment
#work on implementing our found data with the metahit databse.
def parseStrandScores(strand, scores,patient_file, db_cursor):
	metahit_path = "./MetaHit/Pruned"
	start = 0
	end = 0
	with open(metahit_path + patient_file, "r") as g:
		for line in g:
			if ">" in line:
				header_attr = line[1:].split("|")
				start = int(header_attr[2])
				end = int(header_attr[3])
			else:
				pass
	
	
def main():
	strand = []
	is_scores = True
	scores = []
	con = db.connect(user = 'root', password = 'root', host = 'localhost', database = 'metagenome')
	results_metahit_db = glob.glob("results_MH[0-9]+")
	for patient_file in results_metahit_db:
		f = open(patient_file, "r")
		for line in f:
			if ">" in line:
				if "Strand" in line:
					is_scores = False
			else:
				if is_scores:
					scores = map(lambda x: float(x), line.split(" "))
				else:
					strand = map(lambda x: float(x), line.split(" "))
		parseStrandScores(strand, scores, patient_file,con.cursor())
					