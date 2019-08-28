#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Counting on each sequence for a given number of words et ne garde que les séquences supérieures à 27 nucléotides
#Params
#MEFfile = Formatted input file (one sequence by line)
#EnteteFile = Entete Sequence input file
#freqGRFile = frequence reference genome input file
#distFile = Output file with distance calculation

import re, numpy as np, pandas as pd
import ahocorasick as ahc
import scipy.spatial.distance
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

# Create Aho-Corasick automaton
def make_aho_automaton(keywords):
    A = ahc.Automaton()
    for (key, cat) in keywords:
            A.add_word(key, (cat, key)) 
    A.make_automaton() 
    return A

def calcul_dist_prop(MEFfile, EnteteFile, freqGRFile, distFile):
	#Initialization of variables and reading of files
	freqGenomeRef = pd.read_csv(freqGRFile, sep="\t", header=0)
	EnteteSeq = pd.read_csv(EnteteFile, sep="\t", header=None)
	tuples=[]
	allmots = []
	#Creation of tuples with the words to search for
	for elt in freqGenomeRef['mot']:
		mot=elt[:-1].replace(" ", "")
		key=(mot.upper(), 0)
		tuples.append(key)
		allmots.append(mot.upper())
	#Length of a word
	wordLen = len(allmots[0])	
	#Creation of the automaton with tuples
	A = make_aho_automaton(tuples)
	Result = pd.DataFrame(columns=['BrayCurtis', 'Length', 'MeanScore', 'DensityNuc', 'GCpourcent'], dtype=np.float64)
	#Search for words in each sequence (counting)
	for line in open(MEFfile, 'r'):
		find = []
		found_keywords = find + allmots		
		for start_end, (cat, keyw) in A.iter(line):
			if keyw in found_keywords:
				indice = found_keywords.index(keyw)
				found_keywords[indice]=line.count(keyw)	
		for a in range(0, len(found_keywords)):
			if type(found_keywords[a]) != int:
				found_keywords[a]=0
		#Sum of the counts 
		sumcomptageCDScible = sum(found_keywords)		
		#Convert in dataframe format
		found_keywords = pd.DataFrame(data=found_keywords)
		found_keywords.columns = ['found_keywords']
		#Length of the sequence
		lenSeq = len(line)-1 
		#Counting density per nucleotide
		CompNuc = sumcomptageCDScible/lenSeq 
		#Mean of the scores
		TabComptageScore = pd.concat([found_keywords, freqGenomeRef['score']], axis=1)
		TabComptageScore.columns = ['comptage', 'score']
		TabComptageScore = TabComptageScore.replace([np.inf, -np.inf], 600) #600=maximum score -> to be modified 	
		MeanScore = sum(TabComptageScore['comptage'] * TabComptageScore['score'])/len(TabComptageScore.index)	
		#Calculation of the GC percentage
		countGC = line.count('C') + line.count('G')
		pGC = (countGC/lenSeq)*100
		#Transformation in proportions
		comptageall = pd.concat([found_keywords.reset_index(drop=True), freqGenomeRef['count']], axis=1)
		#If no word found then summation is equal to 0 (division by zero impossible) -> maximum value for distance calculation
		if sumcomptageCDScible == 0:
			BC = 1.0
		else :
			comptagetrouve = comptageall
			sumcomptageCSDglobal = sum(comptagetrouve['count'])	
			proportion = comptagetrouve['found_keywords'].apply(lambda x: x*sumcomptageCSDglobal/sumcomptageCDScible)
			proportionall = np.transpose(pd.concat([proportion, comptagetrouve['count']], axis=1))
			#Calculation of the Bray Curtis distance
			BC = scipy.spatial.distance.braycurtis(proportionall.iloc[0], proportionall.iloc[1])
		#Results table with all distance values
		Result = Result.append({'BrayCurtis': BC, 'Length': lenSeq, 'MeanScore': MeanScore, 'DensityNuc': CompNuc, 'GCpourcent': pGC}, ignore_index=True)
	#Concatenate the sequence headers and distance calculation results
	Dist = pd.concat([EnteteSeq.reset_index(drop=True), round(Result, 6)], axis=1)
	np.savetxt(distFile, Dist, fmt='%s', header="Id_Sequence,BrayCurtis,Length,MeanScore,DensityNuc,GCpourcent", delimiter=',', comments='')

calcul_dist_prop(snakemake.input[0], snakemake.input[1], snakemake.input[2], snakemake.output[0])
