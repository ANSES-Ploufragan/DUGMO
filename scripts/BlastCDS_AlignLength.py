#!/usr/bin/env python
import re

def Blast_lengthAlignment(fileTSV, fileEntete):
	fileEntete = open(fileEntete, "a")
	for line in open(fileTSV, 'r'):
		liste=line.split('\t')
		pident=liste[2]
		lengthAlign=liste[3]
		qlen=liste[4]
		slen=liste[5]
		qcovhsp=liste[6]
		if (float(pident)>=98) or (float(qcovhsp)>=98):
			Pourcentq=float(int(qlen)*(15/100))
			Pourcents=float(int(slen)*(15/100))
			if ((int(qlen)-Pourcentq) <= int(lengthAlign) <= (int(qlen)+Pourcentq)) and ((int(slen)-Pourcents) <= int(lengthAlign) <= (int(slen)+Pourcents)) :
				fileEntete.write(liste[0] + '\n')

Blast_lengthAlignment(snakemake.input[0], snakemake.output[0])
