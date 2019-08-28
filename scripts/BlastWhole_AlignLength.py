#!/usr/bin/env python
import re

def Blast_lengthAlignment(fileTSV, fileEntete):
	fileEntete = open(fileEntete, "a")
	for line in open(fileTSV, 'r'):
		liste=line.split('\t')
		lengthAlign=liste[3]
		qlen=liste[4]
		Pourcentq=float(int(qlen)*(15/100))
		if ((int(qlen)-Pourcentq) <= int(lengthAlign) <= (int(qlen)+Pourcentq)) :
			fileEntete.write(liste[0] + '\n')

Blast_lengthAlignment(snakemake.input[0], snakemake.output[0])
