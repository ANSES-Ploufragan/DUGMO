#!/usr/bin/env python
import re

def Blast_lengthAlignment(fileTSV, fileEntete):
	fileEntete = open(fileEntete, "a")
	for line in open(fileTSV, 'r'):
		liste=line.split('\t')
		lengthAlign=liste[2]
		qlen=liste[3]
		Pourcent=float(int(qlen)*(15/100))
		if (int(qlen)-Pourcent) <= int(lengthAlign) <= (int(qlen)+Pourcent):
			fileEntete.write(liste[0] + '\n')

Blast_lengthAlignment(snakemake.input[0], snakemake.output[0])
