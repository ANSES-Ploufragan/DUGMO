#!/usr/bin/env python
import re

#Params
#fileSeq = Input file in fasta to formating
#fileMEF = Formatted output file (one sequence by line) with position 3 of each codon
#fileEntete = Output file with sequence name 
#filePos3codon = Output file with position 3 of each codon on fasta format
#fileMEFPos = Formatted output file (one sequence by line) with position 3 of each codon
#fileEntetePos = Output file with sequence name of seq with position 3 of each codon

def MEF_Entete_pos3codon(fileSeq, fileMEFPos, fileMEF, fileEntete, fileEntetePos, filePos3codon):
	fileMEF = open(fileMEF, "a")
	fileMEFPos = open(fileMEFPos, "a")
	fileEntete = open(fileEntete, "a")
	fileEntetePos = open(fileEntetePos, "a")
	filePos3codon = open(filePos3codon, "a")	
	with open(fileSeq, 'r') as f:
		EnteteSeq=re.split('^>', f.read(), flags=re.MULTILINE)	
		i=1
		while i<len(EnteteSeq):	
			Seqentete=EnteteSeq[i].split('\n')
			entete=Seqentete[0].replace('>','')
			seq=''.join(Seqentete[1:len(Seqentete)])
			fileEntete.write('>' + entete.replace(',','-') + '\n')
			fileMEF.write(seq.upper() + '\n')
			if len(seq) >= 27:
				seq = seq.upper()
				seqlist = list(seq)
				seqlist = seqlist[2::3]
				fileMEFPos.write(''.join(seqlist) + '\n')
				fileEntetePos.write('>' + entete.replace(',','-') + '\n')
				seqfasta = re.findall('.{1,60}', ''.join(seqlist))
				filePos3codon.write('>' + entete.replace(',','-') + '\n' + '\n'.join(seqfasta) + '\n')
			i=i+1
	fileMEF.close()
	fileMEFPos.close() 
	fileEntete.close()
	fileEntetePos.close()
	filePos3codon.close() 
	f.close()

MEF_Entete_pos3codon(snakemake.input[0], snakemake.output[0], snakemake.output[1], snakemake.output[2], snakemake.output[3], snakemake.output[4])
