#!/usr/bin/env python
import re, numpy as np, sys, argparse, pandas as pd, os
import warnings
warnings.filterwarnings("ignore")

#### Deletion of parts of GMO sequences aligned with the reference genome
#### Creating a fasta file without the parts that align

#Params
#GMOdb = Input file in tsv format
#BlastAlignFile = Input file with informations of alignment in BlastDB format
#Output = Output file in fasta format

def DeleteSeqGMOAlignGenomRef (GMOdb, BlastAlignFile, Output):
	
	if os.stat(BlastAlignFile).st_size == 0:
		fileOutput = open(Output, "a")
		fileGMOdb = open(GMOdb, "r")
		fileOutput.write(fileGMOdb.read())
		fileOutput.close()
		fileGMOdb.close()	
	else :
		#Traitement 
		BlastAlign = pd.read_csv(BlastAlignFile, sep="\t", header=None)
		BlastAlign.columns = ['qseqid', 'sseqid', 'pident', 'length', 'sstart', 'send', 'evalue', 'sstrand'] 
		#Deletion of lines containing the same positions with respect to a given GMOdb sequence ID
		alignTable = BlastAlign.sort_values(by='sseqid').drop_duplicates(subset=['sseqid', 'sstart', 'send'], keep='first')
		DeleteSeqNameCount = alignTable['sseqid'].value_counts()
		DeleteSeqPos = pd.DataFrame(columns=['Seq_Id', 'PosStart', 'PosEnd'], dtype=np.int)	
		#Table path containing the number of lines found per sequence
		for i in range (0, len(DeleteSeqNameCount)):
			seqID = DeleteSeqNameCount.index[i]	
			seqIDall = alignTable.loc[alignTable['sseqid']==seqID,:]
			if DeleteSeqNameCount.iloc[i] == 1:
				if seqIDall.iloc[0, 4] > seqIDall.iloc[0, 5]:
					temp = seqIDall.iloc[0, 4]
					seqIDall.iloc[0, 4] = seqIDall.iloc[0, 5]
					seqIDall.iloc[0, 5] = temp
				DeleteSeqPos = DeleteSeqPos.append({'Seq_Id': seqIDall.iloc[0, 1], 'PosStart': seqIDall.iloc[0, 4], 'PosEnd': seqIDall.iloc[0, 5]}, ignore_index=True)
			else: 
				#Sort alignment lengths in descending order
				seqIDallSort = seqIDall.sort_values(by='length', ascending=False)
				#Direction of the strand: the positions must always be from the smallest value to the largest value  
				for k in range(0, len(seqIDallSort)):
					#Change the sstart and send columns so that the smallest value is in the sstart
					if seqIDallSort.iloc[k, 4] > seqIDallSort.iloc[k, 5]:
						temp = seqIDallSort.iloc[k, 4]
						seqIDallSort.iloc[k, 4] = seqIDallSort.iloc[k, 5]
						seqIDallSort.iloc[k, 5] = temp
				#Sort alignments by sstart in ascending orderS
				seqIDallSort=seqIDallSort.sort_values(by='sstart')
				#Set min and max variables for both positions with the longest alignment length					
				minAlign = seqIDallSort.iloc[0, 4]
				maxAlign = seqIDallSort.iloc[0, 5]
				for j in range(0, len(seqIDallSort)-1):
					#If the next positions are not between the two variables min and max then add in a table the positions min and max
					#edit new variables min and max with the new values
					if (minAlign >= seqIDallSort.iloc[j+1, 4] and maxAlign >= seqIDallSort.iloc[j+1, 5]):
						if (seqIDallSort.iloc[j+1, 4] <= minAlign <= seqIDallSort.iloc[j+1, 5]):
							minAlign = seqIDallSort.iloc[j+1, 4]
						else:
							DeleteSeqPos = DeleteSeqPos.append({'Seq_Id': seqIDallSort.iloc[j, 1], 'PosStart': minAlign, 'PosEnd': maxAlign}, ignore_index=True)
							minAlign = seqIDallSort.iloc[j+1, 4]
							maxAlign = seqIDallSort.iloc[j+1, 5]	
					elif (minAlign >= seqIDallSort.iloc[j+1, 4] and maxAlign <= seqIDallSort.iloc[j+1, 5]):
						DeleteSeqPos = DeleteSeqPos.append({'Seq_Id': seqIDallSort.iloc[j, 1], 'PosStart': minAlign, 'PosEnd': maxAlign}, ignore_index=True)
						minAlign = seqIDallSort.iloc[j+1, 4]
						maxAlign = seqIDallSort.iloc[j+1, 5]	
					elif (minAlign <= seqIDallSort.iloc[j+1, 4] and maxAlign <= seqIDallSort.iloc[j+1, 5]):
						if (seqIDallSort.iloc[j+1, 4] <= maxAlign <= seqIDallSort.iloc[j+1, 5]):
							maxAlign = seqIDallSort.iloc[j+1, 5]
						else:
							DeleteSeqPos = DeleteSeqPos.append({'Seq_Id': seqIDallSort.iloc[j, 1], 'PosStart': minAlign, 'PosEnd': maxAlign}, ignore_index=True)
							minAlign = seqIDallSort.iloc[j+1, 4]
							maxAlign = seqIDallSort.iloc[j+1, 5]
				#Table with the positions of the sequence pieces to be removed in the GMO db sequences
				DeleteSeqPos = DeleteSeqPos.append({'Seq_Id': seqIDallSort.iloc[j, 1], 'PosStart': minAlign, 'PosEnd': maxAlign}, ignore_index=True)			
		#Removal of parts of sequences aligned with the ref genome
		with open(GMOdb, "r") as f:	
			Seq=re.split('^>', f.read(), flags=re.MULTILINE)		
			j=0
			while j < len(DeleteSeqPos):
				#Searching for sequence names to be deleted  				
				i=0
				for i, item in enumerate(Seq):
					if DeleteSeqPos.iloc[j,0] in item:
						idx=i				
				#Table with all positions to be deleted for a given sequence
				MultiPosDelete = DeleteSeqPos.loc[DeleteSeqPos['Seq_Id']==DeleteSeqPos.iloc[j,0],:]
				#Processing to delete the \n and select the header (for which it is not included in the position count)
				Seq[idx] = Seq[idx].replace('\n','')
				trouve = re.search('[ATCGN]{6}|[atcgn]{6}', Seq[idx])
				entete = Seq[idx][0:trouve.start()]
				sequence = Seq[idx][trouve.start():len(Seq[idx])]
				#Deletion of all alignments with their corresponding positions for a given sequence
				for k in range (0, len(MultiPosDelete)):
					if int(MultiPosDelete.iloc[k,1]) == 1 :					
						seqNew = sequence[int(MultiPosDelete.iloc[k,2]):len(Seq[idx])]
					else :
						seqNew = sequence[0:int(MultiPosDelete.iloc[k,1])] + sequence[int(MultiPosDelete.iloc[k,2]):len(Seq[idx])]	
					if k < (len(MultiPosDelete)-1):
						sequence = seqNew
						posStart = int(MultiPosDelete.iloc[k,1]) + (int(MultiPosDelete.iloc[k+1,1]) - int(MultiPosDelete.iloc[k,2]))					
						posEnd = posStart + (int(MultiPosDelete.iloc[k+1,2]) - int(MultiPosDelete.iloc[k+1,1]))	
						MultiPosDelete.iloc[k+1,1] = posStart
						MultiPosDelete.iloc[k+1,2] = posEnd					
					k=k+1
				#Replacement by the new sequence without the alignment portions to the format in the list
				if seqNew == "" or len(seqNew) < 27 :
					Seq[idx] = ''
				else : 				
					seqNewfasta = re.findall('.{1,60}', ''.join(list(seqNew)))
					Seq[idx] = entete.replace(',','-') + '\n' + '\n'.join(seqNewfasta) + '\n'
				j=j+len(MultiPosDelete)
		#Writing new sequences (without the parts aligned with the ref genome) in the output file
		fichier = open(Output, "a")
		for i in range(1, len(Seq)):
			if Seq[i] != '':
				fichier.write('>' + Seq[i]) 
			i=i+1
		#Closing files
		fichier.close()
		f.close()

DeleteSeqGMOAlignGenomRef (snakemake.input[0], snakemake.input[1], snakemake.output[0])
