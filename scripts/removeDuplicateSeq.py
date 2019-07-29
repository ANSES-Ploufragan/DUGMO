#!/usr/bin/env python
from itertools import groupby

#Remove duplicate sequences in fasta file

def RemoveDuplicateSeq (fileinput, fileoutput):
	ishead = lambda x: x.startswith('>')
	all_seqs = set()
	with open(fileinput, 'r') as handle:
		with open(fileoutput, 'w') as outhandle:
			head = None
			for h, lines in groupby(handle, ishead):
				if h:
					head = lines.__next__()
				else:
					seq = ''.join(lines)
					if seq not in all_seqs:
						all_seqs.add(seq)
						outhandle.write('%s%s' % (head, seq))
	
RemoveDuplicateSeq(snakemake.input[0], snakemake.output[0])
