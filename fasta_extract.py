#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

fasta_file = sys.argv[1]  # Input fasta file
extract_file = sys.argv[2] # Input wanted file, one gene name per line
result_file = sys.argv[3] # Output fasta file

extract = set()
with open(extract_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            extract.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
f =  open(result_file, "w")
seqs = []
i=0
for seq in fasta_sequences:
        i+=1
        if i % 10000 == 0:
            print i
        nuc = seq.id
	if nuc in extract :
		seqs.append(seq)
		if len(seqs) > 10000:
			SeqIO.write(seqs, f, "fasta")
			seqs = []
SeqIO.write(seqs, f, "fasta")




