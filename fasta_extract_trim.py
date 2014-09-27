#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

fasta_file = sys.argv[1]  # Input fasta file
extract_file = sys.argv[2] # Input wanted file, one gene name per line
result_file = sys.argv[3] # Output fasta file

extract = {}
with open(extract_file) as f:
    for line in f:
        seq, nada, start, end, length = line.strip().split("\t")
        start = int(start)
        end = int(end)
        length = int(length)
        if start != 1 or end != length:
            extract[seq] = (start, end)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
f =  open(result_file, "w")
seqs = []
for seq in fasta_sequences:
        nuc = seq.id
	if nuc in extract :
		seqs.append(seq[extract[nuc][0]-1: extract[nuc][1]-1])
		if len(seqs) > 1000:
			SeqIO.write(seqs, f, "fasta")
			seqs = []
SeqIO.write(seqs, f, "fasta")




