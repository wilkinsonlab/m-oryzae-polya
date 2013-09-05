from Bio import SeqIO
import sys
import os.path

infile = open(sys.argv[1], 'r')
for record in SeqIO.parse(infile, "fasta"):
    SeqIO.write(record, record.id + ".fasta", "fasta")
