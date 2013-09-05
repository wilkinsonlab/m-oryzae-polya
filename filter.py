import sys, re
import os
from Bio import SeqIO
import pysam


infile = pysam.Samfile(sys.argv[1], "rb")
fasta_file = sys.argv[2]
offset = int(sys.argv[3])
outfile = sys.argv[4]

oligoT = 16

tempname = "_" + sys.argv[1] + ".filter.temp"
tempfile = pysam.Samfile(tempname, "wb", template=infile)

fasta_seqs = {}
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    fasta_seqs[seq_record.id] = seq_record
    
total = 0
mapped = 0
fake = 0
for read in infile.fetch():
    total += 1
    if not read.is_unmapped and read.mapq >= 30 and ((read.seq.count('A') + read.seq.count('T')) / float(len(read.seq)) < 0.80):
        chrx = infile.getrname(read.tid)
        if read.is_reverse:
            pos = read.pos + offset
            seq = str(fasta_seqs[chrx].seq[pos + read.rlen: pos + read.rlen + oligoT])
        else:
            pos = read.pos - offset
            seq = str(fasta_seqs[chrx].seq[pos - oligoT: pos].reverse_complement())  
        if (seq.count('A') + seq.count('G')) / float(oligoT) < 0.9:# and not re.search("GGG+", seq)):
            tempfile.write(read)
            mapped += 1
        else:
            fake += 1
            
pysam.sort(tempname, outfile)
pysam.index(outfile + ".bam")
os.remove(tempname)
print "Mapped " + str(mapped) + " out of " + str(total) + " reads,", round(mapped / float(total) * 100, 2), "%", "fakes:", round(fake / float(total) * 100, 2), "%"


tempfile.close()
infile.close()
