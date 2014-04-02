import sys, re
import os
from Bio import SeqIO
import pysam


infile = pysam.Samfile(sys.argv[1], "rb")
fasta_file = sys.argv[2]
offset = int(sys.argv[3])
outfile = sys.argv[4]

oligoT = 8
priming_len = 16

tempname = "_" + sys.argv[1] + ".filter.temp"
tempfile = pysam.Samfile(tempname, "wb", template=infile)

fasta_seqs = {}
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    fasta_seqs[seq_record.id] = seq_record
    
total = 0
mapped = 0
fake = 0
unmapped = 0
low_quality = 0
AT_high = 0
for read in infile.fetch():
    total += 1
    if read.is_unmapped: 
        unmapped += 1
    elif read.mapq < 30: 
        low_quality += 1
    elif ((read.seq.count('A') + read.seq.count('T')) / float(len(read.seq)) >= 0.80): 
        AT_high += 1
    else:
        chrx = infile.getrname(read.tid)
        if read.is_reverse:
            pos = read.pos + offset
            seq = str(fasta_seqs[chrx].seq[pos + read.rlen: pos + read.rlen + priming_len])
        else:
            pos = read.pos - offset
            seq = str(fasta_seqs[chrx].seq[pos - priming_len: pos].reverse_complement())  
        #if (seq.count('A') + seq.count('G')) / float(oligoT) < 0.75:
        if seq[0:oligoT].count('A') < oligoT and (seq.count('A') / float(priming_len)) < 0.75:
            tempfile.write(read)
            mapped += 1
        else:
            fake += 1
            
pysam.sort(tempname, outfile)
pysam.index(outfile + ".bam")
os.remove(tempname)
 
print "Mapped " + str(mapped) + " out of " + str(total) + " reads,", round(mapped / float(total) * 100, 2), "%", "low quality:", round(low_quality / float(total) * 100, 2), "%", "AT_high:", round(AT_high / float(total) * 100, 2), "%", "fakes:", round(fake / float(total) * 100, 2), "%"


tempfile.close()
infile.close()


