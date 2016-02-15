import sys, re, sets
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
total_reads = sets.Set()
mapped_reads = sets.Set()
for read in infile.fetch():
    total += 1
    total_reads.add(read.qname)
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
        sense = ["+", "-"][read.is_reverse]     
        #if (seq.count('A') + seq.count('G')) / float(oligoT) < 0.75:
        if seq[0:oligoT].count('A') < oligoT and (seq.count('A') / float(priming_len)) < 0.75:
            tempfile.write(read)
            mapped += 1
            mapped_reads.add(read.qname)
        else:
            fake += 1
            #internals.write(read)
tempfile.close()
infile.close()


pysam.sort(tempname, outfile)
pysam.index(outfile + ".bam")
os.remove(tempname)
 
print"Total reads: " + str(len(total_reads)), "Mapped reads: " + str(len(mapped_reads)), "Total mappings: " + str(total), "Final mappings: " + str(mapped), "low quality mappings: " + str(low_quality), "AT_high mappings: " +  str(AT_high), "fakes alignment: " +  str(fake) 



