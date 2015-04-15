import numpy, pysam, sys, math
from Bio import SeqIO
from Bio.Seq import Seq

fasta = {}
for sequence in SeqIO.parse(sys.argv[1], 'fasta'):
   fasta[str(sequence.id)] = str(sequence.seq)

samfile = pysam.Samfile( sys.argv[2],'r' )
i = 0
for read in samfile.fetch():
  i += 1
  if i % 1000000 == 0: sys.stderr.write(str(i) + '\n')
  name, val = read.qname.split('-')
  print read.qname  + '\t' + fasta[read.qname] + '\t' + str(read.flag) + '\t' + samfile.getrname(read.tid) + '\t' + str(read.pos+1) + '\t' + str(read.reference_end)  + '\t' + str(int(val) / float(read.get_tag('NH:I'))) + '\t' + read.cigarstring + '\t' + str(read.tags)
samfile.close()
