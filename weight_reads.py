import numpy, pysam, sys, math
from Bio import SeqIO
from Bio.Seq import Seq

fasta = {}
for sequence in SeqIO.parse(sys.argv[1], 'fasta'):
   fasta[str(sequence.id)] = str(sequence.seq)

counts = {}
samfile = pysam.Samfile(sys.argv[2],'r' )
i=0
for read in samfile.fetch():
  break
  i+=1
  if i % 1000000 == 0: sys.stderr.write(str(i) + '\n')
  if read.flag == 16 or read.flag == 272:
    seq = str(Seq(read.query_sequence).reverse_complement())
  elif read.flag == 0 or read.flag == 256:
    seq = read.query_sequence
  name = read.qname + '_' + seq
  if not counts.has_key(name): counts[name] = 0.0
  counts[name] += 1
samfile.close()

samfile = pysam.Samfile( sys.argv[2],'r' )
i = 0
for read in samfile.fetch():
  i += 1
  if i % 1000000 == 0: sys.stderr.write(str(i) + '\n')
  name, val = read.qname.split('-')
  print read.qname  + '\t' + fasta[read.qname] + '\t' + str(read.flag) + '\t' + samfile.getrname(read.tid) + '\t' + str(read.pos+1) + '\t' + str(read.reference_end)  + '\t' + str(int(val) / float(read.get_tag('NH:I'))) 
samfile.close()
