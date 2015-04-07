import numpy, pysam, sys, math
from Bio import SeqIO

fasta = {}
for seq in SeqIO.parse(sys.argv[1], 'fasta'):
   fasta[str(seq.id)] = str(seq.seq)

counts = {}
samfile = pysam.Samfile(sys.argv[2],'r' )
i=0
for read in samfile.fetch():
  i+=1
  if i % 1000000 == 0: sys.stderr.write(str(i) + '\n')
  name = read.qname + '_' + str(read.flag)
  if not counts.has_key(name): counts[name] = 0.0
  counts[name] += 1
samfile.close()

samfile = pysam.Samfile( sys.argv[2],'r' )
i = 0
for read in samfile.fetch():
  i += 1
  if i % 1000000 == 0: sys.stderr.write(str(i) + '\n')
  name, val = read.qname.split('-')
  name += '-' + val + '_' + str(read.flag)
  print name  + '\t' + fasta[name.split('_')[0]] + '\t' + samfile.getrname(read.tid) + '\t' + str(read.pos+1) + '\t' + str(read.reference_end)  + '\t' + str(int(val) / counts[name]) 
samfile.close()

