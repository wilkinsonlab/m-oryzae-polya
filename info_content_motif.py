import sys, re, numpy, math
from Bio import SeqIO

g = open(sys.argv[1], 'r')
f = open(sys.argv[2], 'r')
m = sys.argv[3]
length = int(sys.argv[4])

m = m.replace('R', '[GA]').replace('Y', '[TC]').replace('S', '[GC]').replace('W', '[TA]').replace('K', '[GT]').replace('M', '[AC]').replace('D', '[GTA]').replace('H', '[TAC]').replace('B', '[GTC]').replace('V', '[GAC]').replace('N', '[ATGC]')

A_genome = 0
T_genome = 0
C_genome = 0
G_genome = 0
N_genome = 0

for seq_record in SeqIO.parse(g, 'fasta'):
    A_genome += str(seq_record.seq).upper().count('A')
    T_genome += str(seq_record.seq).upper().count('T') 
    C_genome += str(seq_record.seq).upper().count('C') 
    G_genome += str(seq_record.seq).upper().count('G')

tot = float(A_genome + T_genome + C_genome + G_genome)
A_genome = A_genome / tot
T_genome = T_genome / tot
C_genome = C_genome / tot
G_genome = G_genome / tot
genome = (A_genome, C_genome, G_genome, T_genome)

lines = 0.0
count = []
for x in range(length): count.append([1,1,1,1])
count = numpy.array(count)          
                     
for line in f:
  l = line.strip().upper()
  find = re.findall(m, l)
  if len(find) > 0:
    for i, x in enumerate(find[-1]):
      if x == 'A':
        count[i][0] += 1
      elif x == 'C':
        count[i][1] += 1 
      elif x == 'G': 
        count[i][2] += 1
      elif x == 'T':
        count[i][3] += 1   
  lines += 1 
  
count = count / lines
count /=  count.sum(axis=1)[:,numpy.newaxis]

info = 0
for i, l in enumerate(count):
  for k, n in enumerate(l):
    info += n * math.log(n / genome[k], 2)
    
print info