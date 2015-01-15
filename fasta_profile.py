import sys
import numpy
import matplotlib.pyplot as plt
from Bio import SeqIO


genome_file = sys.argv[1]
fasta_file = open(sys.argv[2], "r")
length = int(sys.argv[3])
end = int(sys.argv[4])
out_file = sys.argv[5]

genome_seqs = {}
A_genome = 0
T_genome = 0
C_genome = 0
G_genome = 0
for seq_record in SeqIO.parse(genome_file, "fasta"):
    genome_seqs[seq_record.id] = seq_record
    A_genome += seq_record.seq.count("A")
    T_genome += seq_record.seq.count("T") 
    C_genome += seq_record.seq.count("C") 
    G_genome += seq_record.seq.count("G")
    
tot = float(A_genome + T_genome + C_genome + G_genome)
A_genome = A_genome / tot
T_genome = T_genome / tot
C_genome = C_genome / tot
G_genome = G_genome / tot

N = {
    'A': numpy.array([0] * (length-end)),
    'T': numpy.array([0] * (length-end)),
    'C': numpy.array([0] * (length-end)),
    'G': numpy.array([0] * (length-end)),
}
L = numpy.array([0.0] * (length-end))

fasta_file.seek(0)
for line in fasta_file:
    if line[0] == ">": continue
    i = 0
    for x in line.strip()[:-end]:
        if x in ('A', 'T', 'G', 'C'):
            N[x][i] += 1
        L[i] += 1
        i += 1

A = list((N['A'] / L / A_genome / 4 ).astype(numpy.float16))
T = list((N['T'] / L / T_genome / 4 ).astype(numpy.float16))
G = list((N['G'] / L / G_genome / 4 ).astype(numpy.float16))
C = list((N['C'] / L / C_genome / 4 ).astype(numpy.float16))

r = range(-length,-end)

#print "%s,%s,%s,%s,%s" % ("pos", "A", "T", "G", "C")
#for i, a, t, g, c in zip(r, A, T, G, C):
    #print "%d,%f,%f,%f,%f" % (i, a, t, g, c)
#print "CG content:", C_genome + G_genome, "%"

#plt.ylim((0,0.8))
plt.plot(r, A, label='A', color='red', linewidth=2.0)
plt.plot(r, T, label='T', color='orange', linewidth=2.0)
plt.plot(r, G, label='G', color='blue', linewidth=2.0)
plt.plot(r, C, label='C', color='gray', linewidth=2.0)
plt.legend(loc=2, ncol=1, fancybox=True, shadow=True)
plt.xlabel('Position relative to intron 3\'')
plt.ylabel('Percentage')
plt.savefig(out_file + ".png")
