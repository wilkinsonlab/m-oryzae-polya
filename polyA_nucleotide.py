import sys
import numpy
import matplotlib.pyplot as plt
from Bio import SeqIO


fasta_file = sys.argv[1]
polyA_file = open(sys.argv[2], "r")
start = int(sys.argv[3])
end = int(sys.argv[4])
opt = sys.argv[5]
out = open(sys.argv[6], 'w')

fasta_seqs = {}
A_genome = 0
T_genome = 0
C_genome = 0
G_genome = 0
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    fasta_seqs[seq_record.id] = seq_record
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
    'A': numpy.array([0] * (abs(start - end) + 1)),
    'T': numpy.array([0] * (abs(start - end) + 1)),
    'C': numpy.array([0] * (abs(start - end) + 1)),
    'G': numpy.array([0] * (abs(start - end) + 1)),
}
L = numpy.array([0.0] * (abs(start - end) + 1))

for line in polyA_file:
    items = line.strip().split(" ")
    pos = int(items[1])
    chrx = items[2]
    value = int(items[0])
    sense = items[3]
    if len(items) >= 5:
        gene = items[4]
    else:
        gene = ""    

    # if value < 10: continue

    if sense == '+':
        stream = fasta_seqs[chrx].seq[
            pos - end - 1: pos - start].reverse_complement()
    else:
        stream = fasta_seqs[chrx].seq[pos + start - 1: pos + end]

    stream = str(stream)

    i = 0
    for x in stream:
        if x in ('A', 'T', 'G', 'C'):
            N[x][i] += 1
        L[i] += 1
        i += 1

    if opt == "print":
        # for x in range(value):
        out.write(">" + chrx + ":" + str(pos) + ":" + str(value) + ":" + sense + ":" + gene + "\n")
        out.write(stream + "\n")

if opt == "view":
    A = list((N['A'] / L / A_genome / 4).astype(numpy.float16))
    T = list((N['T'] / L / T_genome / 4).astype(numpy.float16))
    G = list((N['G'] / L / G_genome / 4).astype(numpy.float16))
    C = list((N['C'] / L / C_genome / 4).astype(numpy.float16))
    r = range(start, end + 1)
    print "%s,%s,%s,%s,%s" % ("pos", "A", "T", "G", "C")
    for i, a, t, g, c in zip(r, A, T, G, C):
        print "%d,%f,%f,%f,%f" % (i, a, t, g, c)
    # print 'A: %.6f %%, T: %.6f %%, G: %.6f %%, C: %.6f %%, ' %
    # (sum(A)/len(A), sum(T)/len(T),sum(G)/len(G),sum(C)/len(C))
    print "CG content:", C_genome + G_genome, "%"

    plt.plot(r, A, label='A', color='red', marker='+')
    plt.plot(r, T, label='T', color='orange', marker='*')
    plt.plot(r, G, label='G', color='blue', marker='o')
    plt.plot(r, C, label='C', color='gray', marker='^')
    plt.legend(loc='upper center', bbox_to_anchor=(
        0.5, 1.05), ncol=3, fancybox=True, shadow=True)
    plt.legend()
    plt.xlabel('Position relative to polyA site')
    plt.ylabel('Percentage')
    plt.show()
