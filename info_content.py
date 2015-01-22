import sys, math
from Bio import SeqIO

genome_file = sys.argv[1]
seq_file = open(sys.argv[2], "r")
length = len(seq_file.readline())

seq_file.seek(0)
fasta_seqs = {}
A_genome = 0
T_genome = 0
C_genome = 0
G_genome = 0
N_genome = 0

for seq_record in SeqIO.parse(genome_file, "fasta"):
    A_genome += str(seq_record.seq).upper().count("A")
    T_genome += str(seq_record.seq).upper().count("T") 
    C_genome += str(seq_record.seq).upper().count("C") 
    G_genome += str(seq_record.seq).upper().count("G")

tot = float(A_genome + T_genome + C_genome + G_genome)
A_genome = A_genome / tot
T_genome = T_genome / tot
C_genome = C_genome / tot
G_genome = G_genome / tot

count = [{'A':1, 'T':1,'G':1,'C':1} for i in range(length)]
lines = float(0)
for line in seq_file:
    if 'N' in line: continue 
    for i, x in enumerate(line.strip()):
        #if x == 'N': continue 
        count[i][x] += 1
    lines += 1

info_content = 0    
for l in range(length):
    info_content += count[l]['A'] / lines * math.log(count[l]['A'] / lines / A_genome, 2)
    info_content += count[l]['T'] / lines * math.log(count[l]['T'] / lines / T_genome, 2)
    info_content += count[l]['G'] / lines * math.log(count[l]['G'] / lines / G_genome, 2)
    info_content += count[l]['C'] / lines * math.log(count[l]['C'] / lines / C_genome, 2) 

print info_content        