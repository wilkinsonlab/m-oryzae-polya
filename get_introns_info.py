import sys
from Bio import SeqIO
from operator import itemgetter

def median(l):
    half = len(l) / 2
    l.sort()
    if len(l) % 2 == 0:
        return (l[half-1] + l[half]) / 2.0
    else:
        return l[half]

class Gene:
    def __init__(self):
        self.id = ""
        self.chrx = ""
        self.sense = ""
        self.exons = []
        self.introns = []
        
class Intron:
    def __init__(self):
        self.gene_id = ""
        self.chrx = ""
        self.sense = ""
        self.start = 0
        self.end = 0
        self.seq = ""


genome_file =open(sys.argv[1], 'r') 
gff_file = open(sys.argv[2], 'r')
spacer = sys.argv[3]
ID = sys.argv[4]

genome_seqs = {}
for seq_record in SeqIO.parse(genome_file, "fasta"):
    genome_seqs[seq_record.id] = seq_record.upper()

genes = {}
for line in gff_file:
    if line[0] == '#' or line[0] == '\n':
        continue
    items = line.split('\t')
    if items[2] == "exon" :
        for x in items[8].split(';'):
            if x.strip().split(spacer)[0] == ID:
                gene_id = x.strip().split(spacer)[1].replace("\"", "")
        if not genome_seqs.has_key(items[0]):
            continue      
        if genes.has_key(gene_id):
            genes[gene_id].exons.append((int(items[3]), int(items[4])))
        else:     
            gene = Gene()
            gene.id = gene_id
            gene.chrx = items[0]
            gene.sense = items[6]
            gene.exons.append((int(items[3]), int(items[4])))
            genes[gene_id] = gene

for gene in  genes.values():
    gene.exons = sorted(gene.exons)

introns = []
for gene_id, gene in genes.items():
    for i, (s, e) in enumerate(gene.exons):
        if i == len(gene.exons)-1: break
        start, end = e, gene.exons[i+1][0]
        start += 1
        end -= 1
        if abs(start-end)<10: continue
        intron = Intron()
        intron.gene_id = gene_id
        intron.sense = gene.sense
        intron.chrx = gene.chrx
        if intron.sense == "+":
            intron.start = start
            intron.end = end
            intron.seq = genome_seqs[intron.chrx].seq[start-4:end+1]
        elif intron.sense == "-":
            intron.start = end
            intron.end = start 
            intron.seq = genome_seqs[intron.chrx].seq[start-2:end+3].reverse_complement()
        introns.append(intron)
        genes[gene_id].introns.append(intron)
        
count = 0
num = []
length = []
donors = {}
acceptors = {}
for gene in genes.values():
    if len(gene.exons) > 1:
        count += 1            
        num.append(len(gene.exons)-1)
for intron in introns:           
    length.append(len(intron.seq) ) 
    donor = str(intron.seq[-4:-1])
    if donors.has_key(donor):
        donors[donor] += 1
    else:
        donors[donor] = 1   
    acceptor = str(intron.seq[3:9])
    if acceptors.has_key(acceptor):
        acceptors[acceptor] += 1
    else:
        acceptors[acceptor] = 1  


# i=0
# for intron in introns:
#         i+=1
#         if intron.start > intron.end :
#             t = intron.start
#             intron.start = intron.end
#             intron.end = t
# 
#         print intron.chrx + "\t" + "protein_coding" + "\t" + "intron" + "\t" + str(intron.start) + "\t" + str(intron.end) + "\t" + "." + "\t" + intron.sense + "\t" + "." + "\t" + "gene_id \"" + intron.gene_id + "\"; ID intron_" + str(i) + "_"
# # # #        
# print " \t" +  sys.argv[2]      
# print "protein coding genes:\t", len(genes.keys())   
# print "protein coding genes containg introns:\t%d" % (count)
# print "average number of introns per gene:\t%.1f" % (sum(num) / float(len(num)))
# print "average intron length:\t%d" % (median(length))  
# print "number of introns:" + str(len(introns))  
for donor, value in donors.items():
     print donor + "\t" + str(value / float(len(introns)) )
#for acceptor, value in acceptors.items():
#         print acceptor + "\t" + str(value / float(len(introns)) )
#for intron in introns:
        #print ">" + intron.gene_id
        #print intron.seq #"N" * (97 - len(intron.seq)) + intron.seq[-97:]
# ratios = []        
# for gene in genes.values():
#     if len(gene.introns) == 0: continue
#     exonic = 0
#     intronic = 0
#     for s, e in gene.exons:
#         exonic += abs(s-e)
#     for intron in gene.introns:
#         intronic += abs(intron.start-intron.end)
#     ratios.append(  round(intronic / float(intronic+exonic), 2))   
# print sum(ratios) / float(len(ratios))