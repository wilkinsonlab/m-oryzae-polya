import sys
import numpy as np
from Bio import SeqIO
from operator import itemgetter

def median(l):
    half = len(l) / 2
    l.sort()
    if len(l) % 2 == 0:
        return (l[half-1] + l[half]) / 2.0
    else:
        return l[half]

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() 
    med = np.median(arr)
    return np.median(np.abs(arr - med))

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
			print line
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
            intron.seq = genome_seqs[intron.chrx].seq[start-1-3:end]
        elif intron.sense == "-":
            intron.start = end
            intron.end = start 
            intron.seq = genome_seqs[intron.chrx].seq[start-1:end+3].reverse_complement()
        introns.append(intron)
        genes[gene_id].introns.append(intron)
        
count = 0
num = []
length = []
donors = {}
acceptors = {}
for gene_id, gene in genes.items():
    #if len(gene.exons) > 1:
    if genes[gene_id].introns != []: 
        count += 1            
        num.append(len(gene.exons)-1)
for intron in introns:           
    length.append(len(intron.seq) ) 
    donor = str(intron.seq[0:9])
    if donors.has_key(donor):
        donors[donor] += 1
    else:
        donors[donor] = 1   
    acceptor = str(intron.seq[-12:])
    #print acceptor
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
## # # #        
#print "number of introns:\t" + str(len(introns))
#if len(num) != 0:
	#print "protein coding genes:\t", len(genes.keys())   
	#print "protein coding genes containg introns:\t%d" % (count)
	#print "number of introns:\t" + str(len(introns))
	#print "average number of introns per gene:\t%.1f" % (sum(num) / float(len(num)))
	#print "median intron length:\t%d" % (median(length))
	#print "average intron length:\t%d" % (sum(length) / len(length)) 
	#print "median absolute deviation:\t" +  str(mad(length)) 
	#print "standard deviation:\t" +  str(np.std(length))
#else:
	#print "protein coding genes:\t", len(genes.keys())   
	#print "protein coding genes containg introns:\t%d" % (count)
	#print "number of introns:\t" + str(len(introns))
	#print "average number of introns per gene:\t%.1f" % 0
	#print "average intron length:\t%d" % 0
	#print "median intron length:\t%d" % 0
	#print "median absolute deviation:\t%.2f" %  0
	#print "standard deviation:\t%.2f" %  0

#for donor, value in donors.items():
#      print donor + "\t" + str(value / float(len(introns)) )
#for acceptor, value in acceptors.items():
##          print acceptor + "\t" + str(value / float(len(introns)) )
i=0
for intron in introns:
	##if len(intron.seq) > 250: continue
	print ">s_" + str(i)
	print "N" * (50 - len(intron.seq)) + intron.seq[-50:-3]
        i+=1
	if i > 10000: break
         

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
# print "ratio:", sum(ratios) / float(len(ratios))

#count = .0
#for intron in introns:
	  #if str(intron.seq[3:5]) == "AT" and str(intron.seq[-2:]) == "AC":
		 #count += 1 
#print count / len(introns) 
