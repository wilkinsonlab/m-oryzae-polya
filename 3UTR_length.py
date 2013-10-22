import sys
import re
import pysam
import matplotlib

gff_file = open(sys.argv[1], 'r')
polyA_file = open(sys.argv[2], 'r')

stops = {}
for line in gff_file:
    if line[0] == '#' or line[0] == '\n':
        continue
    items = line.split('\t')
    if items[2] == "stop_codon":
        chrx = items[0]
        for x in items[8].split(';'):
            if x.split('=')[0] == "Parent":
                name = x.split('=')[1].strip()
        name = re.sub(r'T.', "", name)
        sense = items[6]
        if sense == '+':
            pos = int(items[4])
        else:
            pos = int(items[3])
        stops[name] = pos

distance = [0] * 500
dists = 0
count = 0.0
for line in polyA_file:
    items = line.strip().split(" ")
    pos = int(items[1])
    sense = items[3]
    transcript = items[4]
    value = int(items[0])
    if not stops.has_key(transcript):
        continue
    if sense == '-' and pos > stops[transcript]:
        dist = pos - stops[transcript]
    elif sense == '+' and pos < stops[transcript]:
        dist = stops[transcript] - pos
    else:
        continue
    if dist < 500:
        distance[dist] += 1
    dists += dist
    count += 1

print "%.2f" % (dists / count)
for x in distance:
    pass
    print round(x / count * 100, 5)


gff_file.close()
polyA_file.close()


python -c "

file = open('_n', 'r')
prox1 = {}
dist1 = {}
curr1 = {}
prox2 = {}
dist2 = {}
for line in file:
    polyA, a1, a2, a3, b1, b2, b3 = line.strip().split('\t')
    chrx, pos, sense, gene, start, end = polyA.split(':')
    a1 = int(a1);a2 = int(a2);a3 = int(a3);b1 = int(b1);b2 = int(b2);b3 = int(b3);
    pos = int(pos);
    if not prox1.has_key(gene):
        prox1[gene] = a1 + a2 + a3 
        prox2[gene] = b1 + b2 + b3 
        dist1[gene] = a1 + a2 + a3 
        dist2[gene] = b1 + b2 + b3 
        curr1[gene] = pos
    else:
        prox1[gene] += a1 + a2 + a3 
        prox2[gene] += b1 + b2 + b3 
        if sense == '-':
            if pos > curr1[gene]:    
                dist1[gene] = a1 + a2 + a3
                dist2[gene] = b1 + b2 + b3
                curr1[gene] = pos
        else:
            if pos < curr1[gene]:    
                dist1[gene] = a1 + a2 + a3
                dist2[gene] = b1 + b2 + b3
                curr1[gene] = pos
    
                                 
for gene,val in dist1.items():
    prox1[gene] -= val
for gene,val in dist2.items():
    prox2[gene] -= val

for k,v in prox1.items():

        print prox1[k]/(dist1[k]+0.1) ,  prox2[k]/(dist2[k]+0.1)
     
"

