
import sys
file = open(sys.argv[1], 'r')
prox1 = {}
dist1 = {}
curr1 = {}
prox2 = {}
dist2 = {}
for line in file:
    polyA, a1, a2, a3, b1, b2, b3 = line.strip().split('\t')
    gene, data = polyA.split(':')
    chrx, pos, sense, start, end = data.split("@")
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

nochange = 0
todistal  = 0
toproximal = 0
count = 0.0
for k,v in prox1.items():
    if (prox1[k]-dist1[k] > 0 and prox2[k]-dist2[k] < 0):
        todistal += 1
    elif (prox1[k]-dist1[k] < 0 and prox2[k]-dist2[k] > 0):
        toproximal += 1
        #print k
    else:
        nochange += 1
    count += 1    
print   nochange/count,  toproximal/count, todistal/count            
