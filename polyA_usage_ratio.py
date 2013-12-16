
import sys
file = open(sys.argv[1], 'r')
distal = {}
curr = {}
for line in file:
    gene, exon, disp, pvalue, padjust, meanBase, fold = line.strip().split(',')
    chrx, pos, sense, start, end = exon.split("@")
    pos = int(pos);
    fold = float(fold)
    if not distal.has_key(gene):
        distal[gene] = fold
        curr[gene] = pos
    if sense == '-':
        if pos > curr[gene]:    
            distal[gene] = fold
            curr[gene] = pos
    else:
        if pos < curr[gene]:    
            distal[gene] = fold
            curr[gene] = pos
    
                                 
todistal  = 0
toproximal = 0
count = 0.0
for gene, fold in distal.items():
    if fold > 0:
        todistal += 1 
    else:
        toproximal += 1
    count += 1
        
print   toproximal/count, todistal/count            
