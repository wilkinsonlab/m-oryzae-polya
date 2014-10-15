
import sys
#file = open(sys.stdin, 'r')
all = {}
distal = {}
proximal = {}
curr = {}
for line in sys.stdin:
    #print line
    gene, exon, exonBaseMean, dispersion, stat, pvalue, padj, a, b, fold, genomicData, countData_a1, countData_a2, countData_a3, countData_b1, countData_b2, countData_b3  = line.strip().split(',')
    chrx, pos, sense, start, end = exon.split("@")
    pos = int(pos);
    if fold == "NA" : fold = 0
    fold = float(fold)
    if not distal.has_key(gene):
        distal[gene] = fold, pos
        curr[gene] = pos
    if sense == '-':
        if pos > curr[gene]:    
            distal[gene] = fold, pos
            curr[gene] = pos
    else:
        if pos < curr[gene]:    
            distal[gene] = fold, pos
            curr[gene] = pos
            
    if not proximal.has_key(gene):
        proximal[gene] = fold, pos
        curr[gene] = pos
    if sense == '-':
        if pos < curr[gene]:    
            proximal[gene] = fold, pos
            curr[gene] = pos
    else:
        if pos > curr[gene]:    
            proximal[gene] = fold, pos
            curr[gene] = pos
            
    if not all.has_key(gene):
        all[gene] = []        
    if not all.has_key(gene):
        all[gene] = []       
    all[gene].append((fold, pos))
    
    
todistal  = 0
toproximal = 0
count = 0.0
for gene, (fold_d, pos_d) in distal.items():
    #print gene, pos_d, fold_d
    if fold_d > 0:
        for (fold_p, pos_p) in all[gene]:
            if fold_p < 0:
                print gene, abs(pos_p - pos_d), fold_d
                break
        todistal += 1 
    else:
        for (fold_p, pos_p) in all[gene]:
            if fold_p > 0:
                print gene, abs(pos_p - pos_d), fold_d
                break
        toproximal += 1
    count += 1
#print   toproximal/count, todistal/count



