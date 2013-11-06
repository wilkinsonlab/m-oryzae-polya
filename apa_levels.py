import sys
import re
import pysam
import matplotlib

count_file = open(sys.argv[1], "r")
expr_file = open(sys.argv[2], "r")


expr = [{}, {}, {}]
polyA_table = {}
for line in expr_file:
    gene, val1, val2, val3 = line.strip().split("\t")
    expr[0][gene] = int(val1)
    expr[1][gene] = int(val2)
    expr[2][gene] = int(val3)
    polyA_table[gene] = []

for line in count_file:    
    polyA, val1, val2, val3 = line.strip().split("\t")
    items = polyA.split(":")
    gene = items[3]
    sense = items[2]
    pos = int(items[1])
    if sense == '-':    
        polyA_table[gene].append((pos, ( int(val1) / (float(expr[0][gene]) + 0.001) + int(val2) / (float(expr[1][gene]) + 0.001) + int(val3) / (float(expr[2][gene]) + 0.001) ) / 3))
    elif sense == '+':
        polyA_table[gene].insert(0, (pos, ( int(val1) / (float(expr[0][gene]) + 0.001) + int(val2) / (float(expr[1][gene]) + 0.001) + int(val3) / (float(expr[2][gene]) + 0.001) ) / 3))        

levels = [[] for x in range(6)]
proximal = [0 for x in range(100)]
distal = [0 for x in range(100)]
for gene, polyAs in polyA_table.items():
    if len(polyAs) < 2: 
        continue
    for i, (pos, val) in enumerate(polyAs):
        levels[i].append(val)
    proximal[int(polyAs[0][1]*100)] += 1
    distal[int(polyAs[-1][1]*100)] += 1
    
for x in levels:
    if len(x) > 0:
        print sum(x) / len(x)   
print        

for x in proximal: print x / float(sum(proximal))
print
for x in distal: print x / float(sum(distal))