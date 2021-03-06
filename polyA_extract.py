
import sys
import numpy as np
import scipy.special


gff_file = open(sys.argv[1], 'r')
polyA_file = open(sys.argv[2], 'r')
expr_file = open(sys.argv[3], "r")
gap = int(sys.argv[4])
p_val_sites = float(sys.argv[5])
min_expr = int(sys.argv[6])
opt = sys.argv[7]

table = {}
polyA = {}
expr = {}
lines = {}

for line in gff_file:
    if line[0] == '#' or line[0] == '\n':
        continue
    items = line.split('\t')
    sense = items[6]
    # version 18
    #if items[2] == "gene":
    # version 21
    if items[2] in ("gene", "protein_coding_gene", "pseudogene", "pseudogenic_tRNA", "rRNA_gene", "RNA", "snoRNA_gene", "snRNA_gene", "tRNA_gene"):
        # version 21
        if items[8].find("Parent") != -1: continue
        chrx = items[0]
        start = int(items[3])
        end = int(items[4])
        sense = items[6]
        for x in items[8].split(';'):
            if x.split('=')[0] == "ID":
                transcript = x.split('=')[1].strip()
        table[transcript] = {}
        polyA[transcript] = []
        expr[transcript] = 1

for line in polyA_file:
    items = line.strip().split(" ")
    pos = int(items[1])
    chrx = items[2]
    sense = items[3]
    value = int(items[0])
    transcript = items[4]
    start = int(items[5])
    end = int(items[6])
    table[transcript][pos] = value
    lines[transcript] = (chrx, sense, start, end)


for transcript, array in table.items():
    if not lines.has_key(transcript):
        table[transcript] = {}
        continue
    chrx, sense, start, end = lines[transcript]
    arr = np.array(array.values())
    mean = np.mean(arr)
    std = np.std(arr)
    if std == 0.0:
        #std = 0.1
        table[transcript] = {}
        continue
    for pos, val in array.items():
	if scipy.special.ndtr(-((val - mean) / std)) >= p_val_sites or val < min_expr:
	#if float(val) / np.sum(arr) < 0.1 or val < 5:
	        del table[transcript][pos]

for transcript, array in table.items():
    if array == {}:
        continue
    keys = sorted(list(array.keys()))
    curr_k = keys[0]
    polyA[transcript].append([[curr_k, ], array[curr_k]])
    for iter_k in keys[1:]:
        if abs(iter_k - curr_k) <= gap:
            polyA[transcript][-1][0].append(iter_k)
            polyA[transcript][-1][1] += array[iter_k]
        else:
            polyA[transcript].append([[iter_k, ], array[iter_k]])
        curr_k = iter_k


tmp = []
for line in expr_file:
    transcript, value = line.strip().split("\t")
    tmp.append(int(value))
    expr[transcript] = int(value)
expr_all = np.array(tmp)

expr_std = np.std(expr_all)
expr_mean = np.mean(expr_all)


for transcript, polis in polyA.items():
    # skip low confident genes
    # if scipy.special.ndtr(-((expr[transcript] - expr_mean) / expr_std)) > p_val_genes:
    #    continue
    # Extract the confident polyA sites (sigle or APA), taking the highest
    for poli in polis:
        max_curr = poli[0][0]
        for pos in poli[0][1:]:
            if table[transcript][pos] > table[transcript][max_curr]:
                max_curr = pos
        # change  == 1 or > 1
        if opt == 'sgl' and len(polis) == 1:
                print table[transcript][max_curr], max_curr, lines[transcript][0], lines[transcript][1], transcript, lines[transcript][2], lines[transcript][3]
        elif opt == 'apa' and len(polis) > 1:
                print table[transcript][max_curr], max_curr, lines[transcript][0], lines[transcript][1], transcript, lines[transcript][2], lines[transcript][3]
        elif opt == 'all':
                print table[transcript][max_curr], max_curr, lines[transcript][0], lines[transcript][1], transcript, lines[transcript][2], lines[transcript][3]



