import sys
import numpy as np
import scipy.special
import sys
import numpy
import scipy.stats as stat


import numpy.ma as ma

def MAD(a, c=0.6745, axis=None):


    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)

    return m



gff_file = open(sys.argv[1], 'r')
polyA_file = open(sys.argv[2], 'r')
expr_file = open(sys.argv[3], "r")
gap = int(sys.argv[4])
z_score = float(sys.argv[5])
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

norms = {}
linee = expr_file.readlines()
ref = float(linee[0].strip().split('\t')[1])
for line in linee:
  gene, val = line.strip().split('\t')
  val = int(val)
  if val == 0:
    new = 0
  else:
    new = ref / val

  norms[gene] = float(new)

vals = []
for line in polyA_file:
    items = line.strip().split(" ")
    pos = int(items[1])
    chrx = items[2]
    sense = items[3]
    value = int(items[0])
    transcript = items[4]
    start = int(items[5])
    end = int(items[6])
    if int(value) <= min_expr: continue
    val_norm = int(value) * norms[transcript]
    table[transcript][pos] = value
    lines[transcript] = (chrx, sense, start, end)
    vals.append(val_norm)

arr = numpy.array(vals)
mean = numpy.mean(arr) 
median = numpy.median(arr)
std = numpy.std(arr) 
mad = MAD(arr)
#print mean, median, mad
for transcript, array in table.items():
    if not lines.has_key(transcript):
        table[transcript] = {}
        continue
    chrx, sense, start, end = lines[transcript]
    for pos, val in array.items():
        val_norm = int(val) * norms[transcript]
	#print "DEBUG", table[transcript][pos], pos, lines[transcript][0], lines[transcript][1], transcript, (val_norm - median) / mad
	if (val_norm - median) / mad < z_score:
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

