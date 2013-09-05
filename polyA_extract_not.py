import sys
import numpy as np
import scipy.special


notpolyA_file = open(sys.argv[1], 'r')
gap = int(sys.argv[2])
threshold = int(sys.argv[3])

table_sense = {}
table_antisense = {}
polyA_sense = {}
polyA_antisense = {}

for line in notpolyA_file:
    items = line.strip().split(" ")
    value = int(items[0])
    if value < threshold:
        continue
    pos = int(items[1])
    chrx = items[2]
    sense = items[3]
    if not table_sense.has_key(chrx):
        table_sense[chrx] = {}
    if not table_antisense.has_key(chrx):
        table_antisense[chrx] = {}
    if not polyA_sense.has_key(chrx):
        polyA_sense[chrx] = []
    if not polyA_antisense.has_key(chrx):
        polyA_antisense[chrx] = []
    if sense == "-":
        table_sense[chrx][pos] = value
    else:
        table_antisense[chrx][pos] = value

for chrx, array in table_sense.items():
    if array == {}:
        continue
    keys = sorted(list(array.keys()))
    curr_k = keys[0]
    polyA_sense[chrx].append([curr_k, ])
    for iter_k in keys[1:]:
        if abs(iter_k - curr_k) <= gap:
            polyA_sense[chrx][-1].append(iter_k)
        else:
            polyA_sense[chrx].append([iter_k, ])
        curr_k = iter_k

for chrx, array in table_antisense.items():
    if array == {}:
        continue
    keys = sorted(list(array.keys()))
    curr_k = keys[0]
    polyA_antisense[chrx].append([curr_k, ])
    for iter_k in keys[1:]:
        if abs(iter_k - curr_k) <= gap:
            polyA_antisense[chrx][-1].append(iter_k)
        else:
            polyA_antisense[chrx].append([iter_k, ])
        curr_k = iter_k


for chrx, polis in polyA_sense.items():
    for poli in polis:
        max_curr = poli[0]
        for pos in poli[1:]:
            if table_sense[chrx][pos] > table_sense[chrx][max_curr]:
                max_curr = pos
        print table_sense[chrx][max_curr], max_curr, chrx, "-"

for chrx, polis in polyA_antisense.items():
    for poli in polis:
        max_curr = poli[0]
        for pos in poli[1:]:
            if table_antisense[chrx][pos] > table_antisense[chrx][max_curr]:
                max_curr = pos
        print table_antisense[chrx][max_curr], max_curr, chrx, "+"
