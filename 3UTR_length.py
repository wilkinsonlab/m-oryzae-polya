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
    #print round(x / count * 100, 5)


gff_file.close()
polyA_file.close()
