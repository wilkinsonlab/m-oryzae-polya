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
    if items[2] == "three_prime_UTR":
        chrx = items[0]
        for x in items[8].split(';'):
            if x.split('=')[0] == "Parent":
                name = x.split('=')[1].strip()
        name = re.sub(r'T.', "", name)
        sense = items[6]
        if sense == '+':
            pos = int(items[3])
        else:
            pos = int(items[4])
        stops[name] = pos

distance = [0] * 500
count = 0.0
dists = 0
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


for x in distance:
    print round(x / count * 100, 5)    

print "%.2f" % (dists / count)


gff_file.close()
polyA_file.close()



python -c "
import sys

class Thing:
    name=''
    seq=''
    struct=''
    
file='_res_WT-CM_100'
l = []
for line in open(file, 'r'):
    if line[0] == '>':
        t = Thing()
        t.name = line
    if line[0] in ['A', 'U', 'G', 'C']:
        t.seq = line
    if line[0] in ['.', '(', ')']:
        t.struct = line    
        l.append(t)
for x in l:
    print x.name.strip()
    for i, e in enumerate(x.struct):
        if e == '.':
            sys.stdout.write(x.seq[i])
        elif e == ' ':
            break    
    print            

"





