import sys
import re
import pysam
import matplotlib

gff_file = open(sys.argv[1], 'r')
bam_file = pysam.Samfile(sys.argv[2], 'rb')
feature = sys.argv[3] # start, stop, init, term, 3utr

starts = {}
stops = {}
inits = {}
terms = {}
utr = {}
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
        stops[name] = (pos, sense)
    elif items[2] == "start_codon":
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
        starts[name] = (pos, sense)
    elif items[2] == "gene":
        chrx = items[0]
        for x in items[8].split(';'):
            if x.split('=')[0] == "ID":
                name = x.split('=')[1].strip()
        sense = items[6]
        if sense == '+':
            inits[name] = (int(items[3]), sense)
            terms[name] = (int(items[4]), sense)
        else:
            inits[name] = (int(items[4]), sense)
            terms[name] = (int(items[3]), sense)
    elif items[2] == "three_prime_UTR":
        chrx = items[0]
        for x in items[8].split(';'):
            if x.split('=')[0] == "Parent":
                name = x.split('=')[1].strip()
        sense = items[6]
        if sense == '+':
            utr[name] = (int(items[3]), sense)
        else:
            utr[name] = (int(items[4]), sense)
    

         
if feature == "start":
    select = starts        
elif feature == "stop":
    select = stops
elif feature == "init":
    select = inits
elif feature == "term":
    select = terms
elif feature == "3utr":
    select = utr   

offset = 7
distance = [0] * 500
dists = 0
count = 0.0
already = {}
for read_copy in bam_file.fetch():
    read = read_copy
    if read.is_reverse:
        read.pos += offset
        pos = read.pos + read.rlen
    else:
        read.pos -= offset
        pos = read.pos
    chrx = bam_file.getrname(read.tid)
    if already.has_key((chrx, pos)):
        continue
    else:
        already[(chrx, pos)] = True    
    sense = ('-', '+')[read.is_reverse]
    
    dist = 10000000
    for gene, (loc, g_sense) in select.items():
        if sense == g_sense:
            if abs(pos - loc) < dist:
                dist = abs(pos - loc)
                

    if dist < 500 and dist >= 0:
        distance[dist] += 1
        dists += dist
        count += 1
        if count > 100000: 
            break

print "%.2f" % (dists / count)
for x in distance:
    print round(x / count * 100, 5)


gff_file.close()
bam_file.close()
