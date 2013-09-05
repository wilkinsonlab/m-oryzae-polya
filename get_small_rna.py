import pysam
import sys


samfile = pysam.Samfile(sys.argv[1], "r")
distance = int(sys.argv[2])
offset = int(sys.argv[3])

mates = {}
hits = {}
dists = {}
for read in samfile.fetch():
    if not read.is_unmapped and read.mapq >= 30 and ((read.seq.count('A') + read.seq.count('T')) / float(len(read.seq)) < 0.80):
        if not mates.has_key(read.qname):
            mates[read.qname] = (read.pos, len(read.seq))
        else:
            if read.pos < mates[read.qname][0]:
                dist = mates[read.qname][0] + mates[read.qname][1] - read.pos
            else:
                dist = read.pos + len(read.seq) - mates[read.qname][0]

            label = pysam.Samfile.getrname(samfile, read.tid) + " " + str((mates[read.qname][0] + mates[
                                                                           read.qname][1], mates[read.qname][0])[read.is_reverse]) + " " + ('-', '+')[read.is_reverse]
            if not hits.has_key(label):
                hits[label] = []
                dists[label] = []
            if dist < distance:
                hits[label].append(True)
                dists[label].append(dist)
                del mates[read.qname]
            else:
                hits[label].append(False)


for k, v in hits.items():
    chrx = k.split(" ")[0]
    sense = k.split(" ")[2]
    if sense == "+":
        pos = int(k.split(" ")[1]) - offset
    else:
        pos = int(k.split(" ")[1]) + offset
    print "%s\t%d\t%s\t%d\t%f\t%f\t%s" % (chrx, pos, sense, len(v), v.count(True) / float(len(v)) * 100, sum(dists[k]) / (len(dists[k])+0.1), dists[k])
