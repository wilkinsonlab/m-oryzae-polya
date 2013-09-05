import sys
import numpy
f = open(sys.argv[1], 'r')
g = open(sys.argv[2], 'r')
dist = dict([(i, 0) for i in range(1, 10000)])
for line in f:
    items = line.strip().split(" ")
    chrx = items[2]
    pos = int(items[1])
    sense = items[3]
    for line_ in g:
        items_ = line_.strip().split(" ")
        chrx_ = items_[0]
        pos_ = int(items_[1])
        sense_ = items_[2]
        if chrx == chrx_:
            if sense_ == "+" and sense == "-":
                if pos > pos_:
                    if pos - pos_ >= 10000:
                        continue
                    dist[pos - pos_] += 1
                    break
            if sense_ == "-" and sense == "+":
                if pos < pos_:
                    if pos_ - pos >= 10000:
                        continue
                    dist[pos_ - pos] += 1
                    break
    g.seek(0)

for k, v in dist.items():
    print k, ",", v
