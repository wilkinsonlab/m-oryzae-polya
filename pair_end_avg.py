import sys

infile = open(sys.argv[1], "r")
table = {}
num = 0
count = 0
dists = dict.fromkeys(range(1000), 0)

for line in infile:
    if line[0] == '@':
        continue
    items = line.strip().split("\t")
    name = items[0]
    pos = int(items[3])
    qual = int(items[4])
    seq = items[9]
    if qual < 40:
        continue
    if not table.has_key(name):
        table[name] = (pos, len(seq))
    else:
        if pos < table[name][0]:
            dist = table[name][0] - (pos + len(seq))
        else:
            dist = pos - (table[name][0] + table[name][1])
        if dist < 1000 and dist > 0:
            dists[dist] += 1
        del table[name]
        # num += dist
        # count += 1

# print num / float(count)
for i, x in dists.items():
    print i, x

infile.close()
