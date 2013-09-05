import sys
count_expr = open(sys.argv[1], 'r')
count_poly = open(sys.argv[2], 'r')

norm = {}
for line in count_expr:
    (gene, a1, a2, a3, b1, b2, b3) = line.strip().split("\t")
    norm[gene] = (float(a1) / (float(b1) + 0.001), float(
        a2) / (float(b2) + 0.001), float(a3) / (float(b3) + 0.001))
for line in count_poly:
    (data, a1, a2, a3, b1, b2, b3) = line.strip().split("\t")
    gene = data.split(':')[3]
    print data + "\t" + a1 + "\t" + a2 + "\t" + a3 + "\t" + str(int(int(b1) * norm[gene][0])) + "\t" + str(int(int(b2) * norm[gene][1])) + "\t" + str(int(int(b3) * norm[gene][2]))
