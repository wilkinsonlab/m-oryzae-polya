import sys
import re
import pysam

#
# Given the bad annotation, I only consider start and stop codons
# Result specificy APA only when there a polyA at least in the 3'UTR
#


gff_file = open(sys.argv[1], 'r')
polyA_file = open(sys.argv[2], 'r')
# opt = sys.argv[3]

table = {}
for line in gff_file:
    if line[0] == '#' or line[0] == '\n':
        continue
    items = line.split('\t')
    feature = items[2]
    if feature in ('start_codon', 'stop_codon'):
        for x in items[8].split(';'):
            if x.split('=')[0] == "Parent":
                name = x.split('=')[1].strip()
        name = re.sub(r'T.', "", name)
        if name not in table:
            table[name] = {}
        sense = items[6]
        if sense == '+':
            if feature == 'start_codon':
                table[name][feature] = int(items[3])
            elif feature == 'stop_codon':
                table[name][feature] = int(items[4])
        elif sense == '-':
            if feature == 'start_codon':
                table[name][feature] = int(items[4])
            elif feature == 'stop_codon':
                table[name][feature] = int(items[3])

location = {key: [0, 0, 0] for key in table.keys()}
for line in polyA_file:
    items = line.strip().split(" ")
    pos = int(items[1])
    chrx = items[2]
    value = int(items[0])
    transcript = items[4]
    sense = items[3]
    if transcript not in table or "stop_codon" not in table[transcript] or "start_codon" not in table[transcript]:
            continue
    if sense == '-':
        if pos > table[transcript]["stop_codon"]:
            location[transcript][2] += 1
        elif pos < table[transcript]["start_codon"]:
            location[transcript][0] += 1
        elif pos > table[transcript]["start_codon"] and pos < table[transcript]["stop_codon"]:
            location[transcript][1] += 1
    elif sense == '+':
        if pos < table[transcript]["stop_codon"]:
            location[transcript][2] += 1
        elif pos > table[transcript]["start_codon"]:
            location[transcript][0] += 1
        elif pos < table[transcript]["start_codon"] and pos > table[transcript]["stop_codon"]:
            location[transcript][1] += 1

three_utr = 0
five_utr = 0
cds = 0
count = 0
for transcript, loc in location.items():
    if loc[1] == 1:
        print transcript
    three_utr += loc[2]
    five_utr += loc[0]
    cds += loc[1]
    count += sum(loc)


#sys.stdout.write("%.2f%%,%.2f%%,%.2f%%\n" % (three_utr / float(count) * 100, five_utr / float(count) * 100, cds / float(count) * 100))


gff_file.close()
polyA_file.close()
