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
    if feature in ('start_codon', 'stop_codon', 'exon'):
        for x in items[8].split(';'):
            if x.split('=')[0] == "Parent":
                name = x.split('=')[1].strip()
        name = re.sub(r'T.', "", name)
        if name not in table:
            table[name] = {}
            table[name]["exon"] = []
        sense = items[6]
        if sense == '+':
            if feature == 'start_codon':
                table[name][feature] = int(items[3])
            elif feature == 'stop_codon':
                table[name][feature] = int(items[4])
            elif feature == 'exon':
                table[name][feature].append((int(items[3]), int(items[4])))    
        elif sense == '-':
            if feature == 'start_codon':
                table[name][feature] = int(items[4])
            elif feature == 'stop_codon':
                table[name][feature] = int(items[3])
            elif feature == 'exon':
                table[name][feature].append((int(items[3]), int(items[4])))    

location = {key: {"five_utr": 0, "exon": 0, "intron": 0, "three_utr": 0} for key in table.keys()}
not_annotated = 0

for line in polyA_file:
    items = line.strip().split(" ")
    pos = int(items[1])
    chrx = items[2]
    value = int(items[0])
    transcript = items[4]
    sense = items[3]
    name = chrx + ":" + str(pos) + ":" + sense 
    if transcript not in table or "stop_codon" not in table[transcript] or "start_codon" not in table[transcript]:
            not_annotated += 1
            #print name +  " \t" + "NA"+ "\t" + transcript
            continue
    
    if sense == '-':
        if pos >= table[transcript]["stop_codon"]:
            location[transcript]["three_utr"] += 1
            #print  name + " \t" + "3'UTR" + "\t" + transcript
        elif pos <= table[transcript]["start_codon"]:
            location[transcript]["five_utr"] += 1
            #print name +  " \t" + "5'UTR" + "\t" + transcript
        elif pos > table[transcript]["start_codon"] and pos < table[transcript]["stop_codon"]:
            flag = False
            for exon in table[transcript]["exon"]:
                if pos >= exon[0] and pos <= exon[1]:
                    location[transcript]["exon"] += 1
                    #print name +  " \t" + "exon" + "\t" + transcript
                    flag = True
                    break
            if not flag:
                location[transcript]["intron"] += 1
                #print name +  " \t" + "intron" + "\t" + transcript
    elif sense == '+':
        if pos <= table[transcript]["stop_codon"]:
            location[transcript]["three_utr"] += 1
           # print name +  " \t" + "3'UTR" + "\t" + transcript
        elif pos >= table[transcript]["start_codon"]:
            #print name +  " \t" + "5'UTR" + "\t" + transcript
            location[transcript]["five_utr"] += 1
        elif pos < table[transcript]["start_codon"] and pos > table[transcript]["stop_codon"]:
            flag = False
            for exon in table[transcript]["exon"]:
                if pos >= exon[0] and pos <= exon[1]:
                    location[transcript]["exon"] += 1
                    #print name +  " \t" + "exon" + "\t" + transcript
                    flag = True
                    break
            if not flag:
                location[transcript]["intron"] += 1
                #print name +  " \t" + "intron" + "\t" + transcript

# where the polyA are located
three_utr = 0
five_utr = 0
exon = 0
intron = 0
# how many genes have each kind of APA
apa_three_utr = 0
apa_five_utr = 0
apa_cds = 0
count_poly = 0.0
count_apa = 0.0
for transcript, loc in location.items():
    if loc["intron"] > 0 or loc["exon"] > 0 :
        pass
        #print transcript
    three_utr += loc["three_utr"]
    five_utr += loc["five_utr"]
    exon += loc["exon"]
    intron += loc["intron"]
    if loc["three_utr"] > 1:
        apa_three_utr += 1
    if loc["three_utr"] >= 1 and loc["five_utr"] >= 1:
        apa_five_utr += 1
    if loc["three_utr"] >= 1 and (loc["exon"] >= 1 or loc["intron"] >= 1):
        apa_cds += 1
    count_poly += sum(loc.values())
    if sum(loc.values()) > 1:
        count_apa += 1

sys.stdout.write("%d,%d,%d\n" % (apa_three_utr, apa_five_utr, apa_cds))
#ys.stdout.write("%d,%d,%d,%d,%d\n" % (three_utr, five_utr, exon, intron, not_annotated))


gff_file.close()
polyA_file.close()
