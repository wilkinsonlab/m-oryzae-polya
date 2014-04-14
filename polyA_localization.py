import sys
import re
import pysam

# REMBERER: ANNOTATION CAN BE REDUNDANT


gff_file = open(sys.argv[1], 'r')
polyA_file = open(sys.argv[2], 'r')
# opt = sys.argv[3]
features = ("gene", "protein_coding_gene", "pseudogene", "pseudogenic_tRNA", "rRNA_gene", "RNA", "snoRNA_gene", "snRNA_gene", "tRNA_gene", 'three_prime_UTR', 'CDS', 'five_prime_UTR')
table = {}
for line in gff_file:
    if line[0] == '#' or line[0] == '\n':
        continue
    items = line.split('\t')
    feature = items[2]
    if feature in features:
        if feature in ("gene", "protein_coding_gene", "pseudogene", "pseudogenic_tRNA", "rRNA_gene", "RNA", "snoRNA_gene", "snRNA_gene", "tRNA_gene"):
            if items[8].find("Parent") != -1: continue
            for x in items[8].split(';'):
                if x.split('=')[0] == "ID":
                    name = x.split('=')[1].strip()    
            feature = "gene"            
        else:    
            for x in items[8].split(';'):
                if x.split('=')[0] == "Parent":
                    name = x.split('=')[1].strip()
            name = re.sub(r'T.', "", name)
            
        if name not in table:
            table[name] = {}
            table[name]["sense"] = items[6]
            for x in features:
                table[name][x] = []
        if feature not in table[name]:    
            table[name][feature] = []   
        table[name][feature].append((int(items[3]), int(items[4])))    

location = {key: {"five_prime_UTR": 0, "CDS": 0, "not_annotated": 0, "intron" :0, "three_prime_UTR": 0} for key in table.keys()}
not_annotated = 0

for line in polyA_file:
    items = line.strip().split(" ")
    pos = int(items[1])
    chrx = items[2]
    value = int(items[0])
    gene = items[4]
    sense = items[3]
    flag = False
    for start, end in table[gene]['five_prime_UTR']:
        if pos >= start and pos <= end:
            location[gene]['five_prime_UTR'] += 1
            flag = True
            print line.strip(), "5'UTR"
    for start, end in table[gene]['CDS']:
       if pos >= start and pos <= end:
            location[gene]['CDS'] += 1           
            flag = True  
            print line.strip(), "CDS"
    if not flag:        
        if table[gene]["three_prime_UTR"]:
            if table[gene]["sense"] == "+" and pos >= table[gene]["three_prime_UTR"][0][0] or table[gene]["sense"] == "-" and pos <= table[gene]["three_prime_UTR"][0][1]:
                location[gene]['three_prime_UTR'] += 1
                print line.strip(), "3'UTR"
            elif pos >= table[gene]["gene"][0][0] and pos <= table[gene]["gene"][0][1]:
                location[gene]['intron'] += 1
                print line.strip(), "intron"
            else:
                location[gene]['not_annotated'] += 1
                print line.strip(), "not_annotated"
        else:
                location[gene]['not_annotated'] += 1
                print line.strip(), "not_annotated"
         
           
# where the polyA are located
three_utr = 0
five_utr = 0
cds = 0
intron = 0
not_annotated = 0
count_poly = 0.0
for transcript, loc in location.items():
    three_utr += loc["three_prime_UTR"]
    five_utr += loc["five_prime_UTR"]
    cds += loc["CDS"]
    intron += loc["intron"]
    not_annotated += loc["not_annotated"]
    count_poly += sum(loc.values())
    #print transcript,loc["three_prime_UTR"],loc["five_prime_UTR"],loc["CDS"],loc["intron"],loc["not_annotated"] 
if __name__ == "__main__":
  pass
  sys.stdout.write("%d,%d,%d,%d,%d\n" % (three_utr, five_utr, cds, intron, not_annotated))


gff_file.close()
polyA_file.close()
