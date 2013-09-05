from Bio import SeqIO
import sys
import re
import os.path

output = open(sys.argv[2], 'w')
species = {}
for record in SeqIO.parse(sys.argv[1], "fasta"):
    sp = re.search(r'\[[^>]+\]$', record.description)
    if sp == None:
        continue
    sp = sp.group(0)
    sp = sp.replace("[", "")
    sp = sp.replace("]", "")
    sp = sp[:sp.find(" ")][:4] + "_" + sp[sp.find(" ") + 1:][:4]
    if not species.has_key(sp):
        output.write("> " + sp + "\n")
        output.write(str(record.seq) + "\n")
        print "s/" + sp + "/" + re.search("gi\|[0-9]+", record.description).group(0).replace("|", "\|") + "/g"
        species[sp] = True
