import sys, re
import math
import pysam

gff_file = open(sys.argv[1], 'r')
bam_file = pysam.Samfile(sys.argv[2], "rb")
out_file = open(sys.argv[3], 'w')
notout_file = open(sys.argv[4], 'w')
offset = int(sys.argv[5])

# RULES FOR 3' EXTRA SPACE
# WE USE THE "GENE" FEATURE END, BECAUSE IS THE ONLY ONE ALWAYS PRESENT IN ANNOTATION
# WE CONSIDER ALl THE INTERGENIC SPACE BETWEEN 3' OF A TRANSCRIPT AND 5' OF THE NEXT ONE UP TO "SPAN" NUCLEOTIDES
# IF 3' OVERLAPS WITH 5' WE DO NOT ADD EXTRA SPACE TO THE ANNOTATED 3' END

span = 400

table = {}
genes_sense = {}
genes_antisense = {}
utr = {}

for line in gff_file:
    if line[0] == '#' or line[0] == '\n':
        continue
    items = line.split('\t')

    # version 18
    #if items[2] == "gene":
    # version 21
    if items[2] in ("gene", "protein_coding_gene", "pseudogene", "pseudogenic_tRNA", "rRNA_gene", "RNA", "snoRNA_gene", "snRNA_gene", "tRNA_gene"):
        # version 21
        if items[8].find("Parent") != -1: continue
        for x in items[8].split(';'):
            if x.split('=')[0] == "ID":
                transcript = x.split('=')[1].strip()
        chrx = items[0]
        sense = items[6]
        start = int(items[3])
        end = int(items[4])

        if not table.has_key(chrx):
            table[chrx] = {}

        if not table[chrx].has_key(transcript):
            table[chrx][transcript] = [None, None, None]

        if not genes_sense.has_key(chrx):
            genes_sense[chrx] = []

        if not genes_antisense.has_key(chrx):
            genes_antisense[chrx] = []

        if sense == '+':
            genes_sense[chrx].append((transcript, start, end))
        elif sense == '-':
            genes_antisense[chrx].append((transcript, start, end))
       
def goo(sorted_list, i, sense):
    (transcript, start, end) = sorted_list[i]
    if sense == '+':
        # version 18
#         dist = sorted_list[i + 1][1] - end
#         if dist > span:
#             table[chrx][transcript] = [start, end + span, sense]
#         elif dist <= span and dist >= 0:
#             table[chrx][transcript] = [start, sorted_list[i + 1][1] - 1, sense]
#         elif dist < 0:
#             table[chrx][transcript] = [start, end, sense]        
        # version 21
        table[chrx][transcript] = [start, end + span, sense] 
    elif sense == '-':
        # version 18
#         dist = start - sorted_list[i - 1][2]
#         if dist > span:
#             table[chrx][transcript] = [start - span, end, sense]
#         elif dist <= span and dist >= 0:
#             table[chrx][transcript] = [sorted_list[i - 1][2] + 1, end, sense]
#         elif dist < 0:
#             table[chrx][transcript] = [start, end, sense]
        # version 21
        if start - span < 0: 
            table[chrx][transcript] = [0, end, sense]
        else:
            table[chrx][transcript] = [start - span, end, sense]    

for chrx, v in genes_sense.items():
    sorted_list_sense = sorted(v, key=lambda x: x[2])
    for i in range(1, len(sorted_list_sense) - 1):
        goo(sorted_list_sense, i, '+')
    goo(sorted_list_sense, 0, '+')
    goo(sorted_list_sense, -1, '+')

for chrx, v in genes_antisense.items():
    sorted_list_antisense = sorted(v, key=lambda x: x[1])
    for i in range(1, len(sorted_list_antisense) - 1):
        goo(sorted_list_antisense, i, '-')
    goo(sorted_list_antisense, 0, '-')
    goo(sorted_list_antisense, -1, '-')

block_span = 10000

table_blocks = dict.fromkeys(table.keys())
for chrx in table.keys():
    table_blocks[chrx] = dict((block, [])
                              for block in range(0, 10000000, block_span))
for chrx in table.keys():
    for transcript, (start, end, sense) in table[chrx].items():
        block_start = int(math.floor(start / float(block_span)) * block_span)
        block_end = int(math.floor(end / float(block_span)) * block_span)
        for block in range(block_start, block_end + block_span, block_span):
            table_blocks[chrx][block].append((transcript, start, end, sense))

count = 0
mapped = 0
for read_copy in bam_file.fetch():
    read = read_copy
    if read.is_reverse:
        read.pos += offset
    else:
        read.pos -= offset
    chrx = bam_file.getrname(read.tid)
    correct = None
    block = int(math.floor((read.pos, read.pos + read.rlen)[
                read.is_reverse] / float(block_span)) * block_span)
    for (transcript, start, end, sense) in table_blocks[chrx][block]:
        if not read.is_reverse and sense == '-' and (read.pos >= start and read.pos <= end) or \
                read.is_reverse and sense == '+' and (read.pos + read.rlen >= start and read.pos + read.rlen <= end):
                if correct == None:
                    correct = (transcript, start, end, sense)
                else:
                    if sense == '+' and abs(end - read.pos + read.rlen) < abs(correct[2] - read.pos + read.rlen):
                        correct = (transcript, start, end, sense)
                    elif sense == '-' and abs(start - read.pos) < abs(correct[1] - read.pos):
                        correct = (transcript, start, end, sense)

    if correct != None:
        out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.qname, chrx, read.pos, read.pos + read.rlen, (
            '+', '-')[read.is_reverse], correct[0], correct[1], correct[2], correct[3]))
        mapped += 1
    else:
        notout_file.write("%s\t%s\t%s\t%s\t%s\n" % (read.qname, chrx, read.pos, read.pos + read.rlen, (
            '+', '-')[read.is_reverse]))
    count += 1
    #if count % 1000000 == 0:
    #    print "Processed " + str(count) + " reads"


print "Assigned %d out of %d reads, %.2f %%" % (mapped, count, mapped / float(count) * 100)

gff_file.close()
bam_file.close()
out_file.close()
