import sys

from Bio.Blast import NCBIXML
if len(sys.argv) > 1:
    blast = NCBIXML.parse(open(sys.argv[1], 'rU'))
else:
    blast = NCBIXML.parse(sys.stdin)

for record in blast:
    for align in record.alignments:
        if (align.hsps[0].frame[0] >= 0) and (align.hsps[0].frame[1] >= 0):
            print record.query, "\t", align.title, "\t", align.hsps[0].expect
            break
