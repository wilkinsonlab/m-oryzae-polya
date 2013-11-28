import sys
import polyA_localization

# A GENE CAN BE APA IN MANY WAYS

# how many genes have each kind of APA
apa_three_utr = 0
apa_five_utr = 0
apa_cds = 0
apa_intron = 0
count_poly = 0.0
count_apa = 0.0
for gene, loc in polyA_localization.location.items():
    if loc["three_prime_UTR"] > 1:
        apa_three_utr += 1
    if loc["three_prime_UTR"] >= 1 and loc["five_prime_UTR"] >= 1:
        apa_five_utr += 1
    if loc["three_prime_UTR"] >= 1 and loc["CDS"] >= 1:
        apa_cds += 1
    if loc["three_prime_UTR"] >= 1 and loc["intron"] >= 1:
        apa_intron += 1
    count_poly += sum(loc.values())
    if sum(loc.values()) > 1:
        count_apa += 1

sys.stdout.write("%d,%d,%d,%d\n" % (apa_three_utr, apa_five_utr, apa_cds, apa_intron))
