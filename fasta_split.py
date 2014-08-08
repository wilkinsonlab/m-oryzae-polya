from Bio import SeqIO
import sys

main = SeqIO.parse(open(sys.argv[1]),'fasta')
for i in range(int(sys.argv[2])):
	print "Writing: " + sys.argv[1] + "_" + str(int(i+1)	)
	f =  open(sys.argv[1] + "_" + str(int(i+1)), "w")
	seqs = []
	i = 0
	for seq in main:
        	seqs.append(seq)
        	if len(seqs) > 100000 :
                        SeqIO.write(seqs, f, "fasta")
                        seqs = []
		i += 1
		if i > int(sys.argv[2]):
			SeqIO.write(seqs, f, "fasta")
			seqs = []
			break
	SeqIO.write(seqs, f, "fasta")
	seqs = []
	


