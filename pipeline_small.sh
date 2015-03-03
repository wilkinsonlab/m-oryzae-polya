

# trimming 3' adaptor
for f in `ls WT*.fastq`; 
do
	fastx_clipper -l 26 -a TGGAATTCTCGGGTGCCAAGG -i $f -o $f".trimmed" &
done

# trimming HD ->4 4<- nucleotides
for f in `ls *.fastq.trimmed`;
do
python -c "
import sys
from Bio import SeqIO
fasta_file = sys.argv[1]  
result_file = sys.argv[2] 
fasta_sequences = SeqIO.parse(open(fasta_file),'fastq')
f =  open(result_file, 'w')
seqs = []
for seq in fasta_sequences:
  edit = seq[4:-4]
  seqs.append(edit)
  if len(seqs) > 100000:
      SeqIO.write(seqs, f, 'fastq')
      seqs = []
SeqIO.write(seqs, f, 'fastq')
" $f $f".x" &
done

# aligning
for f in `ls *fastq.trimmed.x`;
do
	bowtie -S -M 1 -v 0 -p 8 --best --strata magna $f > ${f/.*/.sam}
done


# converting to fasta
for f in `ls *.trimmed.x`; do 
	fastq_to_fasta -i $f -o ${f/fastq/fasta}; 
done

# assembly with inchworm
for f in `ls *fasta.trimmed.x`; do
	~/Downloads/trinityrnaseq-2.0.2/Inchworm/bin/inchworm --reads $f --run_inchworm -K 18 -L 18 --num_threads 8 |  awk '{if(substr($1,1,1) == ">") print ">s"++count; else print $1}' > ${f/fasta.trimmed.x/inch.fa}
done

