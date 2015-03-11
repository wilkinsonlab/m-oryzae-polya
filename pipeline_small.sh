

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
	bowtie -S -M 1 -v 0 -p 8 --best --strata --al ${f/.trimmed.x/.al} --max  ${f/.trimmed.x/.max} --un ${f/.trimmed.x/.un}  magna $f | samtools view -Sh -F 4 - >  ${f/.*/.sam}
done

# sort and index
for f in `ls *.sam`; do
	samtools view -bhS $f | samtools sort - ${f/.sam/.sorted}
	samtools index ${f/.sam/.sorted.bam}
done

# create bedgraphs
for f in `ls *sorted.bam`; do genomeCoverageBed -bg -scale RPM -ibam $f > ${f/.sorted.bam/.bedgraph} & done

# converting aligned reads to fasta
for f in `ls *fastq.al`; do 
	fastq_to_fasta -i $f -o ${f/fastq/fasta}; 
done

### annotation diffential expression
for f in `ls ../[Wer]*sorted.bam`; do v=${f/..\//};   htseq-count -a 0 -s yes -r pos -f bam -t ncRNA -i ID $f ../Magnaporthe_oryzae.MG8.25.ncrna.gff3 > ${v/sorted.bam/count} & done
# or 
for f in `ls ../[Wer]*sorted.bam`; do v=${f/..\//};   htseq-count -a 0 -s yes -r pos -f bam -t exon -i gene_id $f ../Magnaporthe_oryzae.MG8.25.gff3 > ${v/sorted.bam/count} & done

### transcripts assembly diffential expression
# in diff_expr dir
# join all fastq into a single fasta
cat ../*fastq.trimmed.x | fastq_to_fasta -o all_reads.fa
# assemble all the reads
~/Downloads/trinityrnaseq-2.0.2/Inchworm/bin/inchworm --reads all_reads.fa --run_inchworm -K 18 -L 18 --num_threads 8 > all_reads.fa.inch
# prepare and run rsem
rsem-prepare-reference --bowtie all_reads.fa.inch all
for f in `ls ../*fastq.trimmed.x`; do v=${f/..\//}; rsem-calculate-expression --bowtie-n 0 --seed-length 25 -p 8 $f all ${v/.fastq.al/}; done
# DE on transcriptome
s=`ls *isoforms.results`; ~/Downloads/trinityrnaseq-2.0.2/util/abundance_estimates_to_matrix.pl --est_method RSEM  --out_prefix trans $s
s=`ls *genes.results`; ~/Downloads/trinityrnaseq-2.0.2/util/abundance_estimates_to_matrix.pl --est_method RSEM  --out_prefix genes $s
~/Downloads/trinityrnaseq-2.0.2/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix trans.counts.matrix --method edgeR --samples_file ../desc.txt ; ~/Downloads/trinityrnaseq-2.0.2/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix genes.counts.matrix --method edgeR --samples_file ../desc.txt


### peaks diffential expression
for f in `ls ../[Wer]*.sorted.bam`; do bedtools bamtobed -i $f | awk '{print $1, $2+2, $3, $4, $5, $6}' | sed 's/ /\t/g' > ${f/.sorted.bam/.bed}; done
mv ../[Wer]*bed .
for  f in `ls *.bed`; do macs14 -t $f -g 40949933 -n $f; done
cat *_peaks.bed | sort -k1,1 -k2,2n | bedtools merge -i - | awk '{print $1,"marco","peak",$2+1,$3,".",".",".","ID=peak_"++count"_"}' | sed 's/ /\t/g' > peaks.gff3
for f in `ls ../[Wer]*sorted.bam`; do v=${f/..\//};   htseq-count -a 0 -s no -r pos -f bam -t peak -i ID $f peaks.gff3 > ${v/sorted.bam/peak_count} & done
#...now use DESeq2

### sequences differential expression
#*fasta.collapsed.count were created after collapsing fastq sequences
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls ../*fasta.collapsed.count`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
cat $tmp > seq.count
cut -f 1,2 seq.count > exp_1.count;cut -f 1,3 seq.count > exp_2.count;cut -f 1,4 seq.count > exp_3.count; cut -f 1,5 seq.count > rbp_1.count; cut -f 1,6 seq.count > rbp_2.count; cut -f 1,7 seq.count > rbp_3.count; cut -f 1,8 seq.count > wt_1.count; cut -f 1,9 seq.count > wt_2.count; cut -f 1,10 seq.count > wt_3.count
#...now use DESeq2


# get info, in db folder
for f in `ls ../*fastq.trimmed.x`;
do
rm "_"${f/..\/}
bowtie -S -p 8 -v 0 -k 1 ncrna $f --un _un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 0 -k 1 cds _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 0 -k 1 utr __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 0 -k 1 unspliced _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 2 -k 1 rrna __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 2 -k 1 retro _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 0 -k 1 ../magna __un --un _unknown 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
cut -done 

for f in `ls ../*fasta.collapsed`;
do
rm "_"${f/..\/}
bowtie -S -p 8 -v 3 -k 1 ncrna -f $f --un _un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 3 -k 1 retro  -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 3 -k 1 rrna  -f __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 3 -k 1 cds -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 3 -k 1 utr  -f __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 3 -k 1 unspliced  -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 3 -k 1 ../magna -f __un --un _unknown 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
done 

# assembly aligned reads with inchworm
#for f in `ls *fasta.al`; do
#	~/Downloads/trinityrnaseq-2.0.2/Inchworm/bin/inchworm --reads $f --run_inchworm -K 18 -L 18 --num_threads 8 |  awk '{if(substr($1,1,1) == ">") print ">s"++count; else print $1}' > ${f/fasta.al/inch.fa}
#done

# alignes assembled reads back to the genome
#for f in `ls *inch.fa`; do 
#	gmap -D . -d Magnaporthe_oryzae.MG8.25.dna.genome.fa.gmap -f 3 -n 1  -t 8 -B 5 --nosplicing $f > ${f/inch.fa/inch.gff3} 2> /dev/null; 
#done

# merge all the alignments and create an annotation
#gt gff3 -sort *inch.gff3 | bedtools merge -d 51 -s -i - | awk '{print $1,"marco","smallrna",$2+1,$3,".",$4,".","ID=small_"++count"_"}' | sed 's/ /\t/g' > annotation.gff3

# extract only unique alignments, this is for htseq
#for f in `ls *.sorted.bam`; do
#	samtools view -h $f | awk '{if($5==255 || substr($0, 1, 1)=="@") print $0}' | samtools view -bSh - > $f".uniq"
#done

# htseq-count
#for f in `ls *.sorted.bam.uniq`; do
#	htseq-count -a 0 -s yes -r pos -f bam -t smallrna -i ID  $f annotation.gff3 > $f".count"
#done
