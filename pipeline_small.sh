

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
for f in `ls ../[Wer]*sorted.bam.uniq`; do v=${f/..\//};   htseq-count -a 0 -s yes -r pos -f bam -t ncRNA -i ID $f ../Magnaporthe_oryzae.MG8.25.ncrna.gff3 > ${v/sorted.bam/ncrna.count} & done
# or 
for f in `ls ../[Wer]*sorted.bam.uniq`; do v=${f/..\//};   htseq-count -a 0 -s yes -r pos -f bam -t exon -i gene_id $f ../Magnaporthe_oryzae.MG8.25.gtf > ${v/sorted.bam/all.count} & done
#...now use DESeq2

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
~/Downloads/trinityrnaseq-2.0.2/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix trans.counts.matrix --method edgeR --samples_file ../../desc.txt ; ~/Downloads/trinityrnaseq-2.0.2/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix genes.counts.matrix --method edgeR --samples_file ../../desc.txt


### clusters diffential expression
for f in `ls ../*sorted.bam.uniq`; do v=${f/..\//}; v=${v/sorted.bam.uniq/cov}; echo $v ; genomeCoverageBed -ibam $f -g genome.txt -d > $v & done
for f in `ls ../*sorted.bam.uniq`; do v=${f/..\//}; v=${v/sorted.bam.uniq/cov} ;  l=`wc -l $f | cut -f 1 -d " "`; awk -v l=$l '{printf( "%s\t%d\t%f\n", $1, $2, $3*1000/l)}' < $v  > $v.norm & done
for f in `ls *cov`; do python /media/marco/Elements/m-oryzae-polya/cluster.py $f  > ${f/cov/bed} & done
cat *.bed | sort -k1,1 -k2,2n | bedtools merge -i - | awk '{print $1,"marco","peak",$2,$3,".",".",".","ID=peak_"++count"_"}' | sed 's/ /\t/g' > peaks.gff3
#...now use DESeq2

### sequences differential expression
#*fasta.collapsed.count were created after collapsing fastq sequences
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls ../*fasta.collapsed.count`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
cat $tmp > seq.count
cut -f 1,2 seq.count > exp_1.count;cut -f 1,3 seq.count > exp_2.count;cut -f 1,4 seq.count > exp_3.count; cut -f 1,5 seq.count > rbp_1.count; cut -f 1,6 seq.count > rbp_2.count; cut -f 1,7 seq.count > rbp_3.count; cut -f 1,8 seq.count > wt_1.count; cut -f 1,9 seq.count > wt_2.count; cut -f 1,10 seq.count > wt_3.count
#...now use DESeq2


# get info, in db folder
for f in `ls ../*fasta.collapsed`;
do
rm "_"${f/..\/}
bowtie -S -p 8 -v 0 -k 1 ../db/ncrna -f $f --un _un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 0 -k 1 ../db/retro  -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 0 -k 1 ../db/rrna  -f __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 0 -k 1 ../db/cds -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 0 -k 1 ../db/utr  -f __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 0 -k 1 ../db/unspliced  -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
bowtie -S -p 8 -v 0 -k 1 ../magna -f __un --un _unknown 1> /dev/null 2> _res; grep reported _res >> "_"${f/..\/}
done 

# for cluster etc...
for f in `ls *_vs_*.fa`;
do
rm _*	
bowtie2  -f $f --un=_un -x ../../db/bowtie2/ncrna --local -k 1 --quiet --no-head   | awk '{if($2!=4)print $0}'  > _ncrna_out 
bowtie2  -f _un --un=__un -x ../../db/bowtie2/retro --local -k 1 --quiet --no-head   | awk '{if($2!=4)print $0}'  > _retro_out 
bowtie2  -f __un --un=_un -x ../../db/bowtie2/utr --local -k 1 --quiet --no-head   | awk '{if($2!=4)print $0}' > _cds_out 
bowtie2  -f _un --un=__un -x ../../db/bowtie2/cds --local -k 1 --quiet --no-head   | awk '{if($2!=4)print $0}'  > _utr_out 
bowtie2  -f __un --un=_un -x ../../db/bowtie2/unspliced --local -k 1 --quiet --no-head   | awk '{if($2!=4)print $0}' > _intron_out
bowtie2  -f _un --un=__un -x ../../db/bowtie2/magna --local -k 1 --quiet --no-head   | awk '{if($2!=4)print $0}'  > _intergenic_out
cmsearch --cpu 4 -E 1e-2 --tblout _intergenic_rfam --noali ~/Downloads/Rfam.cm _un > /dev/null &
cmsearch --cpu 4 -E 1e-2 --tblout _unaligned_rfam --noali ~/Downloads/Rfam.cm __un > /dev/null; 
echo $f >  "_"$f
for g in _ncrna_out _retro_out _cds_out _utr_out _intron_out _intergenic_out; do echo ${g/_out/} > "_"$g; awk '{if($2==0)s="+"; else s="-"; print $1,$3,s}' < $g >> "_"$g ; done
echo "_intergenic_rfam" > __intergenic_rfam
if [ -e _intergenic_rfam ] ; then  sed -i -e 's/^#.*//' -e '/^$/d'  _intergenic_rfam; cat _intergenic_rfam >> __intergenic_rfam; fi
echo "_unaligned_rfam" > __unaligned_rfam
if [ -e _unaligned_rfam ] ; then  sed -i -e 's/^#.*//' -e '/^$/d'  _unaligned_rfam; cat _unaligned_rfam >> __unaligned_rfam; fi
paste  __ncrna_out __retro_out __cds_out __utr_out __intron_out __intergenic_out >> "_"$f
cat __intergenic_rfam >> "_"$f
cat __unaligned_rfam >> "_"$f
done


