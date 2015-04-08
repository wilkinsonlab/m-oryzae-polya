# trimming 3' adaptor, quality and artifacts filtering
for f in `ls *.fastq`; 
do
	fastx_clipper -c -l 26 -a TGGAATTCTCGGGTGCCAAGG -i $f | fastq_quality_filter -q 30 -p 70 | fastx_artifacts_filter -o $f".trimmed" &
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

### bowtie1
# aligning not unique
for f in `ls *fastq.trimmed.x`;
do
	bowtie -S -M 1 -v 0 -p 8 --best --strata --al ${f/.trimmed.x/.al} --max  ${f/.trimmed.x/.max} --un ${f/.trimmed.x/.un}  db/bowtie/genome $f | samtools view -Sh -F 4 - >  ${f/.*/.sam}
done
# aligning unique
for f in `ls *fastq.trimmed.x`;
do
	bowtie -S -m 1 -v 0 -p 8 --best --strata db/bowtie/genome $f | samtools view -Sh -F 4 - >  ${f/.*/.uniq.sam}
done


### bowtie2 --end-to-end OR --local
for f in `ls *fastq.trimmed.x`; 
do 
	bowtie2 -p 8 --end-to-end --al genomic_based/bowtie2_end_to_end/${f/.trimmed.x/.al}  --un genomic_based/bowtie2_end_to_end/${f/.trimmed.x/.un}  -x db/bowtie2/genome -q $f | samtools view -bSh -F 4 - | samtools sort - genomic_based/bowtie2_end_to_end/${f/.*/.sorted}; 
done

# sort and index
for f in `ls *.sam`; do
	samtools view -bhS $f | samtools sort - ${f/.sam/.sorted}
	samtools index ${f/.sam/.sorted.bam}
done

# create bedgraphs
for f in `ls *sorted.bam`; do genomeCoverageBed -bg -ibam $f -g genome.txt > ${f/.sorted.bam/.bedgraph} & done

# collapse reads
for f in `ls *fastq.trimmed.x`; do fastx_collapser -i $f -o $f.collapsed & done

### transcripts assembly diffential expression
# in diff_expr dir
# join all fastq into a single fasta
cat ../*fastq.trimmed.x | fastq_to_fasta -o all_reads.fa
# assemble all the reads
~/Downloads/trinityrnaseq-2.0.2/Inchworm/bin/inchworm --reads all_reads.fa --run_inchworm -K 18 -L 18 --num_threads 8 > all_reads.fa.inch
# prepare and run rsem
~/Downloads/rsem-1.2.19/rsem-prepare-reference --bowtie all_reads.fa.inch all
for f in `ls ../*fastq.trimmed.x`; do v=${f/..\//}; ~/Downloads/rsem-1.2.19/rsem-calculate-expression --bowtie-n 0 --seed-length 18 -p 8 $f all ${v/.fastq.trimmed.x/}; done
# DE on transcriptome
s=`ls *isoforms.results`; ~/Downloads/trinityrnaseq-2.0.2/util/abundance_estimates_to_matrix.pl --est_method RSEM  --out_prefix trans $s
s=`ls *genes.results`; ~/Downloads/trinityrnaseq-2.0.2/util/abundance_estimates_to_matrix.pl --est_method RSEM  --out_prefix genes $s
~/Downloads/trinityrnaseq-2.0.2/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix trans.counts.matrix --method edgeR --samples_file ../desc.txt 
~/Downloads/trinityrnaseq-2.0.2/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix genes.counts.matrix --method edgeR --samples_file ../desc.txt

### transcripts assembly diffential expression with DESeq2
cat ../W*fastq.trimmed.x | fastq_to_fasta -o all_reads_WT.fa
cat ../e*fastq.trimmed.x | fastq_to_fasta -o all_reads_exp5.fa
cat ../r*fastq.trimmed.x | fastq_to_fasta -o all_reads_rbp35.fa
for f in `ls *fa`; do ~/Downloads/trinityrnaseq-2.0.2/Inchworm/bin/inchworm --reads $f --run_inchworm -K 32 -L 32 --minKmerCount 2 --num_threads 8 > $f.inch; done
cat *inch | fasta_formatter | fastx_collapser -o all_reads.fa.inch
bowtie-build all_reads.fa.inch all
for f in `ls ../*fastq.trimmed.x`; do v=${f/..\//}; bowtie -S -p 8 -v 0 -m 1 --sam-nohead  all $f  | awk '{if($2==0)print $3}' | sort | uniq -c | awk '{print $2"\t"$1}' > "_"$v ; done
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls _*x`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
cat $tmp > _count
cut -f 1,2 _count > exp5_1.count;cut -f 1,3 _count > exp5_2.count;cut -f 1,4 _count > exp5_3.count; cut -f 1,5 _count > rbp35_1.count; cut -f 1,6 _count > rbp35_2.count; cut -f 1,7 _count > rbp35_3.count; cut -f 1,8 _count > WT_1.count; cut -f 1,9 _count > WT_2.count; cut -f 1,10 _count > WT_3.count
Rscript ../diff.R WT_1.count WT_2.count  WT_3.count exp5_1.count exp5_2.count exp5_3.count WT_vs_EXP5.csv
Rscript ../diff.R WT_1.count WT_2.count  WT_3.count rbp35_1.count rbp35_2.count rbp35_3.count WT_vs_RBP35.csv

### sequences differential expression
for f in `ls *.fastq.trimmed.x`; do fastx_collapser -i $f -o ${f/fastq/fasta}.collapsed ; done
for f in `ls *fasta.trimmed.x.collapsed`; do awk '{if (substr($0, 1, 1) == ">"){ split($0, arr, "-"); val=arr[2];} else { print $0"\t"val; val="-" }  }' < $f > $f.count; done
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls ../*fasta.collapsed.count`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
cat $tmp > seq.count
cut -f 1,2 seq.count > exp_1.count;cut -f 1,3 seq.count > exp_2.count;cut -f 1,4 seq.count > exp_3.count; cut -f 1,5 seq.count > rbp_1.count; cut -f 1,6 seq.count > rbp_2.count; cut -f 1,7 seq.count > rbp_3.count; cut -f 1,8 seq.count > wt_1.count; cut -f 1,9 seq.count > wt_2.count; cut -f 1,10 seq.count > wt_3.count
Rscript ../diff.R WT_1.count WT_2.count  WT_3.count exp5_1.count exp5_2.count exp5_3.count WT_vs_EXP5.csv
Rscript ../diff.R WT_1.count WT_2.count  WT_3.count rbp35_1.count rbp35_2.count rbp35_3.count WT_vs_RBP35.csv

### HD adapter differential expression
for f in `ls ../*fasta.trimmed`; do v=${f/..\//}; grep -v ">" $f | awk '{print substr($0,1,4)}' | sort | uniq -c | awk '{print $2"\t"$1}' > ${v/fasta.trimmed/5prime.count}; done
for f in `ls ../*fasta.trimmed`; do v=${f/..\//}; grep -v ">" $f | awk '{print substr($0,length($0)-3,4)}' | sort | uniq -c | awk '{print $2"\t"$1}' > ${v/fasta.trimmed/3prime.count}; done
Rscript /media/marco/Elements/m-oryzae-polya/adapters_diff.R WT_1.5prime.count WT_2.5prime.count WT_3.5prime.count exp5_1.5prime.count exp5_2.5prime.count exp5_3.5prime.count WT_vs_EXP5.5prime.csv
Rscript /media/marco/Elements/m-oryzae-polya/adapters_diff.R WT_1.3prime.count WT_2.3prime.count WT_3.3prime.count exp5_1.3prime.count exp5_2.3prime.count exp5_3.3prime.count WT_vs_EXP5.3prime.csv
Rscript /media/marco/Elements/m-oryzae-polya/adapters_diff.R WT_1.5prime.count WT_2.5prime.count WT_3.5prime.count rbp35_1.5prime.count rbp35_2.5prime.count rbp35_3.5prime.count WT_vs_RBP35.5prime.csv
Rscript /media/marco/Elements/m-oryzae-polya/adapters_diff.R WT_1.3prime.count WT_2.3prime.count WT_3.3prime.count rbp35_1.3prime.count rbp35_2.3prime.count rbp35_3.3prime.count WT_vs_RBP35.3prime.csv



# extraction of adapter-related sequences
#sample oriente
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $6<0.05 && $6!="NA" && $3>0) print $1}' < WT_vs_RBP35_3prime.csv | sed 's/"//g'`; do grep "^"$f ../[Wr]*.fasta.trimmed.collapsed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_RBP35_3prime_up.fa & 
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $6<0.05 && $6!="NA" && $3<0) print $1}' < WT_vs_RBP35_3prime.csv | sed 's/"//g'`; do grep "^"$f ../[Wr]*.fasta.trimmed.collapsed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_RBP35_3prime_down.fa &
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $6<0.05 && $6!="NA" && $3>0) print $1}' < WT_vs_RBP35_5prime.csv | sed 's/"//g'`; do grep "^"$f ../[Wr]*.fasta.trimmed.collapsed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_RBP35_5prime_up.fa & 
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $6<0.05 && $6!="NA" && $3<0) print $1}' < WT_vs_RBP35_5prime.csv | sed 's/"//g'`; do grep "^"$f ../[Wr]*.fasta.trimmed.collapsed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_RBP35_5prime_down.fa &
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $6<0.05 && $6!="NA" && $3>0) print $1}' < WT_vs_EXP5_3prime.csv | sed 's/"//g'`; do grep "^"$f ../[We]*.fasta.trimmed.collapsed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_EXP5_3prime_up.fa & 
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $6<0.05 && $6!="NA" && $3<0) print $1}' < WT_vs_EXP5_3prime.csv | sed 's/"//g'`; do grep "^"$f ../[We]*.fasta.trimmed.collapsed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_EXP5_3prime_down.fa &
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $6<0.05 && $6!="NA" && $3>0) print $1}' < WT_vs_EXP5_5prime.csv | sed 's/"//g'`; do grep "^"$f ../[We]*.fasta.trimmed.collapsed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_EXP5_5prime_up.fa & 
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $6<0.05 && $6!="NA" && $3<0) print $1}' < WT_vs_EXP5_5prime.csv | sed 's/"//g'`; do grep "^"$f ../[We]*.fasta.trimmed.collapsed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_EXP5_5prime_down.fa &
#adapter oriented
 for f in `cat *5prime.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $6<0.05 && $6!="NA" ) print $1}' | sed 's/"//g'`; do  grep "^"$f ../[Wer]*.fasta.trimmed.collapsed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }' > "_"$f ; done
 for f in `cat *3prime.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $6<0.05 && $6!="NA" ) print $1}' | sed 's/"//g'`; do  grep $f"$" ../[Wer]*.fasta.trimmed.collapsed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }' > "__"$f ; done


### for the avr-regions (modify cluster.py)
for f in `ls ../[Wer]*sorted.bam`; do v=${f/..\//}; v=${v/sorted.bam/cov}; echo $v ; genomeCoverageBed -ibam $f -g ../genome.txt -d > $v & done
for f in `ls *cov`; do python /media/marco/Elements/m-oryzae-polya/cluster.py $f  > ${f/cov/bed} & done
cat *.bed | sort -k1,1 -k2,2n | bedtools merge -i - | awk '{print $1,"marco","region",$2,$3,".",".",".","ID=region_"++count"_"}' | sed 's/ /\t/g' > regions.gff3
 
 


# get info unique
rm _*
for f in `ls *[123].fasta.trimmed.x`;
do
bowtie -S -p 8 -v 1 -k 1 db/bowtie/ncrna -f $f --un _un 1> /dev/null 2> _res; grep reported _res >> "_"$f
bowtie -S -p 8 -v 1 -k 1 db/bowtie/rrna  -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"$f
bowtie -S -p 8 -v 1 -k 1 db/bowtie/retro  -f __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"$f
bowtie -S -p 8 -v 1 -k 1 db/bowtie/transcripts  -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"$f
bowtie -S -p 8 -v 1 -k 1 db/bowtie/unspliced  -f __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"$f
bowtie -S -p 8 -v 1 -k 1 db/bowtie/genome -f _un --un _unknown 1> /dev/null 2> _res; grep reported _res >> "_"$f
done 
rm _*
for f in `ls *.fasta.trimmed.x.collapsed`; do 
bowtie2  -p 2  -x db/bowtie2/ncrna -f $f --un _un 1> /dev/null 2> _res; grep exactly _res >> "_"$f; 
bowtie2  -p 2  -x db/bowtie2/rrna  -f _un --un __un 1> /dev/null 2> _res; grep exactly _res >> "_"$f; 
bowtie2  -p 2  -x db/bowtie2/retro -f __un --un _un 1> /dev/null 2> _res; grep exactly _res >> "_"$f; 
bowtie2  -p 2  -x db/bowtie2/transcripts  -f _un --un __un 1> /dev/null 2> _res; grep exactly _res >> "_"$f; 
bowtie2  -p 2  -x db/bowtie2/unspliced  -f __un --un _un 1> /dev/null 2> _res; grep exactly _res >> "_"$f; 
bowtie2  -p 2  -x db/bowtie2/genome -f _un --un _unknown 1> /dev/null 2> _res; grep exactly _res >> "_"$f; 
done
#get info expression


# for cluster, milRNA, assemblies
rm _*
for f in `ls all_siRNA.cluster.fa `;
do	
db_dir="../../../../db"
touch 	_ncrna_out _rrna_out _retro_out _transcripts_out _gene_out _intergenic_out
blastn -dust no  -num_threads 4  -task  blastn -query $f -db $db_dir/ncrna.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2> /dev/null| awk '{if ($4/$5 > 0.9 || $4/$6 > 0.9)  print $0}' > _ncrna_out
blastn -dust no  -num_threads 4  -task  blastn -query $f -db $db_dir/rrna.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2>  /dev/null | awk '{if ($4/$5 > 0.9 || $4/$6 > 0.9) print $0}' > _rrna_out
blastn -dust no  -num_threads 4  -task  blastn -query $f -db $db_dir/retro.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2> /dev/null | awk '{if ($4/$5 > 0.9 || $4/$6 > 0.9) print $0}' > _retro_out
blastn -dust no  -num_threads 4  -task  blastn -query $f -db $db_dir/transcripts.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2> /dev/null  | awk '{if ($4/$5 > 0.9 || $4/$6 > 0.9) print $0}' > _transcripts_out
blastn -dust no  -num_threads 4  -task  blastn -query $f -db $db_dir/unspliced.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2> /dev/null |  awk '{if ($4/$5 > 0.9 || $4/$6 > 0.9) print $0}' > _gene_out
blastn -dust no  -num_threads 4  -task  blastn -query $f -db $db_dir/genome.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2> /dev/null  | awk '{if ($4/$5 > 0.9 || $4/$6 > 0.9) print $0}' > _intergenic_out
for g in `ls _*`; do sort -k1,1 -k12,12nr -k11,11n $g | sort -u -k1,1 --merge | sort -o $g; done
## print numbers
#for g in _ncrna_out _rrna_out _retro_out _transcripts_out _gene_out _intergenic_out; do awk -v g=$g '{print g"\t"$0}' < $g ; done |  awk '{if ($2 in arr == 0)  arr[$2]=$0}END{for (k in arr) print arr[k] }' | cut -f 1 | sort | uniq -c | sed -e 's/_out//' -e 's/_//' | awk '{print $2"\t"$1}'> "__"$f
## print details	
echo -e "\n"$f
python -c "
yes = []
for file in ('_ncrna_out', '_rrna_out', '_retro_out','_transcripts_out','_gene_out','_intergenic_out'):
  out = open('_'+file, 'w')
  for line in open(file):
    items = line.strip().split('\t')
    if items[0] not in yes:
      out.write(items[1] + '\n')
      yes.append(items[0])
  out.close()       
" 
for f in `ls __*`; do sort $f | uniq -c | sort -o $f; done
paste __*
done
## print numbers
#tmp=$(mktemp);tmp2=$(mktemp);for file in `ls __*`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
#cat $tmp


rm _*
for f in `ls *.fasta.trimmed.x.collapsed`;
do	
bowtie2 --end-to-end -x db/bowtie2/ncrna  -f $f --un _un --quiet --no-head -p 2  | awk '{if($2!=4)print $0}'  > _ncrna_out 
bowtie2 --end-to-end -x db/bowtie2/rrna -f _un --un __un --quiet --no-head -p 2   | awk '{if($2!=4)print $0}'  > _rrna_out 
bowtie2 --end-to-end -x db/bowtie2/retro -f __un --un _un --quiet --no-head -p 2   | awk '{if($2!=4)print $0}' > _retro_out 
bowtie2 --end-to-end -x db/bowtie2/transcripts -f _un --un __un --quiet --no-head -p 2   | awk '{if($2!=4)print $0}'  > _transcripts_out 
bowtie2 --end-to-end -x db/bowtie2/unspliced  -f __un --un _un --quiet --no-head -p 2   | awk '{if($2!=4)print $0}' > _gene_out
bowtie2 --end-to-end -x db/bowtie2/genome -f _un --al _intergenic --un=_unaligned  --quiet --no-head -p 2   | awk '{if($2!=4)print $0}'  > _intergenic_out
for g in _ncrna_out _rrna_out _retro_out _transcripts_out _gene_out _intergenic_out ;
do
  awk '{split($1, arr, "-"); count+=arr[2]}END{print count}' < $g > "_"$g
done
echo -ne $f","
paste  -d "," __ncrna_out __rrna_out __retro_out __transcripts_out __gene_out __intergenic_out
done

echo -ne $f" "  
#for g in _ncrna_out _retro_out _cds_out _utr_out _gene_out _intergenic_out; do echo ${g/_out/} > "_"$g; awk '{if($2==0)s="+"; else s="-"; print $1,$3,s}' < $g >> "_"$g ; done
for g in _ncrna_out _rrna_out _retro_out _transcripts_out _gene_out _intergenic_out; do  wc -l $g | cut -f 1 -d " " >> "_"$g ; done
paste  -d "," __ncrna_out __rrna_out __retro_out __transcripts_out __gene_out __intergenic_out
continue


#num=`grep -c ">" $f`
#for g in _ncrna_out _retro_out _cds_out _utr_out _gene_out _intergenic_out; do echo ${g/_out/} > "_"$g; wc -l $g | cut -f 1 -d " " | awk -v num=$num '{print $0/num}' >> "_"$g ; done
#echo -e "adapter\n"$f |  paste  -d "," - __ncrna_out __retro_out __cds_out __utr_out __gene_out __intergenic_out
done


#################### 
# CLUSTERS, BOTH GENOMIC AND TRANSCRIPTOME
#################### 
# get all the alignments
for f in `ls *fasta.trimmed.x.collapsed`;  do  
bowtie2 -p 8 -a --end-to-end  -x db/bowtie2/genome -f $f | samtools view -bSh -F 4 - | samtools sort - genomic_based/cluster/bowtie/${f/.*/.sorted};  
done

# compute weighted expression of reads
for f in WT_1  WT_2 WT_3 exp5_1 exp5_2 exp5_3 rbp35_1  rbp35_2 rbp35_3; 
do 
python /media/marco/Elements/m-oryzae-polya/weight_reads.py ../../../$f.fasta.trimmed.x.collapsed $f.sorted.bam > $f.weighted 
done

# calculate coverage on weighted reads
for f in `ls *.weighted`;
do
python -c "
import numpy, pysam, sys
chrom = {}
for line in open(sys.argv[1], 'r'):
  chrx, length = line.strip().split('\t')
  chrom[chrx] = numpy.array([0.0] * (int(length)+1))

i = 0
for line in open(sys.argv[2], 'r'):
  name, seq, rchrx, rstart, rend, val = line.strip().split('\t')  
  i += 1
  if i % 1000000 == 0: sys.stderr.write(str(i) + '\n')
  chrom[rchrx][int(rstart):int(rend)+1] += float(val)

for chrx, array in chrom.items():
  for pos, val in enumerate(array):
     print chrx + '\t' + str(pos) + '\t' + str(val)
" ../../../genome.txt $f > ${f/weighted/cov} &
done



# detect clusters
for f in `ls *cov`; do python /media/marco/Elements/m-oryzae-polya/cluster.py $f | sed 's/ /\t/g' > ${f/cov/bed}  done
cat *.bed | sort -k1,1 -k2,2n | bedtools merge -i - | awk '{print $1,"marco","cluster",$2,$3,".",".",".","ID=cluster_"++count"_"}' | sed 's/ /\t/g' > clusters.gff3
cat clusters.gff3 | bedtools getfasta -bed - -fi ../../../db/genome.fa -fo clusters.fa
cat clusters.gff3 | while read c x x s e x x x d; do  n=$c":"$((s-1))"-"$e ; sed  's/'$n'/'${d/ID=/}'/' clusters.fa > _o ;  mv _o clusters.fa; done

# assign weighted reads to clusters
for f in `ls *.weighted`;
do
python -c "
import numpy, pysam, sys, math
from Bio import SeqIO

clusters = []
reads = {}
for line in open(sys.argv[1], 'r'):
  items = line.strip().split('\t')
  clusters.append((items[8].replace('ID=', ''),  items[0], int(items[3]), int(items[4])))
  reads[items[8].replace('ID=', '')] = {}

#block_span = 100000
#chroms = ('supercont8.8', '7', '6', '5', '4', '3', '1', '2')
#table_blocks = dict.fromkeys(chroms)
#for chrx in chroms:
#   table_blocks[chrx] = dict((block, [])
#                              for block in range(0, 10000000, block_span))
#for chrx in chroms:
#    for cluster, chrx, start, end in clusters:
#        block_start = int(math.floor(start / float(block_span)) * block_span)
#        block_end = int(math.floor(end / float(block_span)) * block_span)
#        for block in range(block_start, block_end + block_span, block_span):
##            table_blocks[chrx][block].append((cluster, chrx, start, end))
### trascriptome
table_blocks = {}
for cluster, chrx, start, end in clusters:
  if not table_blocks.has_key(chrx): table_blocks[chrx] = []
  table_blocks[chrx].append((cluster, chrx, start, end)) 
###
	 

i = 0
for line in open(sys.argv[2], 'r'):
  name, seq, rchrx, rstart, rend, val = line.strip().split('\t')
  i += 1
  if i % 1000000 == 0: sys.stderr.write(str(i) + '\n')
### transcriptome
  if not table_blocks.has_key(rchrx): continue
  for (cluster, chrx, start, end) in table_blocks[rchrx]:
###    
#  block = int(math.floor(int(rstart) / float(block_span)) * block_span)
#  for (cluster, chrx, start, end) in table_blocks[rchrx][block]:
    if rchrx == chrx and int(rstart) >= start and int(rend) <= end:
       print cluster + '\t' + chrx + '\t' + str(start) + '\t' + str(end) + '\t' + name+ '\t' + seq + '\t' + rchrx + '\t' + rstart + '\t' + rend + '\t' + val 
       break
" clusters.gff3 $f > ${f/weighted/assign} &
done
wait


# compute clusters differential expression
for f in `ls *assign`; do awk '{arr[$1]+=$10}END{for (k in arr) print k"\t"arr[k]}' < $f > ${f/assign/expr} ; done
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls *expr`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
awk '{print $1"\t"int($2)}' < $tmp > exp5_1.expr; awk '{print $1"\t"int($3)}' < $tmp > exp5_2.expr; awk '{print $1"\t"int($4)}' < $tmp > exp5_3.expr;
awk '{print $1"\t"int($5)}' < $tmp > rbp35_1.expr; awk '{print $1"\t"int($6)}' < $tmp > rbp35_2.expr; awk '{print $1"\t"int($7)}' < $tmp > rbp35_3.expr;
awk '{print $1"\t"int($8)}' < $tmp > WT_1.expr; awk '{print $1"\t"int($9)}' < $tmp > WT_2.expr; awk '{print $1"\t"int($10)}' < $tmp > WT_3.expr;
Rscript ../../../../diff.R WT_1.expr WT_2.expr WT_3.expr exp5_1.expr exp5_2.expr exp5_3.expr WT_vs_EXP5.csv
Rscript ../../../../diff.R WT_1.expr WT_2.expr WT_3.expr rbp35_1.expr rbp35_2.expr rbp35_3.expr WT_vs_RBP35.expr.csv

#############
# for DEXSeq
#############
for f in `ls *assign`;
do
awk '{if (($4-$3)<=300) {split($5, arr, "_");  print $1":E"$6"_"arr[2]"\t"int($10)}}' < $f | sort -u > "_"$f 
done
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls _*assign`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
awk '{print $1"\t"int($8)}' < $tmp > _out1 ; awk '{print $1"\t"int($9)}' < $tmp > _out2 ; awk '{print $1"\t"int($10)}' < $tmp > _out3 ;
awk '{print $1"\t"int($2)}' < $tmp > _out4 ; awk '{print $1"\t"int($3)}' < $tmp > _out5 ; awk '{print $1"\t"int($4)}' < $tmp > _out6 ;
awk '{print $1"\t"int($5)}' < $tmp > _out7 ; awk '{print $1"\t"int($6)}' < $tmp > _out8 ; awk '{print $1"\t"int($7)}' < $tmp > _out9 ;

############# 
# siRNA discovery >300 discarted for MEMORY reasons
############# 
for g in WT_1  WT_2 WT_3 exp5_1 exp5_2 exp5_3 rbp35_1  rbp35_2 rbp35_3; 
do 
python -c "
import sys
class siRNA:
  def __init__(self, cluster, start, end):
    self.cluster = cluster
    self.start = start
    self.end = end
    self.expr = 0.0
    self.seqs = {}
    self.vals = {}
    self.senses = {}

siRNAs = {}    
i = 0
for line in open(sys.argv[1], 'r'):
  cluster, chrx, start, end, name, seq, rchrx, rstart, rend, val = line.strip().split('\t')
  i += 1
  if i % 1000000 == 0: sys.stderr.write(str(i) + '\n')
  if int(end)-int(start) > 300: continue
  if not siRNAs.has_key(cluster): siRNAs[cluster] = siRNA(cluster, start, end)
  flag = name.split('_')[1]
  s=str(bin(int(flag))).replace('0b', '')
  x='0' * (16-len(s)) + s
  if x[-5] == '0': 
    siRNAs[cluster].senses[name] = '+'
  else:
    siRNAs[cluster].senses[name] = '-' 
  siRNAs[cluster].expr += float(val)
  siRNAs[cluster].vals[name] = float(val)
  siRNAs[cluster].seqs[name] = seq 
  
for cluster, sirna in siRNAs.items():
  for name, val in sirna.vals.items():
    if val / sirna.expr > 0.5:
      print sirna.cluster + '\t' + str(sirna.start) + '\t' + str(sirna.end) + '\t' + str(sirna.expr) + '\t' + sirna.seqs[name] + '\t' + sirna.senses[name] + '\t' + str(val)
      break
" $g.assign > $g.siRNA 
done

#| awk -v count=$count '{norm=(($2/(40949933/1000))/(count/1000000)); if (norm > 0.000020) print $1}' | while read f; do  grep -m 1 $f clusters.gff3 | awk '{if($5-$4<=100) print $0}' ; done | sort > $g.siRNA
done

sort exp5_[123].siRNA | cut -f 1,2,3,5,6 | uniq -c | awk '{if ($1==3)print $2,$3,$4,$5,$6}'  | sed -e 's/ID=//' -e 's/ /\t/g' > exp5_siRNA
sort rbp35_[123].siRNA | cut -f 1,2,3,5,6 | uniq -c | awk '{if ($1==3)print $2,$3,$4,$5,$6}'  | sed -e 's/ID=//' -e 's/ /\t/g' > rbp35_siRNA
sort WT_[123].siRNA | cut -f 1,2,3,5,6 | uniq -c | awk '{if ($1==3)print $2,$3,$4,$5,$6}'  | sed -e 's/ID=//' -e 's/ /\t/g' > WT_siRNA
sort WT_siRNA rbp35_siRNA exp5_siRNA | uniq > all_siRNA
awk '{print ">"$1"\n"$4}' < exp5_siRNA > exp5_siRNA.fa
awk '{print ">"$1"\n"$4}' < rbp35_siRNA > rbp35_siRNA.fa
awk '{print ">"$1"\n"$4}' < WT_siRNA > WT_siRNA.fa
awk '{print ">"$1"\n"$4}' < all_siRNA > all_siRNA.fa
cut -f 1 all_siRNA | sort -u | grep -f - ../WT_vs_EXP5.expr.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3<0  ) print $1}' | sed 's/"//g' > _exp_down
cut -f 1 all_siRNA | sort -u | grep -f - ../WT_vs_EXP5.expr.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3>0  ) print $1}' | sed 's/"//g' > _exp_up
cut -f 1 all_siRNA | sort -u | grep -f - ../WT_vs_RBP35.expr.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3<0  ) print $1}' | sed 's/"//g' > _rbp_down
cut -f 1 all_siRNA | sort -u | grep -f - ../WT_vs_RBP35.expr.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3>0  ) print $1}' | sed 's/"//g' > _rbp_up
python /media/marco/Elements/m-oryzae-polya/fasta_extract.py all_siRNA.fa _exp_down WT_vs_EXP5.siRNA.down.fa
python /media/marco/Elements/m-oryzae-polya/fasta_extract.py all_siRNA.fa _exp_up WT_vs_EXP5.siRNA.up.fa
python /media/marco/Elements/m-oryzae-polya/fasta_extract.py all_siRNA.fa _rbp_down WT_vs_RBP35.siRNA.down.fa
python /media/marco/Elements/m-oryzae-polya/fasta_extract.py all_siRNA.fa _rbp_up WT_vs_RBP35.siRNA.up.fa  

### retrotransposons map
for f in `ls ../[Wer]*fastq.trimmed.x`; 
do 
	v=${f/..\//}; v=${v/fastq.trimmed.x/sam}; bowtie --best -S -v 0 -p 8 -m 1 ../db/bowtie/retro $f | awk '{if($2!=4)print $0}'  > $v; 
done

for f in `ls [Wer]*sam`; 
do 
grep "SN:.*LN:[0-9]*" -o < $f | sed -e 's/SN://' -e 's/LN://' > "_"$f
cat "_"$f | while read a b; 
do 
grep $a $f | grep -v "^@" | awk '{if($2==0) printf "%d\n", $4/10}' | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k 1,1 > "_"$f"_"$a"_pos" &
grep $a $f | grep -v "^@" | awk '{if($2==16) printf "%d\n", $4/10}' | sort | uniq -c | awk '{print $2"\t"(-$1)}' | sort -k 1,1 > "_"$f"_"$a"_neg" &
done
done

grep "SN:.*LN:[0-9]*" -o < exp5_IP17_GTAGAG_L006_R1.sam | sed -e 's/SN://' -e 's/LN://' | cut -f 1,2 > "__"$f
while read a b
do 
rm "__"$a; 
for i in `seq 1 $b`
do 
echo $i | awk '{printf "%d\t0\n", $0/10}' 
done | sort -u > "__"$a
done < __exp5_IP17_GTAGAG_L006_R1.sam

# for all the retros......
reset;tmp=$(mktemp);tmp2=$(mktemp);for file in `ls *retro5*`; do sort -k 1,1 $file -o $file ;  if [ -s "$tmp" ];   then   join -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";   else     cp "$file" "$tmp2";   fi;   cp "$tmp2" "$tmp"; done; cat $tmp | sort -nk1


#####################
# OLD STUFF
#####################

### annotation diffential expression
for f in `ls ../[Wer]*uniq.sorted.bam`; do v=${f/..\//};   htseq-count -a 0 -s yes -r pos -f bam -t gene -i ID $f Magnaporthe_oryzae.MG8.25.gff3 > ${v/sorted.bam/count} & done
Rscript ../diff.R WT_1.uniq.count WT_2.uniq.count  WT_3.uniq.count exp5_1.uniq.count exp5_2.uniq.count exp5_3.uniq.count WT_vs_EXP5_csv
Rscript ../diff.R WT_1.uniq.count WT_2.uniq.count  WT_3.uniq.count rbp35_1.uniq.count rbp35_2.uniq.count rbp35_3.uniq.count WT_vs_RBP35_csv

### clusters diffential expression
for f in `ls ../[Wer]*uniq.sorted.bam`; do v=${f/..\//}; v=${v/uniq.sorted.bam/cov}; genomeCoverageBed -ibam $f -g ../genome.txt -d > $v & done
for f in `ls *cov`; do python /media/marco/Elements/m-oryzae-polya/cluster.py $f | sed 's/ /\t/g' > ${f/cov/bed} & done
cat *.bed | sort -k1,1 -k2,2n | bedtools merge -i - | awk '{print $1,"marco","cluster",$2,$3,".",".",".","ID=cluster_"++count"_"}' | sed 's/ /\t/g' > clusters.gff3
for f in `ls ../[Wer]*uniq.sorted.bam`; do v=${f/..\//};   htseq-count -a 0 -s no -r pos -f bam -t cluster -i ID $f clusters.gff3 > ${v/sorted.bam/count} & done
Rscript ../diff.R WT_1.uniq.count WT_2.uniq.count  WT_3.uniq.count rbp35_1.uniq.count rbp35_2.uniq.count rbp35_3.uniq.count WT_vs_RBP35.csv
Rscript ../diff.R WT_1.uniq.count WT_2.uniq.count  WT_3.uniq.count exp5_1.uniq.count exp5_2.uniq.count exp5_3.uniq.count WT_vs_EXP5.csv

## adapters diff expr
for f in `ls ../*fasta.trimmed`; do v=${f/..\//}; grep -v ">" $f | awk '{print substr($0,1,4)}' | sort | uniq -c | awk '{print $2"\t"$1}' > ${v/fasta.trimmed/5prime.count}; done
for f in `ls ../*fasta.trimmed`; do v=${f/..\//}; grep -v ">" $f | awk '{print substr($0,length($0)-3,4)}' | sort | uniq -c | awk '{print $2"\t"$1}' > ${v/fasta.trimmed/3prime.count}; done
Rscript ../diff.R WT_1.5prime.count WT_2.5prime.count WT_3.5prime.count exp5_1.5prime.count exp5_2.5prime.count exp5_3.5prime.count WT_vs_EXP5.5prime.csv
Rscript ../diff.R WT_1.3prime.count WT_2.3prime.count WT_3.3prime.count exp5_1.3prime.count exp5_2.3prime.count exp5_3.3prime.count WT_vs_EXP5.3prime.csv


##### differential coverage
python -c "
window = 3
for line in open('clusters.gff3'):
  items = line.strip().split('\t')
  print items[0] + '\t' + 'dexseq_prepare_annotation.py' + '\t' + 'aggregate_gene' + '\t' + items[3] + '\t' + items[4] + '\t' + \
        items[5] + '\t' +  items[6] + '\t' +  items[7] + '\t' + 'gene_id ' + '\"' + items[8].replace('ID=', '') + '\"'
  left = int(items[3])
  right = int(items[4])
  if right-left > 300: continue
  k = 1
  for i in range(0, right-left, window):
    start = left+i
    if left + i + window - 1 > right:
       end = right
    else:
       end =  left + i + window - 1
    if k < 10: 
      exon = '00' + str(k)
    else:
      exon = '0' + str(k)       
    print items[0] + '\t' + 'dexseq_prepare_annotation.py' + '\t' + 'exonic_part' + '\t' + str(start) + '\t' + str(end) + '\t' + \
        items[5] + '\t' +  items[6] + '\t' +  items[7] + '\t' + ' transcripts \"' + items[8].replace('ID=', '') + '.1\"'+ '; exonic_part_number \"'  + exon + '\";' +  ' gene_id \"' + items[8].replace('ID=', '') + '\"'
    k+=1
"

for f in `ls ../*.uniq.sam`
do
v=${f/..\//}
python -c "
import sys
class Exon:
  def __init__(self, name ,chrx, start, end):
    self.name = name
    self.chrx = chrx
    self.start = start
    self.end = end
    self.count = 0
    
exons = []
for line in open('exons.gtf'):
    items = line.strip().split('\t')
    if items[2] != 'exonic_part': continue
    name = items[8].replace('exon \"', '').replace('\";', '')
    
    exons.append(    Exon(name , items[0], int(items[3]), int(items[4]) )    )


for line in open(sys.argv[1]):
  if line[0] == '@': continue
  items = line.strip().split('\t')
  chrx = items[2]
  if items[1] == '0':
    pos =  int(items[3])
  elif items[1] == '16':
    pos =  int(items[3]) + len(items[9])
  for exon in exons:
    if pos <= exon.end and pos >= exon.start:
       exon.count += 1 
       break
       
for exon in exons:
  print exon.name + '\t' + str(exon.count)
" $f > "_"$v.count &
done

