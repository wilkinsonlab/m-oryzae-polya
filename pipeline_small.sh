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


# aligning not unique
for f in `ls *fastq.trimmed.x`;
do
	bowtie -S -M 1 -v 0 -p 8 --best --strata --al ${f/.trimmed.x/.al} --max  ${f/.trimmed.x/.max} --un ${f/.trimmed.x/.un}  db/bowtie/magna $f | samtools view -Sh -F 4 - >  ${f/.*/.sam}
done

# aligning unique
for f in `ls *fastq.trimmed.x`;
do
	bowtie -S -m 1 -v 0 -p 8 --best --strata db/bowtie/magna $f | samtools view -Sh -F 4 - >  ${f/.*/.uniq.sam}
done


# sort and index
for f in `ls *.sam`; do
	samtools view -bhS $f | samtools sort - ${f/.sam/.sorted}
	samtools index ${f/.sam/.sorted.bam}
done

# create bedgraphs
for f in `ls *sorted.bam`; do genomeCoverageBed -bg -ibam $f -g genome.txt > ${f/.sorted.bam/.bedgraph} & done

# converting aligned reads to fasta
for f in `ls *fastq.al`; do 
	fastq_to_fasta -i $f -o ${f/fastq/fasta}; 
done

### annotation diffential expression
for f in `ls ../[Wer]*uniq.sorted.bam`; do v=${f/..\//};   htseq-count -a 0 -s yes -r pos -f bam -t ncRNA -i ID $f Magnaporthe_oryzae.MG8.25.ncrna.gff3 > ${v/sorted.bam/ncrna.count} & done
Rscript ../diff.R WT_1.uniq.all.count WT_2.uniq.all.count  WT_3.uniq.all.count exp5_1.uniq.all.count exp5_2.uniq.all.count exp5_3.uniq.all.count WT_vs_EXP5_all.csv
Rscript ../diff.R WT_1.uniq.all.count WT_2.uniq.all.count  WT_3.uniq.all.count rbp35_1.uniq.all.count rbp35_2.uniq.all.count rbp35_3.uniq.all.count WT_vs_RBP35_all.csv

## ncRNA  diffential expression (ncrna is db made of unique ncRNA sequences)
for f in `ls ../*[123].fastq.trimmed.x`; do v=${f/..\//}; bowtie -S -p 8 -m 1 -v 0 ncrna --quiet --sam-nohead $f | awk '{if($2==0)print $3}' | sort | uniq -c | awk '{print $2"\t"$1}' > "_"${v/.fastq.trimmed.x/.count}; done
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls _*count`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
cat $tmp > _count
cut -f 1,2 _count > exp5_1.count;cut -f 1,3 _count > exp5_2.count;cut -f 1,4 _count > exp5_3.count; cut -f 1,5 _count > rbp35_1.count; cut -f 1,6 _count > rbp35_2.count; cut -f 1,7 _count > rbp35_3.count; cut -f 1,8 _count > WT_1.count; cut -f 1,9 _count > WT_2.count; cut -f 1,10 _count > WT_3.count
Rscript ../diff.R WT_1.count WT_2.count  WT_3.count rbp35_1.count rbp35_2.count rbp35_3.count WT_vs_RBP35.csv
Rscript ../diff.R WT_1.count WT_2.count  WT_3.count exp5_1.count exp5_2.count exp5_3.count WT_vs_EXP5.csv



### clusters diffential expression
for f in `ls ../[Wer]*uniq.sorted.bam`; do v=${f/..\//}; v=${v/uniq.sorted.bam/cov}; genomeCoverageBed -ibam $f -g ../genome.txt -d > $v & done
for f in `ls *cov`; do python /media/marco/Elements/m-oryzae-polya/cluster.py $f | sed 's/ /\t/g' > ${f/cov/bed} & done
cat *.bed | sort -k1,1 -k2,2n | bedtools merge -i - | awk '{print $1,"marco","cluster",$2,$3,".",".",".","ID=cluster_"++count"_"}' | sed 's/ /\t/g' > cluster.gff3
for f in `ls ../[Wer]*uniq.sorted.bam`; do v=${f/..\//};   htseq-count -a 0 -s no -r pos -f bam -t cluster -i ID $f cluster.gff3 > ${v/sorted.bam/count} & done
Rscript ../diff.R WT_1.uniq.count WT_2.uniq.count  WT_3.uniq.count rbp35_1.uniq.count rbp35_2.uniq.count rbp35_3.uniq.count WT_vs_RBP35.csv
Rscript ../diff.R WT_1.uniq.count WT_2.uniq.count  WT_3.uniq.count exp5_1.uniq.count exp5_2.uniq.count exp5_3.uniq.count WT_vs_EXP5.csv

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
for f in `ls *fa`; do ~/Downloads/trinityrnaseq-2.0.2/Inchworm/bin/inchworm --reads $f --run_inchworm -K 18 -L 18 --num_threads 8 > $f.inch; done
cat *inch | fasta_formatter | fastx_collapser -o all_reads.fa.inch
for f in `ls ../*fastq.trimmed.x`; do v=${f/..\//}; bowtie -S -p 8 -v 0 -m 1 --sam-nohead  all $f  | awk '{if($2!=4)print $3}' | sort | uniq -c | awk '{print $2"\t"$1}' > "_"$v ; done
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
Rscript ../diff.R WT_1.5prime.count WT_2.5prime.count WT_3.5prime.count exp5_1.5prime.count exp5_2.5prime.count exp5_3.5prime.count WT_vs_EXP5.5prime.csv
Rscript ../diff.R WT_1.3prime.count WT_2.3prime.count WT_3.3prime.count exp5_1.3prime.count exp5_2.3prime.count exp5_3.3prime.count WT_vs_EXP5.3prime.csv

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
 
 


# get info, in db folder
# for diff_expr_sequences
rm _*
for f in `ls *[123].fasta.trimmed.x`;
do
bowtie -S -p 8 -v 2 -k 1 db/bowtie/ncrna -f $f --un _un 1> /dev/null 2> _res; grep reported _res >> "_"$f
bowtie -S -p 8 -v 2 -k 1 db/bowtie/retro  -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"$f
bowtie -S -p 8 -v 2 -k 1 db/bowtie/cds -f __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"$f
bowtie -S -p 8 -v 2 -k 1 db/bowtie/utr  -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"$f
bowtie -S -p 8 -v 2 -k 1 db/bowtie/unspliced  -f __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"$f
bowtie -S -p 8 -v 2 -k 1 db/bowtie/magna -f _un --un _unknown 1> /dev/null 2> _res; grep reported _res >> "_"$f
done 

# for cluster, milRNA, assemblies
type="--end-to-end"
for f in *.fa ;
do
rm _*	
bowtie2  -f $f --un=_un -x ../db/bowtie2/ncrna $type -k 1 --quiet --no-head -p 8  | awk '{if($2!=4)print $0}'  > _ncrna_out 
bowtie2  -f _un --un=__un -x ../db/bowtie2/retro $type -k 1 --quiet --no-head -p 8   | awk '{if($2!=4)print $0}'  > _retro_out 
bowtie2  -f __un --un=_un -x ../db/bowtie2/utr $type -k 1 --quiet --no-head -p 8   | awk '{if($2!=4)print $0}' > _utr_out 
bowtie2  -f _un --un=__un -x ../db/bowtie2/cds $type -k 1 --quiet --no-head -p 8   | awk '{if($2!=4)print $0}'  > _cds_out 
bowtie2  -f __un --un=_un -x ../db/bowtie2/unspliced $type -k 1 --quiet --no-head -p 8   | awk '{if($2!=4)print $0}' > _intron_out
bowtie2  -f _un --al _intergenic --un=_unaligned -x ../db/bowtie2/magna $type -k 1 --quiet --no-head -p 8   | awk '{if($2!=4)print $0}'  > _intergenic_out
echo $f 
for g in _ncrna_out _retro_out _cds_out _utr_out _intron_out _intergenic_out; do echo ${g/_out/} > "_"$g; awk '{if($2==0)s="+"; else s="-"; print $1,$3,s}' < $g >> "_"$g ; done
paste  -d "," __ncrna_out __retro_out __cds_out __utr_out __intron_out __intergenic_out
if [ -s _intergenic ] ; then  
	cmsearch --cpu 8 -E 1e-2 --tblout _intergenic_rfam --noali ~/Downloads/Rfam.cm _intergenic > /dev/null
	sed -i -e 's/^#.*//' -e '/^$/d'  _intergenic_rfam; 
	echo "_intergenic_rfam"; cat _intergenic_rfam 
fi
if [ -s _unaligned ] ; then  
	cmsearch --cpu 8 -E 1e-2 --tblout _unaligned_rfam --noali ~/Downloads/Rfam.cm _unaligned > /dev/null
	sed -i -e 's/^#.*//' -e '/^$/d'  _unaligned_rfam;  
	echo "_unaligned_rfam"; cat _unaligned_rfam;
fi
echo ""
done



# for adapters and reads
for f in `ls *.fa`;
do
rm _*	
bowtie -S  -v 0  ../db/bowtie/ncrna  -f $f --un _un  -k 1 --quiet --sam-nohead -p 8  | awk '{if($2!=4)print $0}'  > _ncrna_out 
bowtie -S  -v 0 ../db/bowtie/retro -f _un --un __un   -k 1 --quiet --sam-nohead -p 8   | awk '{if($2!=4)print $0}'  > _retro_out 
bowtie -S  -v 0 ../db/bowtie/utr -f __un --un _un    -k 1 --quiet --sam-nohead -p 8   | awk '{if($2!=4)print $0}' > _utr_out 
bowtie -S  -v 0  ../db/bowtie/cds -f _un --un __un  -k 1 --quiet --sam-nohead -p 8   | awk '{if($2!=4)print $0}'  > _cds_out 
bowtie -S  -v 0 ../db/bowtie/unspliced  -f __un --un _un  -k 1 --quiet --sam-nohead -p 8   | awk '{if($2!=4)print $0}' > _intron_out
bowtie -S -v 0  ../db/bowtie/magna -f _un --al _intergenic --un=_unaligned  -k 1 --quiet --sam-nohead -p 8   | awk '{if($2!=4)print $0}'  > _intergenic_out
echo -ne $f" "  
#for g in _ncrna_out _retro_out _cds_out _utr_out _intron_out _intergenic_out; do echo ${g/_out/} > "_"$g; awk '{if($2==0)s="+"; else s="-"; print $1,$3,s}' < $g >> "_"$g ; done
for g in _ncrna_out _retro_out _cds_out _utr_out _intron_out _intergenic_out; do  wc -l $g | cut -f 1 -d " " >> "_"$g ; done
paste  -d "," __ncrna_out __retro_out __cds_out __utr_out __intron_out __intergenic_out
continue
if [ -s _intergenic ] ; then  
	cmsearch --cpu 8 -E 1e-2 --tblout _intergenic_rfam --noali ~/Downloads/Rfam.cm _intergenic > /dev/null
	sed -i -e 's/^#.*//' -e '/^$/d'  _intergenic_rfam; 
	echo "_intergenic_rfam"; cat _intergenic_rfam 
fi
if [ -s _unaligned ] ; then  
	cmsearch --cpu 8 -E 1e-2 --tblout _unaligned_rfam --noali ~/Downloads/Rfam.cm _unaligned > /dev/null
	sed -i -e 's/^#.*//' -e '/^$/d'  _unaligned_rfam;  
	echo "_unaligned_rfam"; cat _unaligned_rfam;
fi
echo ""

#num=`grep -c ">" $f`
#for g in _ncrna_out _retro_out _cds_out _utr_out _intron_out _intergenic_out; do echo ${g/_out/} > "_"$g; wc -l $g | cut -f 1 -d " " | awk -v num=$num '{print $0/num}' >> "_"$g ; done
#echo -e "adapter\n"$f |  paste  -d "," - __ncrna_out __retro_out __cds_out __utr_out __intron_out __intergenic_out
done


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

