# trimming 3' adaptor, quality and artifacts filtering
for f in `ls *.fastq`; 
do
	fastx_clipper -c  -a TGGAATTCTCGGGTGCCAAGG -i $f | fastq_quality_filter -q 30 -p 70 | fastx_artifacts_filter -o $f".trimmed" &
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
  if len(edit) == 0: continue
  seqs.append(edit)
  if len(seqs) > 100000:
      SeqIO.write(seqs, f, 'fastq')
      seqs = []
SeqIO.write(seqs, f, 'fastq')
" $f $f".x" &
done


# aligning bowtie
for f in `ls *fasta.trimmed.x`;
do
	bowtie -S -M 1 -v 0 -p 8 --best --strata ../db/bowtie/genome $f | samtools view -bSh -F 4 - | samtools sort - visualization/bowtie/${f/.*/.sorted}
done


# convert reads to fasta and collapse 
for f in `ls *fastq.trimmed.x`; do fastq_to_fasta -i $f -o ${f/fastq/fasta} ; done
for f in `ls *fastq.trimmed.x`; do fastx_collapser -i $f -o $f.collapsed ; done

# get reads length
rm __*
for f in  `ls _*coll*ncrna`; do   grep -v ">" $f | awk '{if (length($0)>9) print length($0)}' | sort | uniq -c | awk '{print $2"\t"$1}' > "__"$f  ; done
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls __*`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done ; sort -n $tmp > _tmp
Rscript ../transcriptomic_based/norm.R _tmp; gedit _tmp.norm

# get first nucleotide
rm _*; for f in `ls W*fa` ;  do  grep -v ">" $f | awk '{print substr($0, 1, 1)}' | sort | uniq -c | awk '{print $2"\t"$1}' > "__"$f  ;done
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls __*`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done  
cat $tmp

# get last nucleotide
rm _*; for f in `ls W*fa` ;  do  grep -v ">" $f | awk '{print substr($0, length($0), 1)}' | sort | uniq -c | awk '{print $2"\t"$1}' > "__"$f  ;done
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls __*`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done  
cat $tmp


# get one alignement, more than 1 alignment, no alignments on genome (from genomic_based/clusters)
for f in `ls *sorted.bam`;
do
samtools view $f | cut -f 1,14 | sort | uniq | grep -o NH:i:[0-9]* | sed 's/NH:i://' | awk -v f=$f '{if($0==1) one++; else more++}END{print f" "one" "more}' &
done

### sequences differential expression
for f in `ls *fasta.trimmed.x.collapsed`; do awk '{if (substr($0, 1, 1) == ">"){ split($0, arr, "-"); val=arr[2];} else { print $0"\t"val; val="-" }  }' < $f > $f.count; done
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls *collapsed.count`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
awk '{if ($2>1 || $3>1 || $4>1 || $5>1 || $6>1 || $7>1 || $8>1 || $9>1 || $10>1) print $0 }' <  $tmp > seq.count
cut -f 1,2 seq.count > exp5_1.count;cut -f 1,3 seq.count > exp5_2.count;cut -f 1,4 seq.count > exp5_3.count; cut -f 1,5 seq.count > rbp35_1.count; cut -f 1,6 seq.count > rbp35_2.count; cut -f 1,7 seq.count > rbp35_3.count; cut -f 1,8 seq.count > WT_1.count; cut -f 1,9 seq.count > WT_2.count; cut -f 1,10 seq.count > WT_3.count
Rscript ../../diff.R WT_1.count WT_2.count  WT_3.count exp5_1.count exp5_2.count exp5_3.count WT_vs_EXP5.csv
Rscript ../../diff.R WT_1.count WT_2.count  WT_3.count rbp35_1.count rbp35_2.count rbp35_3.count WT_vs_RBP35.csv
cat WT_vs_EXP5.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3>0  ) print ">"$1"\n"$1}' | sed 's/"//g' > WT_vs_EXP5.up.fa
cat WT_vs_EXP5.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3<0  ) print ">"$1"\n"$1}' | sed 's/"//g' > WT_vs_EXP5.down.fa
cat WT_vs_RBP35.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3>0  ) print ">"$1"\n"$1}' | sed 's/"//g' > WT_vs_RBP35.up.fa
cat WT_vs_RBP35.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3<0  ) print ">"$1"\n"$1}' | sed 's/"//g' > WT_vs_RBP35.down.fa

### HD adapter differential expression
for f in `ls ../*fasta.trimmed`; do v=${f/..\//}; grep -v ">" $f | awk '{if (length($0)>8) print substr($0,1,4)}' | sort | uniq -c | awk '{print $2"\t"$1}' > ${v/fasta.trimmed/5prime.count}; done
for f in `ls ../*fasta.trimmed`; do v=${f/..\//}; grep -v ">" $f | awk '{if (length($0)>8) print substr($0,length($0)-3,4)}' | sort | uniq -c | awk '{print $2"\t"$1}' > ${v/fasta.trimmed/3prime.count}; done
Rscript /media/marco/Elements/m-oryzae-polya/adapters_diff.R WT_1.5prime.count WT_2.5prime.count WT_3.5prime.count exp5_1.5prime.count exp5_2.5prime.count exp5_3.5prime.count WT_vs_EXP5.5prime.csv
Rscript /media/marco/Elements/m-oryzae-polya/adapters_diff.R WT_1.3prime.count WT_2.3prime.count WT_3.3prime.count exp5_1.3prime.count exp5_2.3prime.count exp5_3.3prime.count WT_vs_EXP5.3prime.csv
Rscript /media/marco/Elements/m-oryzae-polya/adapters_diff.R WT_1.5prime.count WT_2.5prime.count WT_3.5prime.count rbp35_1.5prime.count rbp35_2.5prime.count rbp35_3.5prime.count WT_vs_RBP35.5prime.csv
Rscript /media/marco/Elements/m-oryzae-polya/adapters_diff.R WT_1.3prime.count WT_2.3prime.count WT_3.3prime.count rbp35_1.3prime.count rbp35_2.3prime.count rbp35_3.3prime.count WT_vs_RBP35.3prime.csv



# extraction of adapter-related sequences
#sample oriente
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $6!="NA" && $3>0) print $1}' < WT_vs_RBP35.3prime.csv | sed 's/"//g'`; do grep "^"$f ../[Wr]*.fasta.trimmed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_RBP35.3prime_up.fa & 
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $6!="NA" && $3<0) print $1}' < WT_vs_RBP35.3prime.csv | sed 's/"//g'`; do grep "^"$f ../[Wr]*.fasta.trimmed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_RBP35.3prime_down.fa &
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $6!="NA" && $3>0) print $1}' < WT_vs_RBP35.5prime.csv | sed 's/"//g'`; do grep "^"$f ../[Wr]*.fasta.trimmed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_RBP35.5prime_up.fa & 
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $6!="NA" && $3<0) print $1}' < WT_vs_RBP35.5prime.csv | sed 's/"//g'`; do grep "^"$f ../[Wr]*.fasta.trimmed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_RBP35.5prime_down.fa &
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $6!="NA" && $3>0) print $1}' < WT_vs_EXP5.3prime.csv | sed 's/"//g'`; do grep "^"$f ../[We]*.fasta.trimmed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_EXP5.3prime_up.fa & 
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $6!="NA" && $3<0) print $1}' < WT_vs_EXP5.3prime.csv | sed 's/"//g'`; do grep "^"$f ../[We]*.fasta.trimmed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_EXP5.3prime_down.fa &
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $6!="NA" && $3>0) print $1}' < WT_vs_EXP5.5prime.csv | sed 's/"//g'`; do grep "^"$f ../[We]*.fasta.trimmed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_EXP5.5prime_up.fa & 
 for f in `awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $6!="NA" && $3<0) print $1}' < WT_vs_EXP5.5prime.csv | sed 's/"//g'`; do grep "^"$f ../[We]*.fasta.trimmed. -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }'  ; done > WT_vs_EXP5.5prime_down.fa &
#adapter oriented
 for f in `cat *5prime.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $6!="NA" ) print $1}' | sed 's/"//g'`; do  grep "^"$f ../[Wer]*.fasta.trimmed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }' > "_"$f ; done
 for f in `cat *3prime.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $6!="NA" ) print $1}' | sed 's/"//g'`; do  grep $f"$" ../[Wer]*.fasta.trimmed -h | sort -u | awk '{ print ">s_"++count; print substr($1, 5, length($1)-8) }' > "__"$f ; done

### transcripts assembly diffential expression with DESeq2
cat ../../*fasta.trimmed.x > all_reads.fa
fastx_collapser -i all_reads.fa -o all_reads.fa.collapsed
for f in 5 7 9 11 13 15 17 19 21 23 25 27 29 31; do ~/Downloads/trinityrnaseq-2.0.2/Inchworm/bin/inchworm --reads all_reads.fa --run_inchworm -K $f -L $f --min_assembly_coverage 10 --num_threads 7 > all.$f; done
cat all.* | fasta_formatter | fastx_collapser -o all
for f in 5 7 9 11 13 15 17 19 21 23 25 27 29 31; do ~/Downloads/trinityrnaseq-2.0.2/Inchworm/bin/inchworm --reads all --run_inchworm -K $f -L $f --min_assembly_coverage 1 --num_threads 7 > all.inch.$f; done
# choose one
cp all.inch.21 all.inch.final
~/Downloads/segemehl/segemehl.x -x all.inch.final.idx -d all.inch.final
for f in `ls *.fasta.trimmed.x.collapsed`; do  ~/Downloads/segemehl/segemehl.x -A 100 -D 0 -t 8 -m 10 -M 100000 -E 1000 -i assembly/all.inch.final.idx -d assembly/all.inch.final -q $f | samtools view -bS - | samtools sort - assembly/${f/.*/.sorted} ; done
# weight and DE
for f in WT_1  WT_2 WT_3 exp5_1 exp5_2 exp5_3 rbp35_1  rbp35_2 rbp35_3;  do  python /media/marco/Elements/m-oryzae-polya/weight_reads.py ../$f.fasta.trimmed.x.collapsed $f.sorted.bam > $f.weighted ; done
for f in `ls *weighted`; do awk '{arr[$4]+=$7}END{for (k in arr) print k"\t"arr[k]}' < $f > ${f/weighted/expr} ; done
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls *expr`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
awk '{print $1"\t"int($2)}' < $tmp > exp5_1.expr; awk '{print $1"\t"int($3)}' < $tmp > exp5_2.expr; awk '{print $1"\t"int($4)}' < $tmp > exp5_3.expr;
awk '{print $1"\t"int($5)}' < $tmp > rbp35_1.expr; awk '{print $1"\t"int($6)}' < $tmp > rbp35_2.expr; awk '{print $1"\t"int($7)}' < $tmp > rbp35_3.expr;
awk '{print $1"\t"int($8)}' < $tmp > WT_1.expr; awk '{print $1"\t"int($9)}' < $tmp > WT_2.expr; awk '{print $1"\t"int($10)}' < $tmp > WT_3.expr;
Rscript ../../diff.R WT_1.expr WT_2.expr WT_3.expr exp5_1.expr exp5_2.expr exp5_3.expr WT_vs_EXP5.csv
Rscript ../../diff.R WT_1.expr WT_2.expr WT_3.expr rbp35_1.expr rbp35_2.expr rbp35_3.expr WT_vs_RBP35.csv
cat WT_vs_EXP5.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3<0  ) print $1}' | sed -e 's/"//g' > _exp_down &
cat WT_vs_EXP5.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3>0  ) print $1}' | sed -e 's/"//g' > _exp_up  &
cat WT_vs_RBP35.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3<0  ) print $1}' | sed -e 's/"//g' > _rbp_down &
cat WT_vs_RBP35.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3>0  ) print $1}' | sed -e 's/"//g' > _rbp_up
python /media/marco/Elements/m-oryzae-polya/fasta_extract.py all.inch.final _exp_down WT_vs_EXP5.down.fa
python /media/marco/Elements/m-oryzae-polya/fasta_extract.py all.inch.final _exp_up WT_vs_EXP5.up.fa
python /media/marco/Elements/m-oryzae-polya/fasta_extract.py all.inch.final _rbp_down WT_vs_RBP35.down.fa
python /media/marco/Elements/m-oryzae-polya/fasta_extract.py all.inch.final _rbp_up WT_vs_RBP35.up.fa




### for the avr-regions (modify cluster.py)
for f in `ls ../[Wer]*sorted.bam`; do v=${f/..\//}; v=${v/sorted.bam/cov}; echo $v ; genomeCoverageBed -ibam $f -g ../genome.txt -d > $v & done
for f in `ls *cov`; do python /media/marco/Elements/m-oryzae-polya/cluster.py $f  > ${f/cov/bed} & done
cat *.bed | sort -k1,1 -k2,2n | bedtools merge -i - | awk '{print $1,"marco","region",$2,$3,".",".",".","ID=region_"++count"_"}' | sed 's/ /\t/g' > regions.gff3
 
 


# get info unique
rm _*
for f in `ls *.fasta.trimmed.x.collapsed`;
do
db_dir="/media/marco/Elements/EXP5/db/"
~/Downloads/segemehl/segemehl.x -D 0 -A 90 -H 1 -m 10 -M 100000 -E 1000 -t 8 -i $db_dir/segemehl/ncrna.idx -d $db_dir/ncrna.fa -q $f -u _un -nohead > _res; cut -f 1 _res | sort -u | wc -l >  "_"$f; cut -f 1 _res | sort -u | awk '{split($1, arr, "-"); count+=arr[2]}END{print count}'  > "__"$f
~/Downloads/segemehl/segemehl.x -D 0 -A 90 -H 1 -m 10 -M 100000 -E 1000  -t 8 -i $db_dir/segemehl/rrna.idx -d $db_dir/rrna.fa   -q _un -u __un -nohead > _res; cut -f 1 _res | sort -u | wc -l >>   "_"$f; cut -f 1 _res | sort -u | awk '{split($1, arr, "-"); count+=arr[2]}END{print count}'  >> "__"$f
~/Downloads/segemehl/segemehl.x -D 0 -A 90 -H 1 -m 10 -M 100000 -E 1000  -t 8 -i $db_dir/segemehl/retro.idx -d $db_dir/retro.fa   -q __un -u _un -nohead > _res; cut -f 1 _res | sort -u | wc -l >>   "_"$f; cut -f 1 _res | sort -u | awk '{split($1, arr, "-"); count+=arr[2]}END{print count}'  >> "__"$f
~/Downloads/segemehl/segemehl.x -D 0 -A 90 -H 1 -m 10 -M 100000 -E 1000  -t 8 -i $db_dir/segemehl/transcripts.idx -d $db_dir/transcripts.fa   -q _un -u __un -nohead > _res; cut -f 1 _res | sort -u | wc -l >>   "_"$f;cut -f 1 _res | sort -u |  awk '{split($1, arr, "-"); count+=arr[2]}END{print count}'  >> "__"$f
~/Downloads/segemehl/segemehl.x -D 0 -A 90 -H 1 -m 10 -M 100000 -E 1000  -t 8 -i $db_dir/segemehl/unspliced.idx -d $db_dir/unspliced.fa   -q __un -u _un -nohead > _res; cut -f 1 _res | sort -u | wc -l >>  "_"$f; cut -f 1 _res | sort -u | awk '{split($1, arr, "-"); count+=arr[2]}END{print count}'  >> "__"$f
~/Downloads/segemehl/segemehl.x -D 0 -A 90 -H 1 -m 10 -M 100000 -E 1000  -t 8 -i $db_dir/segemehl/genome.idx -d $db_dir/genome.fa  -q _un -u __un -nohead > _res; cut -f 1 _res | sort -u | wc -l >>   "_"$f; cut -f 1 _res | sort -u | awk '{split($1, arr, "-"); count+=arr[2]}END{print count}'  >> "__"$f
~/Downloads/segemehl/segemehl.x -D 0 -A 90 -H 1 -m 10 -M 100000 -E 1000  -t 8 -i $db_dir/segemehl/est.idx -d $db_dir/EST.fa  -q __un -u _unknown -nohead > _res; cut -f 1 _res | sort -u | wc -l >>   "_"$f; cut -f 1 _res | sort -u | awk '{split($1, arr, "-"); count+=arr[2]}END{print count}'  >> "__"$f
done 

rm _*
for f in `ls WT_1*.collapsed`; do 
db_dir="/media/marco/Elements/EXP5/db/"
bowtie -k 1 -p 8  -v 0  $db_dir/bowtie/ncrna -f $f --un _un 1> /dev/null 2> _res; grep reported _res >> "_"$f; 
bowtie -k 1 -p 8  -v 0  $db_dir/bowtie/rrna  -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"$f; 
bowtie -k 1 -p 8  -v 0  $db_dir/bowtie/retro -f __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"$f; 
bowtie -k 1 -p 8  -v 0  $db_dir/bowtie/transcripts  -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"$f; 
bowtie -k 1 -p 8  -v 0  $db_dir/bowtie/unspliced  -f __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"$f; 
bowtie -k 1 -p 8  -v 0  $db_dir/bowtie/genome -f _un --un __un 1> /dev/null 2> _res; grep reported _res >> "_"$f; 
bowtie -k 1 -p 8  -v 0  $db_dir/bowtie/est -f __un --un _un 1> /dev/null 2> _res; grep reported _res >> "_"$f; 
done

# classify siRNA ans single reads
rm _*
for f in `ls *collapsed`;
do
db_dir="/media/marco/Elements/EXP5/db/"
~/Downloads/segemehl/segemehl.x -D 0 -A 100 -m 10 -M 100000 -E 1000 -t 8 -i $db_dir/segemehl/ncrna.idx -d $db_dir/ncrna.fa -q $f -u _un -nohead > _ncrna_out
~/Downloads/segemehl/segemehl.x -D 0 -A 100 -m 10 -M 100000 -E 1000 -t 8 -i $db_dir/segemehl/rrna.idx -d $db_dir/rrna.fa   -q _un -u __un -nohead > _rrna_out
~/Downloads/segemehl/segemehl.x -D 0 -A 100 -m 10 -M 100000 -E 1000 -t 8 -i $db_dir/segemehl/retro.idx -d $db_dir/retro.fa   -q __un -u _un -nohead > _retro_out
~/Downloads/segemehl/segemehl.x -D 0 -A 100 -m 10 -M 100000 -E 1000 -t 8 -i $db_dir/segemehl/transcripts.idx -d $db_dir/transcripts.fa   -q _un -u __un -nohead > _transcripts_out
~/Downloads/segemehl/segemehl.x -D 0 -A 100 -m 10 -M 100000 -E 1000 -t 8 -i $db_dir/segemehl/unspliced.idx -d $db_dir/unspliced.fa   -q __un -u _un -nohead > _introns_out
~/Downloads/segemehl/segemehl.x -D 0 -A 100 -m 10 -M 100000 -E 1000 -t 8 -i $db_dir/segemehl/genome.idx -d $db_dir/genome.fa  -q _un -u __un -nohead > _intergenic_out
~/Downloads/segemehl/segemehl.x -D 0 -A 100 -m 10 -M 100000 -E 1000 -t 8 -i $db_dir/segemehl/est.idx -d $db_dir/EST.fa   -q __un -u _un -nohead > _est_out
#for g in _ncrna_out _rrna_out _retro_out _transcripts_out  _introns_out _intergenic_out _est_out; do cut -f 1 $g | sort -u | wc -l > "_"$g; done
#echo -ne $f"\t" > "_"$f
#paste __ncrna_out __rrna_out __retro_out __transcripts_out    __introns_out __intergenic_out __est_out >> "_"$f
#for g in _ncrna_out _rrna_out _retro_out _transcripts_out _introns_out _intergenic_out _est_out; do cut -f 3 $g | sort | uniq -c > "_"$g; done
#echo $f > "_"$f
#paste __ncrna_out __rrna_out __retro_out __transcripts_out __introns_out __intergenic_out __est_out >> "_"$f
done 



# classify clusters and assemblies
rm _*
for f in `ls WT_vs_EXP5.down.fa`;
do	
db_dir="/media/marco/Elements/EXP5/db/"
# assemblies
awk_filt="\$13<0.0001"
# clusters
#awk_filt="\$4/\$5>0.9||\$4/\$6>0.9"
#awk_filt="\$4==\$5"
touch 	_ncrna_out _rrna_out _retro_out _transcripts_out _gene_out _intergenic_out
blastn -dust no  -num_threads 4 -evalue 100 -task  blastn -query $f -db $db_dir/ncrna.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2> /dev/null| awk '{if ('$awk_filt')  print $0}' > _ncrna_out
blastn -dust no  -num_threads 4 -evalue 100 -task  blastn -query $f -db $db_dir/rrna.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2>  /dev/null | awk '{if ('$awk_filt') print $0}' > _rrna_out
blastn -dust no  -num_threads 4 -evalue 100 -task  blastn -query $f -db $db_dir/retro.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2> /dev/null | awk '{if ('$awk_filt') print $0}' > _retro_out
blastn -dust no  -num_threads 4 -evalue 100 -task  blastn -query $f -db $db_dir/transcripts.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2> /dev/null  | awk '{if ('$awk_filt') print $0}' > _transcripts_out
blastn -dust no  -num_threads 4 -evalue 100 -task  blastn -query $f -db $db_dir/EST.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2> /dev/null  | awk '{if ('$awk_filt') print $0}' > _est_out
blastn -dust no  -num_threads 4 -evalue 100 -task  blastn -query $f -db $db_dir/unspliced.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2> /dev/null |  awk '{if ('$awk_filt') print $0}' > _gene_out
blastn -dust no  -num_threads 4 -evalue 100 -task  blastn -query $f -db $db_dir/genome.fa -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"  -max_target_seqs 1  2> /dev/null  | awk '{if ('$awk_filt') print $0}' > _intergenic_out
for g in  _ncrna_out _rrna_out _retro_out _transcripts_out  _gene_out _intergenic_out _est_out; do sort -k1,1 -k12,12nr -k11,11n $g | sort -u -k1,1 --merge | sort -o $g; done
## print numbers
#for g in _ncrna_out _rrna_out _retro_out _transcripts_out  _gene_out _intergenic_out _est_out; do awk -v g=$g '{print g"\t"$0}' < $g ; done |  awk '{if ($2 in arr == 0)  arr[$2]=$0}END{for (k in arr) print arr[k] }' | cut -f 1 | sort | uniq -c | sed -e 's/_out//' -e 's/_//' | awk '{print $2"\t"$1}'> "__"$f
## print details	
python -c "
yes = []
for file in ('_ncrna_out', '_rrna_out', '_retro_out','_transcripts_out','_gene_out','_intergenic_out','_est_out'):
  out = open('_'+file, 'w')
  for line in open(file):
    items = line.strip().split('\t')
    if items[0] not in yes:
      if file == '_intergenic_out':
        out.write(items[1] + ':' + items[10] + '-' + items[11] + '\n')
      else:  
        out.write(items[1] + '\n')
      yes.append(items[0])
  out.close()       
" 
for k in __ncrna_out __rrna_out __retro_out __transcripts_out  __gene_out __intergenic_out __est_out; do sort $k | uniq -c | sort -o $k; done
paste __ncrna_out __rrna_out __retro_out __transcripts_out  __gene_out __intergenic_out __est_out > "___"$f
done 
## print numbers
#tmp=$(mktemp);tmp2=$(mktemp);for file in `ls __*`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
#cat $tmp


#################### 
# CLUSTERS, BOTH GENOMIC AND TRANSCRIPTOME
#################### 
# get all the alignments
for f in `ls *fasta.trimmed.x.collapsed`;  do  
 #bowtie2 -p 8 -a --end-to-end  -x db/bowtie2/genome -f $f | samtools view -bSh -F 4 - | samtools sort - genomic_based/cluster/bowtie/${f/.*/.sorted};  
 ~/Downloads/segemehl/segemehl.x -D 0 -A 100 -t 8 -m 14 -M 100000  -i ../db/segemehl/genome.idx -d ../db/genome.fa -q $f | samtools view -bSh - | samtools sort - cluster/perfect/${f/.*/.sorted} ;
done

# compute weighted expression of reads
for f in WT rbp35 exp5; do
for g in _1 _2 _3; do   
python /media/marco/Elements/m-oryzae-polya/weight_reads.py ../../$f$g.fasta.trimmed.x.collapsed $f$g.sorted.bam > $f$g.weighted &
done
wait 
done

# for transcripts analysis, we only consider unique alignments
#for f in  *weighted; do  grep "('NH', 1)"  $f > _o; mv _o $f; done

# calculate coverage on weighted reads
for f in WT rbp35 exp5; do
for g in _1 _2 _3; do  
python -c "
import numpy, sys
chrom = {}
for line in open(sys.argv[1], 'r'):
  chrx, length = line.strip().split('\t')
  chrom[chrx] = numpy.array([0.0] * (int(length)+1))
i = 0
for line in open(sys.argv[2], 'r'):
  name, seq, flag, rchrx, rstart, rend, val, cigar, tags = line.strip().split('\t')  
  i += 1
  if i % 1000000 == 0: sys.stderr.write(str(i) + '\n')
  chrom[rchrx][int(rstart):int(rend)+1] += float(val) / (int(rend)+1-int(rstart))
for chrx, array in chrom.items():
  for pos, val in enumerate(array):
     print chrx + '\t' + str(pos) + '\t' + str(val)
" ../../../genome.txt $f$g.weighted > $f$g.cov &
done
wait
done

# detect clusters (change genome.fa for transcripts analysis)
for f in WT rbp35 exp5; do
for g in _1 _2 _3; do  
python /media/marco/Elements/m-oryzae-polya/cluster.py $f$g.cov  | sed 's/ /\t/g' > $f$g.bed &
done
wait
done
cat WT*.bed | sort -k1,1 -k2,2n | bedtools merge -i - -c 5 -o count | awk '{if($4>2)print $0}' > _WT
cat exp*.bed | sort -k1,1 -k2,2n | bedtools merge -i - -c 5 -o count | awk '{if($4>2)print $0}' > _exp
cat rbp*.bed | sort -k1,1 -k2,2n | bedtools merge -i - -c 5 -o count | awk '{if($4>2)print $0}' > _rbp
cat _WT _exp _rbp | sort -k1,1 -k2,2n | bedtools merge -i - | awk '{print $1,"marco","cluster",$2,$3,".",".",".","ID=cluster_"++count"_"}' | sed 's/ /\t/g' > _clusters.gff3
cat clusters.gff3 | bedtools getfasta -bed - -fi ../../../db/genome.fa -fo clusters.fa
awk '{if (substr($0,1,1) == ">") print ">cluster_"++count"_"; else print $0}' clusters.fa > _o ; mv _o clusters.fa

# assign weighted reads to clusters
for t in exp5 rbp35 WT
do
for f in `ls $t*.weighted`;
do
python -c "
import numpy, pysam, sys, math
from Bio import SeqIO

clusters = []
for line in open(sys.argv[1], 'r'):
  items = line.strip().split('\t')
  clusters.append((items[8].replace('ID=', ''),  items[0], int(items[3]), int(items[4])))

### genome
block_span = 100000
chroms = ('supercont8.8', '7', '6', '5', '4', '3', '1', '2')
table_blocks = dict.fromkeys(chroms)
for chrx in chroms:
   table_blocks[chrx] = dict((block, [])
                              for block in range(0, 10000000, block_span))
for chrom in chroms:
    for cluster, chrx, start, end in clusters:
        if chrx == chrom:
          block_start = int(math.floor(start / float(block_span)) * block_span)
          block_end = int(math.floor(end / float(block_span)) * block_span)
          for block in range(block_start, block_end + block_span, block_span):
            table_blocks[chrom][block].append((cluster, chrx, start, end))
###            
### trascriptome
#table_blocks = {}
#for cluster, chrx, start, end in clusters:
#  if not table_blocks.has_key(chrx): table_blocks[chrx] = []
#  table_blocks[chrx].append((cluster, chrx, start, end)) 
###

i = 0
for line in open(sys.argv[2], 'r'):
  name, seq, flag, rchrx, rstart, rend, val, cigar, tags = line.strip().split('\t')
  i += 1
  if i % 1000000 == 0: sys.stderr.write(str(i) + '\n')
### transcriptome
#  if not table_blocks.has_key(rchrx): continue
#  for (cluster, chrx, start, end) in table_blocks[rchrx]:
### 
### genome   
  block = int(math.floor(int(rstart) / float(block_span)) * block_span)
  for (cluster, chrx, start, end) in table_blocks[rchrx][block]:
#  
    if rchrx == chrx and int(rstart) >= start and int(rend) <= end:
       print cluster + '\t' + chrx + '\t' + str(start) + '\t' + str(end) + '\t' + name+ '\t' + seq + '\t' + flag + '\t' + rchrx + '\t' + rstart + '\t' + rend + '\t' + val + '\t' + cigar + '\t' + tags
" clusters.gff3 $f > ${f/weighted/assign} &
done
wait
done


# compute clusters differential expression
for f in `ls *assign`; do awk '{arr[$1]+=$11}END{for (k in arr) print k"\t"arr[k]}' < $f > ${f/assign/expr} ; done
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls *expr`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
awk '{print $1"\t"int($2)}' < $tmp > exp5_1.expr; awk '{print $1"\t"int($3)}' < $tmp > exp5_2.expr; awk '{print $1"\t"int($4)}' < $tmp > exp5_3.expr;
awk '{print $1"\t"int($5)}' < $tmp > rbp35_1.expr; awk '{print $1"\t"int($6)}' < $tmp > rbp35_2.expr; awk '{print $1"\t"int($7)}' < $tmp > rbp35_3.expr;
awk '{print $1"\t"int($8)}' < $tmp > WT_1.expr; awk '{print $1"\t"int($9)}' < $tmp > WT_2.expr; awk '{print $1"\t"int($10)}' < $tmp > WT_3.expr;
Rscript ../../../diff.R WT_1.expr WT_2.expr WT_3.expr exp5_1.expr exp5_2.expr exp5_3.expr WT_vs_EXP5.expr.csv
Rscript ../../../diff.R WT_1.expr WT_2.expr WT_3.expr rbp35_1.expr rbp35_2.expr rbp35_3.expr WT_vs_RBP35.expr.csv
cat WT_vs_EXP5.expr.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3<0  ) print $1}' | sed -e 's/"//g' > _exp_down &
cat WT_vs_EXP5.expr.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3>0  ) print $1}' | sed -e 's/"//g' > _exp_up  &
cat WT_vs_RBP35.expr.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3<0  ) print $1}' | sed -e 's/"//g' > _rbp_down &
cat WT_vs_RBP35.expr.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3>0  ) print $1}' | sed -e 's/"//g' > _rbp_up
python /media/marco/Elements/m-oryzae-polya/fasta_extract.py clusters.fa _exp_down WT_vs_EXP5.expr.down.fa
python /media/marco/Elements/m-oryzae-polya/fasta_extract.py clusters.fa _exp_up WT_vs_EXP5.expr.up.fa
python /media/marco/Elements/m-oryzae-polya/fasta_extract.py clusters.fa _rbp_down WT_vs_RBP35.expr.down.fa
python /media/marco/Elements/m-oryzae-polya/fasta_extract.py clusters.fa _rbp_up WT_vs_RBP35.expr.up.fa


############# 
# siRNA discovery >300 discarted for MEMORY reasons
# only uniquely alignining siRNA
# >= 14 nt
############# 
for f in WT rbp35 exp5; do
for g in _1 _2 _3; do   
count=$(grep $f$g ../../coverage.txt | cut -f 2)
python -c "
import sys,re
class siRNA:
  def __init__(self, cluster):
    self.cluster = cluster
    self.expr = 0.0
    self.seqs = {}
    self.poss = {}
    self.vals = {}
    self.senses = {}
    self.cigars = {}
    self.tagss = {}

siRNAs = {}    
i = 0
for line in open(sys.argv[1], 'r'):
  cluster, chrx, start, end, name, seq, flag, rchrx, rstart, rend, val, cigar, tags = line.strip().split('\t')
  i += 1
  if i % 1000000 == 0: sys.stderr.write(str(i) + '\r')
  if int(end)-int(start) > 300: continue
  #if int(re.findall('\'NH\', ([0-9]*)', tags)[0]) != 1: continue
  if not siRNAs.has_key(cluster): siRNAs[cluster] = siRNA(cluster)
  if flag == '0': 
    siRNAs[cluster].senses[name] = '+'
  elif flag == '16':
    siRNAs[cluster].senses[name] = '-' 
  siRNAs[cluster].expr += float(val)
  siRNAs[cluster].vals[name] = float(val)
  siRNAs[cluster].seqs[name] = seq 
  siRNAs[cluster].poss[name] = (rchrx, rstart, rend)
  siRNAs[cluster].cigars[name] = cigar
  siRNAs[cluster].tagss[name] = tags

  
for cluster, sirna in siRNAs.items():
  for name, val in sirna.vals.items():
    if val / sirna.expr > 0.5 and len(sirna.seqs[name]) >= 17:
      print sirna.cluster + '\t' + name + '\t' + sirna.poss[name][0] + '\t' + sirna.poss[name][1] + '\t' + sirna.poss[name][2] + '\t' + sirna.senses[name] + '\t' + sirna.seqs[name] + '\t' + str(sirna.expr) + '\t' + str(val) + '\t' + sirna.cigars[name] + '\t' + sirna.tagss[name] 
      break
" $f$g.assign | awk -v count=$count '{norm=($9*1000000/count); if (norm > 1) print $0"\t"norm}' > $f$g.siRNA &
done
wait
done

# siRNA confirmed in 2 replicates
sort exp5_[123].siRNA | cut -f 1,3-7 | uniq -c | awk '{if ($1>=2)print $2,$3,$4,$5,$6,$7}'  | sed -e 's/ID=//' -e 's/ /\t/g' > exp5_siRNA
sort rbp35_[123].siRNA | cut -f 1,3-7 | uniq -c | awk '{if ($1>=2)print $2,$3,$4,$5,$6,$7}'  | sed -e 's/ID=//' -e 's/ /\t/g' > rbp35_siRNA
sort WT_[123].siRNA | cut -f 1,3-7 | uniq -c | awk '{if ($1>=2)print $2,$3,$4,$5,$6,$7}'  | sed -e 's/ID=//' -e 's/ /\t/g' > WT_siRNA
sort WT_siRNA rbp35_siRNA exp5_siRNA | uniq > all_siRNA
awk '{print ">"$1"\n"$6}' < exp5_siRNA > exp5_siRNA.fa
awk '{print ">"$1"\n"$6}' < rbp35_siRNA > rbp35_siRNA.fa
awk '{print ">"$1"\n"$6}' < WT_siRNA > WT_siRNA.fa
awk '{print ">"$1"\n"$6}' < all_siRNA > all_siRNA.fa
for f in `cat ../../../transcriptomic_based/single_read/WT_vs_EXP5.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3<0  ) print $1}' | sed -e 's/"//g'`; do grep $'\t'$f$ all_siRNA  ;done | awk '{print ">"$1"\n"$6}' > WT_vs_EXP5.siRNA.down.fa &
for f in `cat ../../../transcriptomic_based/single_read/WT_vs_EXP5.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3>0  ) print $1}' | sed -e 's/"//g'`; do grep $'\t'$f$ all_siRNA  ;done | awk '{print ">"$1"\n"$6}' > WT_vs_EXP5.siRNA.up.fa &
for f in `cat ../../../transcriptomic_based/single_read/WT_vs_RBP35.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3<0  ) print $1}' | sed -e 's/"//g'`; do grep $'\t'$f$ all_siRNA  ;done | awk '{print ">"$1"\n"$6}' > WT_vs_RBP35.siRNA.down.fa &
for f in `cat ../../../transcriptomic_based/single_read/WT_vs_RBP35.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 && $3>0  ) print $1}' | sed -e 's/"//g'`; do grep $'\t'$f$ all_siRNA  ;done | awk '{print ">"$1"\n"$6}' > WT_vs_RBP35.siRNA.up.fa 


# making the table
cut -f 6 all_siRNA | sort -u > _all_seq
awk '{print $6"\t"$1}' all_siRNA > _clusters
cat ../../../transcriptomic_based/single_reads/WT_vs_EXP5.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 ) print $1"\t"$3}' | sed -e 's/"//g'  > _exp_DE &
cat ../../../transcriptomic_based/single_reads/WT_vs_RBP35.csv | awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($7) && $7<0.1 ) print $1"\t"$3}' | sed -e 's/"//g' > _rbp_DE 
for f in `ls [Wer]*.siRNA`; do cut -f 7,10 $f > "_"$f; done 
awk '{print $6"\t"$2,$3,$4,$5}' all_siRNA > _info
python3 -c "
import sys,re
genes = []
for line in open(sys.argv[1]):
  if line[0] == '#': continue
  chrx, x, feature, start, end, x, sense, x, info = line.strip().split('\t')
  if feature in ('exon', 'transcript', 'RNA', 'snoRNA', 'snRNA', 'rRNA'): continue
  genes.append((chrx, x, feature, start, end, x, sense, x, info) )
for line in open(sys.argv[2]):
  cluster, rchrx, rstart, rend, rsense, rseq = line.strip().split('\t')
  sys.stdout.write(rseq + '\t')
  for chrx, x, feature, start, end, x, sense, x, info in genes:
    if rchrx == chrx and int(rstart) >= int(start) and int(rend) <= int(end):
      sys.stdout.write(re.sub('.*[:=]', '', info.split(';')[0]) + ' ' + feature + ',')
  print ()   
"    ../../../db/Magnaporthe_oryzae.MG8.26.gff3 all_siRNA   | sed 's/,$//' > _feature

tmp=$(mktemp);tmp2=$(mktemp);for file in _all_seq _clusters _exp_DE _rbp_DE _WT_1.siRNA _WT_2.siRNA _WT_3.siRNA _exp5_1.siRNA _exp5_2.siRNA _exp5_3.siRNA _rbp35_1.siRNA _rbp35_2.siRNA _rbp35_3.siRNA _info _feature ; do sort -k 1,1 $file -uo $file ;  if [ -s "$tmp" ];   then   join -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file"  > "$tmp2";   else     cp "$file" "$tmp2";   fi;   cp "$tmp2" "$tmp"; done; cat $tmp | sort -nk1 | awk '{if($2!=0) print $0}'

# for visualization purposes
for f in `ls *fasta.trimmed.x`; do bowtie -S -M 1 -v 0 -p 8 --best --strata -f ../db/bowtie/genome $f | samtools view -bSh -F 4 - | samtools sort - visualization/bowtie/${f/.*/.sorted}; done

# convert segmehl alignment for visualization
for f in `ls *sorted.bam`; do
samtools view -h $f | python -c "
import sys
for line in sys.stdin:
  if line[0] == '@':
    print line.strip()
  else:
    #for x in line.strip().split('\t')[11:]:
    #  if x == 'NH:i:1':
        for i in range(int(line.strip().split('\t')[0].split('-')[1])):
          print line.strip()
" | samtools view -bSh - > visualization/$(basename $f) &
done

# R clusters correlation
a=read.table("WT_1.expr", row.names=1)
b=read.table("WT_2.expr", row.names=1)
c=read.table("WT_3.expr", row.names=1)
cor(a,b,method="spearman")
cor(a,c,method="spearman")
cor(b,c,method="spearman")
d=read.table("exp5_1.expr", row.names=1)
e=read.table("exp5_2.expr", row.names=1)
f=read.table("exp5_3.expr", row.names=1)
cor(d,e,method="spearman")
cor(d,f,method="spearman")
cor(e,f,method="spearman")
g=read.table("rbp35_1.expr", row.names=1)
h=read.table("rbp35_2.expr", row.names=1)
i=read.table("rbp35_3.expr", row.names=1)
cor(g,h,method="spearman")
cor(g,i,method="spearman")
cor(h,i,method="spearman")


### expression map

# align first...data.fa
~/Downloads/segemehl/segemehl.x -x data.idx -d data.fa
for f in `ls *.fasta.trimmed.x.collapsed`; do  ~/Downloads/segemehl/segemehl.x -A 100 -D 0 -t 8 -m 10 -M 100000 -E 1000 -i ../maps/data.idx -d ../maps/data.fa -q $f > ../maps/${f/fasta.trimmed.x.collapsed/sam} ; done

# create pos & neg
for f in `ls *sam*`;  do  
grep "SN:.*LN:[0-9]*" -o < $f | sed -e 's/SN://' -e 's/LN://' > "_"$f; 
cat "_"$f | while read a b;  do  
grep $a $f | grep -v "^@" | awk '{if($2==0) {split($1, arr, "-" ); for (i=0;i<arr[2];i++)printf "%d\n", $4/10}}' | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k 1,1 > "_"$f"_"$a"_pos" & 
grep $a $f | grep -v "^@" | awk '{if($2==16) {split($1, arr, "-" ); for (i=0;i<arr[2];i++)printf "%d\n", $4/10}}' | sort | uniq -c | awk '{print $2"\t"(-$1)}' | sort -k 1,1 > "_"$f"_"$a"_neg" & 
done; done

# normalize
for f in exp5_1 exp5_2 exp5_3 rbp35_1 rbp35_2 rbp35_3 WT_1 WT_2 WT_3 ; 
do
  cov=`grep $f coverage.txt | cut -f 2`
  for g in `ls "_"$f*pos`;
  do
    awk -v cov=$cov '{print $1"\t"($2*1000000)/cov}' $g > $g.norm
  done 
  for g in `ls "_"$f*neg`;
  do
    awk -v cov=$cov '{print $1"\t"($2*1000000)/cov}' $g > $g.norm
  done 
done

# plots
for f in `grep ">" data.fa | sed 's/>//'`;
do
echo $f
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls _*$f*pos.norm`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done ; sort -n $tmp > _tmp
Rscript pos.R $f > /dev/null
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls _*$f*neg.norm`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done ; sort -n $tmp > _tmp
Rscript neg.R $f > /dev/null
done

### reads length map

for f in `ls  W*fa`; do  ~/Downloads/segemehl/segemehl.x -A 100 -D 0 -t 8 -m 10 -M 100000 -E 1000 -i ../../maps/data.idx -d ../../maps/data.fa -q $f > ../maps/${f/.fa/.sam} ; done
seq 10 41 > _base
for f in `grep ">" data.fa | sed 's/>//'`;
do
echo $f
echo -e "9\t0" > __a
echo -e "9\t0" > __b
grep $f WT_vs_EXP5.down.sam | grep -v "@" | cut -f 6 | sort | uniq -c | awk '{print $2"\t"$1}'| sed 's/M//' >> __a
grep $f WT_vs_EXP5.up.sam | grep -v "@" | cut -f 6 | sort | uniq -c | awk '{print $2"\t"$1}' | sed 's/M//'>> __b
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls __*`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done ; sort -nk1 $tmp > _tmp
Rscript length.R "EXP5_"$f > /dev/null
echo -e "9\t0" > __a
echo -e "9\t0" > __b
grep $f WT_vs_RBP35.down.sam | grep -v "@" | cut -f 6 | sort | uniq -c | awk '{print $2"\t"$1}'| sed 's/M//' >> __a
grep $f WT_vs_RBP35.up.sam | grep -v "@" | cut -f 6 | sort | uniq -c | awk '{print $2"\t"$1}' | sed 's/M//'>> __b
tmp=$(mktemp);tmp2=$(mktemp);for file in `ls __*`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done ; sort -nk1 $tmp > _tmp
Rscript length.R "RBP35_"$f > /dev/null
done


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

###########################################test snsp

rm *dat
for t in exp5 rbp35 WT
do
for g in `ls ../$t*assign` 
do
python -c "
import sys

i = 0
curr = None
for line in open(sys.argv[1], 'r'):
  cluster, chrx, start, end, name, seq, flag, rchrx, rstart, rend, val = line.strip().split('\t')
  i += 1
  if i % 1000000 == 0: sys.stderr.write(str(i) + '\r')
  if cluster != curr:
    curr = cluster
    fp = open (cluster + sys.argv[1].replace('../', '').replace('.assign','') + '.dat', 'a')
  if int(name.split('-')[1]) == int(float(val)):     
    fp.write(seq + '\t' + val + '\n')
    
" $g &
done
wait
done



for f in `cut -f 9 ../clusters.gff3 | sed 's/ID=//'`; 
do 
  echo -e "-\t0" >> $f"WT_1.dat" 
  echo -e "-\t0" >> $f"WT_2.dat" 
  echo -e "-\t0" >> $f"WT_3.dat" 
  echo -e "-\t0" >> $f"exp5_1.dat" 
  echo -e "-\t0" >> $f"exp5_2.dat" 
  echo -e "-\t0" >> $f"exp5_3.dat" 
  echo -e "-\t0" >> $f"rbp35_1.dat" 
  echo -e "-\t0" >> $f"rbp35_2.dat" 
  echo -e "-\t0" >> $f"rbp35_3.dat" 
done


for f in `cut -f 9 ../clusters.gff3 | sed 's/ID=//'`; 
do 

 sum=`grep $f ../clusters.gff3 | awk  '{if($5-$4 >100 ) print "continue" ;else print "keep"}'`
 if [ "$sum" == "continue" ] ;
 then
   continue     
 fi
  echo $f
 tmp=$(mktemp);tmp2=$(mktemp)	;for file in `ls $f*"_"[123].dat`; do sort -k 1,1 $file -o $file ;    if [ -s "$tmp" ];     then      join  -a 1 -a 2 -e 0 -o auto -t $'\t' "$tmp" "$file" > "$tmp2";     else         cp "$file" "$tmp2";     fi;     cp "$tmp2" "$tmp"; done
 awk '{print $1"\t"int($8)}' < $tmp > $f"_WT_1" ; awk '{print $1"\t"int($9)}' < $tmp > $f"_WT_2" ; awk '{print $1"\t"int($10)}' < $tmp > $f"_WT_3" ;
 awk '{print $1"\t"int($2)}' < $tmp > $f"_exp5_1" ; awk '{print $1"\t"int($3)}' < $tmp > $f"_exp5_2" ; awk '{print $1"\t"int($4)}' < $tmp > $f"_exp5_3" ;
 awk '{print $1"\t"int($5)}' < $tmp > $f"_rbp35_1" ; awk '{print $1"\t"int($6)}' < $tmp > $f"_rbp35_2" ; awk '{print $1"\t"int($7)}' < $tmp > $f"_rbp35_3" ;
done


for f in `cut -f 9 ../clusters.gff3 | sed 's/ID=//'`; 
do 

 sum=`grep $f ../clusters.gff3 | awk  '{if($5-$4 >100 ) print "continue" ;else print "keep"}'`
 if [ "$sum" == "continue" ] ;
 then
   continue     
 fi
  echo $f
 Rscript /media/marco/Elements/m-oryzae-polya/adapters_diff.R $f"_WT_1" $f"_WT_2" $f"_WT_3" $f"_rbp35_1" $f"_rbp35_2" $f"_rbp35_3" $f"WT_vs_RBP35.csv"
done



for f in `cut -f 9 ../clusters.gff3 | sed 's/ID=//'`; 
do 
 echo $f"WT_vs_RBP35.csv"
 awk -F "," 'function isnum(x){return(x==x+0)}  {if(isnum($10) && $10<0.1 ) print $0}' <  $f"WT_vs_RBP35.csv" 
done


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



