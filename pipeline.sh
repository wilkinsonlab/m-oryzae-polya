# trim low quality reads, polyA tails and adapters

# fastx_clipper -v -l 17 -a TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -Q33 -i $f   -o temp  
# fastx_clipper -v -l 17  -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCG      -Q33  -i temp -o "${f%%.*}".trimmed.fastq
    
fastq-mcf -o file1_trimmed_1.fastq -o file2_trimmed_2.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa file1.fastq file2.fastq  

    
# build database

gmap_build -d MG8_21 -D ./MG8_21 $d

# align

gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21  a_1.fastq a_2.fastq > "${f%%.*}".sam

# convert to bam, extract first reads only, sort and index

for f in `ls *.sam`
do
    #samtools view -bS "${f%%.*}".sam > "${f%%.*}".bam
    #samtools view -b -h -f 0x0040 "${f%%.*}".bam > "${f%%.*}"_1.bam
    samtools view -bS -h -f 0x0040 "${f%%.*}".sam > "${f%%.*}"_1.bam
    samtools sort "${f%%.*}"_1.bam "${f%%.*}".sorted
    samtools index "${f%%.*}".sorted.bam
done

# filter out low quality mapping, reads with high A/T content and internal priming

for f in `ls *.sorted.bam`
do
    python ../../m-oryzae-polya/filter.py "${f%%.*}".sorted.bam Magnaporthe_oryzae.MG8.21.dna.toplevel.fa 7 "${f%%.*}".filtered
done

# assign reads to transcripts

for f in `ls *.filtered.bam`
do
    python ../../m-oryzae-polya/assign.py Magnaporthe_oryzae.MG8.21.gff3 "${f%%.*}".filtered.bam "${f%%.*}".assign "${f%%.*}".notassign 7
done

# create bedgraphs  

for f in `ls *.assign`
do
    cat  "${f%%.*}".assign | cut -f 2,3,4,5 | awk '{  if ( $4 == "-" ) print $3,$1  }' | sort -n | uniq -c | awk '{ printf "%s\t%d\t%d\t%d\n",$3,$2-1,$2,$1 }' > "${f%%.*}"_plus.bedgraph 
    cat  "${f%%.*}".assign | cut -f 2,3,4,5 | awk '{  if ( $4 == "+" ) print $2,$1  }' | sort -n | uniq -c | awk '{ printf "%s\t%d\t%d\t%d\n",$3,$2-1,$2,$1 }' > "${f%%.*}"_minus.bedgraph
done

# create not bedgraphs  

for f in `ls *.notassign`
do
    cat  "${f%%.*}".notassign | cut -f 2,3,4,5 | awk '{  if ( $4 == "-" ) print $3,$1  }' | sort -n | uniq -c | awk '{ printf "%s\t%d\t%d\t%d\n",$3,$2-1,$2,$1 }' > "${f%%.*}"_not_plus.bedgraph 
    cat  "${f%%.*}".notassign | cut -f 2,3,4,5 | awk '{  if ( $4 == "+" ) print $2,$1  }' | sort -n | uniq -c | awk '{ printf "%s\t%d\t%d\t%d\n",$3,$2-1,$2,$1 }' > "${f%%.*}"_not_minus.bedgraph
done


# create polyA file

for f in `ls *.assign`
do
    cat "${f%%.*}".assign | cut -f 2,3,4,5,6,7,8 | awk '{  if ( $4 == "+" ) print $2,$1,$4,$5,$6,$7; else print $3,$1,$4,$5,$6,$7  }' | sort -k4 | uniq -c > "${f%%.*}".polyA
done

# create notpolyA file

for f in `ls *.notassign`
do
    cat "${f%%.*}".notassign | cut -f 2,3,4,5,6,7,8 | awk '{  if ( $4 == "+" ) print $2,$1,$4,$5,$6,$7; else print $3,$1,$4,$5,$6,$7  }' | sort -k4 | uniq -c > "${f%%.*}".notpolyA
done

# gene expression count

#cat Magnaporthe_oryzae.MG8.18.gff3 | awk '{if($3 == "gene") print $0}' | grep  "ID=.*;"  -o | sed -e 's/ID=//' -e 's/;//' > _t
cut -f 1 gene_summary.txt  > _t
for f in `ls *.assign`
do
    cut -f 6 $f | cat - _t | sort | uniq -c | awk '{print $2"\t"$1-1}' > "${f%%.*}".expr
done
rm _t 


# extract polyA most significant sites

for f in `ls *.polyA`; do python ../../m-oryzae-polya/polyA_extract.py Magnaporthe_oryzae.MG8.21.gff3 $f "${f%%.*}".expr 33 0.05 5 all > $f"_"all_m; done 

# extract notpolyA most significant sites

for f in `ls *.notpolyA`; do python ../../m-oryzae-polya/polyA_extract_not.py $f 33 100 > $f"_"all_m_high; done 
for f in `ls *.notpolyA`; do python ../../m-oryzae-polya/polyA_extract_not.py $f 33 10 > $f"_"all_m_low; done

# cumulate replicates
for f in  "WT-CM" "WT-MM" "WT--N" "WT--C" "2D4-CM" "2D4-MM" "2D4--N" "2D4--C"
do    
    cat $f"-"1.polyA_all_m $f"-"2.polyA_all_m $f"-"3.polyA_all_m |  sort -k 2,7 | uniq -f 1 -c | awk '{if ($1 >= 2) print 0,$3,$4,$5,$6,$7,$8}' > $f"-"X.polyA_all_m
done
cat WT*X*.polyA_all_m |  sort -k 1,7 | uniq  > WT-ALL-X.polyA_all_m
cat 2D4*X*.polyA_all_m |  sort -k 1,7 | uniq  > 2D4-ALL-X.polyA_all_m


# cumulate not replicates
for f in  "WT-CM" "WT-MM" "WT--N" "WT--C" "2D4-CM" "2D4-MM" "2D4--N" "2D4--C"
do    
    cat $f"-"1.notpolyA_all_m_high $f"-"2.notpolyA_all_m_high $f"-"3.notpolyA_all_m_high |  sort -k 2,7 | uniq -f 1 -c | awk '{if ($1 >= 2) print 0,$3,$4,$5,$6,$7,$8}' > $f"-"X.notpolyA_all_m_high
done
cat WT*X*notpolyA_all_m_high |  sort -k 1,4 | uniq  > WT-ALL-X.notpolyA_all_m_high
cat 2D4*X*notpolyA_all_m_high | sort -k 1,4 | uniq  > 2D4-ALL-X.notpolyA_all_m_high

for f in  "WT-CM" "WT-MM" "WT--N" "WT--C" "2D4-CM" "2D4-MM" "2D4--N" "2D4--C"
do    
    cat $f"-"1.notpolyA_all_m_low $f"-"2.notpolyA_all_m_low $f"-"3.notpolyA_all_m_low |  sort -k 2,7 | uniq -f 1 -c | awk '{if ($1 >= 2) print 0,$3,$4,$5,$6,$7,$8}' > $f"-"X.notpolyA_all_m_low
done
cat WT*X*notpolyA_all_m_low |  sort -k 1,4 | uniq  > WT-ALL-X.notpolyA_all_m_low
cat 2D4*X*notpolyA_all_m_low | sort -k 1,4 | uniq  > 2D4-ALL-X.notpolyA_all_m_low


# gff of polyA
for f in `ls *.polyA_all_m`; do cut $f -d " " -f 1,2,3,4,5 | awk '{ if ($4 == "+") sense = "-"; else sense = "+"; printf "%s\tmarco\tpolyA_site\t%d\t%d\t.\t%s\t.\ttranscript=%s;value=%d\n", $3, $2, $2, sense, $5, $1  }' > "${f%%.*}"_polyA.gff ; done
for f in `ls *.notpolyA_all_m_high`; do cut $f -d " " -f 1,2,3,4,5 | awk '{ if ($4 == "+") sense = "-"; else sense = "+"; printf "%s\tmarco\tpolyA_site\t%d\t%d\t.\t%s\t.\ttranscript=%s;value=%d\n", $3, $2, $2, sense, $5, $1  }' > "${f%%.*}"_notpolyA_high.gff ; done
for f in `ls *.notpolyA_all_m_low`;  do cut $f -d " " -f 1,2,3,4,5 | awk '{ if ($4 == "+") sense = "-"; else sense = "+"; printf "%s\tmarco\tpolyA_site\t%d\t%d\t.\t%s\t.\ttranscript=%s;value=%d\n", $3, $2, $2, sense, $5, $1  }' > "${f%%.*}"_notpolyA_low.gff ; done


# extract single and apa
for f in  `ls *X.polyA_all_m`
do    
    cat $f | sort -k 5 | uniq -u  -f 4 > "${f%%.*}".polyA_sgl_m
    cat $f | sort -k 5 | uniq -D  -f 4 > "${f%%.*}".polyA_apa_m
done 
 

# gene differential expression

Rscript ../../m-oryzae-polya/diff_expr.R WT-CM WT-MM   
Rscript ../../m-oryzae-polya/diff_expr.R WT-CM WT--N    
Rscript ../../m-oryzae-polya/diff_expr.R WT-CM WT--C    
Rscript ../../m-oryzae-polya/diff_expr.R WT-MM WT--N    
Rscript ../../m-oryzae-polya/diff_expr.R WT-MM WT--C   
Rscript ../../m-oryzae-polya/diff_expr.R WT-CM 2D4-CM    
Rscript ../../m-oryzae-polya/diff_expr.R WT-MM 2D4-MM    
Rscript ../../m-oryzae-polya/diff_expr.R WT--N 2D4--N    
Rscript ../../m-oryzae-polya/diff_expr.R WT--C 2D4--C    
Rscript ../../m-oryzae-polya/diff_expr.R 2D4-CM 2D4-MM    
Rscript ../../m-oryzae-polya/diff_expr.R 2D4-CM 2D4--N     
Rscript ../../m-oryzae-polya/diff_expr.R 2D4-CM 2D4--C    
Rscript ../../m-oryzae-polya/diff_expr.R 2D4-MM 2D4--N    
Rscript ../../m-oryzae-polya/diff_expr.R 2D4-MM 2D4--C    
# version 18
#for f in diff_expr/*expr.csv
#do
# cat $f | awk -F "," '{if ($8 < 0.05 && $6 > 0) print $0}' > "${f%%.*}"_up.csv
# cat $f | awk -F "," '{if ($8 < 0.05 && $6 < 0) print $0}' > "${f%%.*}"_down.csv
#done
# version 21
for f in diff_expr/*expr.csv
do
 cat $f | awk -F "," '{if ($7 < 0.05 && $3 > 0) print $0}' > "${f%%.*}"_up.csv
 cat $f | awk -F "," '{if ($7 < 0.05 && $3 < 0) print $0}' > "${f%%.*}"_down.csv
done



# OLD differential polyA (p-value < 0.1)

function old_diff {
 s1=$1
 s2=$2
 c1=$3
 c2=$4
 cat $s1-$c1-X.polyA_all_m $s2-$c2-X.polyA_all_m | sort -k 1,7 | uniq > _k
 cat _k $s1-$c1-1.polyA | sed s'/^ *//' | sort -k 5,5 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3":"$2":"$4":"$5":"$6":"$7"\t"$1}' | sort -k 1,1 > _a
 cat _k $s1-$c1-2.polyA | sed s'/^ *//' | sort -k 5,5 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3":"$2":"$4":"$5":"$6":"$7"\t"$1}' | sort -k 1,1 > _b
 cat _k $s1-$c1-3.polyA | sed s'/^ *//' | sort -k 5,5 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3":"$2":"$4":"$5":"$6":"$7"\t"$1}' | sort -k 1,1 > _c
 cat _k $s2-$c2-1.polyA | sed s'/^ *//' | sort -k 5,5 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3":"$2":"$4":"$5":"$6":"$7"\t"$1}' | sort -k 1,1 > _d
 cat _k $s2-$c2-2.polyA | sed s'/^ *//' | sort -k 5,5 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3":"$2":"$4":"$5":"$6":"$7"\t"$1}' | sort -k 1,1 > _e
 cat _k $s2-$c2-3.polyA | sed s'/^ *//' | sort -k 5,5 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3":"$2":"$4":"$5":"$6":"$7"\t"$1}' | sort -k 1,1 > _f
 join  _a _b -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,2.2 -t $'\t' > _t1
 join  _t1 _c -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,2.2 -t $'\t'  > _t2
 join  _t2 _d -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,1.4,2.2 -t $'\t'  > _t3
 join  _t3 _e -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,1.4,1.5,2.2 -t $'\t'  > _t4
 join  _t4 _f -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,1.4,1.5,1.6,2.2 -t $'\t'  > _t5
 cat _t5 | awk -F "\t" '{print "\""$1"\"""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.count
 python ../../m-oryzae-polya/norm.py diff_expr/$s1"-"$c1"_"vs"_"$s2"-"$c2"_"expr.count $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.count > _n
 Rscript ../../m-oryzae-polya/diff_polyA.R _n $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.csv
 cat $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.csv | awk -F "," '{if ($8 < 0.1 && $6 > 0) print $1}' | awk -F ":" '{print 0,$2,$1,$3,$4,$5,$6}' > $s1"-"$c1"_"vs"_"$s2"-"$c2"_"up.polyA_all_m
 cat $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.csv | awk -F "," '{if ($8 < 0.1 && $6 < 0) print $1}' | awk -F ":" '{print 0,$2,$1,$3,$4,$5,$6}' > $s1"-"$c1"_"vs"_"$s2"-"$c2"_"down.polyA_all_m
 mv $s1"-"$c1"_"vs"_"$s2"-"$c2"_"up.polyA_all_m $s1"-"$c1"_"vs"_"$s2"-"$c2"_"down.polyA_all_m $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.count $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.csv diff_polyA
 rm _*
}

function calc_diff {
 s1=$1
 s2=$2
 c1=$3
 c2=$4
 cat $s1-$c1-X.polyA_all_m $s2-$c2-X.polyA_all_m | sort -k 1,7 | uniq > _k
 cat _k $s1-$c1-1.polyA | sed s'/^ *//' | sort -k 5,5 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $5":"$3"@"$2"@"$4"@"$6"@"$7"\t"$1}' | sort -k 1,1 > _a
 cat _k $s1-$c1-2.polyA | sed s'/^ *//' | sort -k 5,5 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $5":"$3"@"$2"@"$4"@"$6"@"$7"\t"$1}' | sort -k 1,1 > _b
 cat _k $s1-$c1-3.polyA | sed s'/^ *//' | sort -k 5,5 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $5":"$3"@"$2"@"$4"@"$6"@"$7"\t"$1}' | sort -k 1,1 > _c
 cat _k $s2-$c2-1.polyA | sed s'/^ *//' | sort -k 5,5 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $5":"$3"@"$2"@"$4"@"$6"@"$7"\t"$1}' | sort -k 1,1 > _d
 cat _k $s2-$c2-2.polyA | sed s'/^ *//' | sort -k 5,5 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $5":"$3"@"$2"@"$4"@"$6"@"$7"\t"$1}' | sort -k 1,1 > _e
 cat _k $s2-$c2-3.polyA | sed s'/^ *//' | sort -k 5,5 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $5":"$3"@"$2"@"$4"@"$6"@"$7"\t"$1}' | sort -k 1,1 > _f
 join  _a _b -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,2.2 -t $'\t' > _t1
 join  _t1 _c -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,2.2 -t $'\t'  > _t2
 join  _t2 _d -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,1.4,2.2 -t $'\t'  > _t3
 join  _t3 _e -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,1.4,1.5,2.2 -t $'\t'  > _t4
 join  _t4 _f -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,1.4,1.5,1.6,2.2 -t $'\t'  > _t5
 cat _t5 | awk '{print $1"\t"$2}' > _out1
 cat _t5 | awk '{print $1"\t"$3}' > _out2 
 cat _t5 | awk '{print $1"\t"$4}' > _out3
 cat _t5 | awk '{print $1"\t"$5}' > _out4
 cat _t5 | awk '{print $1"\t"$6}' > _out5
 cat _t5 | awk '{print $1"\t"$7}' > _out6
 mv _t5 $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.count
 Rscript ../../m-oryzae-polya/diff_polyA.R $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA
 # version 18
 #cat $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.csv | awk -F "," '{if ($5 < 0.05 && $7 < 0) print $1,$2}' | sed -e 's/EChromosome/Chromosome/' -e 's/@/ /g'  | awk '{print 0,$3,$2,$4,$1,$5,$6}' | sed 's/  / /' | sort -k 1,7 > $s1"-"$c1"_"vs"_"$s2"-"$c2"_"down.polyA_all_m
 #cat $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.csv | awk -F "," '{if ($5 < 0.05 && $7 > 0) print $1,$2}' | sed -e 's/EChromosome/Chromosome/' -e 's/@/ /g'  | awk '{print 0,$3,$2,$4,$1,$5,$6}' | sed 's/  / /' | sort -k 1,7 > $s1"-"$c1"_"vs"_"$s2"-"$c2"_"up.polyA_all_m
 # version 21
 cat $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.csv | awk -F "," '{if ($7 < 0.05 && $10 < 0) print $1,$2}' | sed -e 's/@/ /g' -e 's/ E\(.*\)/ \1/'  | awk '{print 0,$3,$2,$4,$1,$5,$6}' | sed 's/  / /' | sort -k 1,7 > $s1"-"$c1"_"vs"_"$s2"-"$c2"_"down.polyA_all_m
 cat $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.csv | awk -F "," '{if ($7 < 0.05 && $10 > 0) print $1,$2}' | sed -e 's/@/ /g' -e 's/ E\(.*\)/ \1/'  | awk '{print 0,$3,$2,$4,$1,$5,$6}' | sed 's/  / /' | sort -k 1,7 > $s1"-"$c1"_"vs"_"$s2"-"$c2"_"up.polyA_all_m
 mv $s1"-"$c1"_"vs"_"$s2"-"$c2"_"up.polyA_all_m $s1"-"$c1"_"vs"_"$s2"-"$c2"_"down.polyA_all_m $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.count $s1"-"$c1"_"vs"_"$s2"-"$c2"_"polyA.csv diff_polyA
 rm _*
}  

calc_diff "WT" "2D4" "CM" "CM"
calc_diff "WT" "2D4" "MM" "MM"
calc_diff "WT" "2D4" "-N" "-N"
calc_diff "WT" "2D4" "-C" "-C"
calc_diff "WT" "WT" "CM" "MM"
calc_diff "WT" "WT" "CM" "-N"
calc_diff "WT" "WT" "CM" "-C"
calc_diff "WT" "WT" "MM" "-N"
calc_diff "WT" "WT" "MM" "-C"
calc_diff "2D4" "2D4" "CM" "MM"
calc_diff "2D4" "2D4" "CM" "-N"
calc_diff "2D4" "2D4" "CM" "-C"
calc_diff "2D4" "2D4" "MM" "-N"
calc_diff "2D4" "2D4" "MM" "-C"

# differential notpolyA (old)

function old_notdiff {
 s1=$1
 s2=$2
 c1=$3
 c2=$4
 
 cat $s1-$c1-X.notpolyA_all_m $s2-$c2-X.notpolyA_all_m | sort -k 1,4 | uniq > _k
 cat _k $s1-$c1-1.notpolyA | sed s'/^ *//' | sort -k 3,3 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3":"$2":"$4"\t"$1}' | sort -k 1,1 > _a
 cat _k $s1-$c1-2.notpolyA | sed s'/^ *//' | sort -k 3,3 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3":"$2":"$4"\t"$1}' | sort -k 1,1 > _b
 cat _k $s1-$c1-3.notpolyA | sed s'/^ *//' | sort -k 3,3 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3":"$2":"$4"\t"$1}' | sort -k 1,1 > _c
 cat _k $s2-$c2-1.notpolyA | sed s'/^ *//' | sort -k 3,3 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3":"$2":"$4"\t"$1}' | sort -k 1,1  > _d
 cat _k $s2-$c2-2.notpolyA | sed s'/^ *//' | sort -k 3,3 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3":"$2":"$4"\t"$1}' | sort -k 1,1  > _e
 cat _k $s2-$c2-3.notpolyA | sed s'/^ *//' | sort -k 3,3 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3":"$2":"$4"\t"$1}' | sort -k 1,1  > _f
 join  _a _b -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,2.2 -t $'\t' > _t1
 join  _t1 _c -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,2.2 -t $'\t'  > _t2
 join  _t2 _d -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,1.4,2.2 -t $'\t'  > _t3
 join  _t3 _e -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,1.4,1.5,2.2 -t $'\t'  > _t4
 join  _t4 _f -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,1.4,1.5,1.6,2.2 -t $'\t'  > _t5
 cat _t5 | awk -F "\t" '{print "\""$1"\"""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > $s1"-"$c1"_"vs"_"$s2"-"$c2"_"notpolyA.count
 Rscript ../../m-oryzae-polya/diff_polyA.R $s1"-"$c1"_"vs"_"$s2"-"$c2"_"notpolyA.count $s1"-"$c1"_"vs"_"$s2"-"$c2"_"notpolyA.csv
 mv  $s1"-"$c1"_"vs"_"$s2"-"$c2"_"notpolyA.count $s1"-"$c1"_"vs"_"$s2"-"$c2"_"notpolyA.csv notdiff_polyA
 rm _*
}  

function calc_notdiff {
 s1=$1
 s2=$2
 c1=$3
 c2=$4
 type="low"
 cat $s1-$c1-X.notpolyA_all_m_$type $s2-$c2-X.notpolyA_all_m_$type | sort -k 1,4 | uniq > _k
 cat _k $s1-$c1-1.notpolyA | sed s'/^ *//' | sort -k 3,3 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3"@"$2"@"$4"\t"$1}' | sort -k 1,1 > _a
 cat _k $s1-$c1-2.notpolyA | sed s'/^ *//' | sort -k 3,3 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3"@"$2"@"$4"\t"$1}' | sort -k 1,1 > _b
 cat _k $s1-$c1-3.notpolyA | sed s'/^ *//' | sort -k 3,3 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3"@"$2"@"$4"\t"$1}' | sort -k 1,1 > _c
 cat _k $s2-$c2-1.notpolyA | sed s'/^ *//' | sort -k 3,3 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3"@"$2"@"$4"\t"$1}' | sort -k 1,1 > _d
 cat _k $s2-$c2-2.notpolyA | sed s'/^ *//' | sort -k 3,3 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3"@"$2"@"$4"\t"$1}' | sort -k 1,1 > _e
 cat _k $s2-$c2-3.notpolyA | sed s'/^ *//' | sort -k 3,3 -k 2 | uniq -f 1 -D | awk '{if ($1 > 0) print $3"@"$2"@"$4"\t"$1}' | sort -k 1,1 > _f
 join  _a _b -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,2.2 -t $'\t' > _t1
 join  _t1 _c -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,2.2 -t $'\t'  > _t2
 join  _t2 _d -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,1.4,2.2 -t $'\t'  > _t3
 join  _t3 _e -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,1.4,1.5,2.2 -t $'\t'  > _t4
 join  _t4 _f -1 1 -2 1 -a 1 -a 2 -e 0  -o 0,1.2,1.3,1.4,1.5,1.6,2.2 -t $'\t'  > _t5
 cat _t5 | awk '{print $1"\t"$2}' > _out1
 cat _t5 | awk '{print $1"\t"$3}' > _out2 
 cat _t5 | awk '{print $1"\t"$4}' > _out3
 cat _t5 | awk '{print $1"\t"$5}' > _out4
 cat _t5 | awk '{print $1"\t"$6}' > _out5
 cat _t5 | awk '{print $1"\t"$7}' > _out6
 mv _t5 $s1"-"$c1"_"vs"_"$s2"-"$c2"_"notpolyA.count
 Rscript ../../m-oryzae-polya/notdiff_polyA.R  $s1"-"$c1"_"vs"_"$s2"-"$c2"_"notpolyA".csv"
 mv  $s1"-"$c1"_"vs"_"$s2"-"$c2"_"notpolyA.count $s1"-"$c1"_"vs"_"$s2"-"$c2"_"notpolyA.csv notdiff_polyA_$type
 rm _*
}  

calc_notdiff "WT" "2D4" "CM" "CM"
calc_notdiff "WT" "2D4" "MM" "MM"
calc_notdiff "WT" "2D4" "-N" "-N"
calc_notdiff "WT" "2D4" "-C" "-C"
calc_notdiff "WT" "WT" "CM" "MM"
calc_notdiff "WT" "WT" "CM" "-N"
calc_notdiff "WT" "WT" "CM" "-C"
calc_notdiff "WT" "WT" "MM" "-N"
calc_notdiff "WT" "WT" "MM" "-C"
calc_notdiff "2D4" "2D4" "CM" "MM"
calc_notdiff "2D4" "2D4" "CM" "-N"
calc_notdiff "2D4" "2D4" "CM" "-C"
calc_notdiff "2D4" "2D4" "MM" "-N"
calc_notdiff "2D4" "2D4" "MM" "-C"



# kegg enrichment				
function kegg_enrich {
	type=$1
	file_p=$2
    file_t=$3
	echo $file_p
	echo "term"$'\t'"a"$'\t'"b"$'\t'"c"$'\t'"d"$'\t'"p_value"$'\t'"description"$'\t'"genes" > "${file_p%%.*}"_kegg_enrich.tsv
    grep -f $file_p $type | awk '{if ($2 != "") print $2}' | sort | uniq > _de_list
    awk '{if ($2 != "") print $2}' < $type | sort | uniq > _kegg
	cat $file_p $file_p $file_t | sort | uniq -u | cat - _kegg | sort | uniq -d > _nde_list
	cat _de_list | xargs -ipat grep pat $type | cut -f 1  | sort | uniq >  _kegg_list
	for kegg in `cat _kegg_list`; do grep $kegg $type | cut -f 2 | grep MGG | sort | uniq > "_"$kegg; done
	for f in `ls _mgr*`
	do
	    a=$(cat _de_list $f | sort | uniq -d | wc -l )
	    b=$(cat _de_list $f $f | sort | uniq -u | wc -l)
	    c=$(cat _nde_list $f | sort | uniq -d | wc -l)
	    d=$(cat _nde_list $f $f | sort | uniq -u | wc -l)
	    res=$(Rscript ../../m-oryzae-polya/fisher_test.R $a $b $c $d )
	    echo -ne ${f/_/}$'\t'"$res"$'\t' >> "${file_p%%.*}"_kegg_enrich.tsv
	    grep -m 1 "${f/_/}" $type | cut -f 3 | tr '\n' '\t' >> "${file_p%%.*}"_kegg_enrich.tsv
	    cat _de_list $f | sort | uniq -d | tr '\n' ',' | sed 's/,$/\n/' >>  "${file_p%%.*}"_kegg_enrich.tsv
	done
    Rscript ../../m-oryzae-polya/FDR.R "${file_p%%.*}"_kegg_enrich.tsv
	rm _kegg _de_list _nde_list _kegg_list _mgr*
}




# reactome enrichment (s.cerevisiae only)
function reactome_enrich {
	type=$1
	file_p=$2
        file_t=$3
	echo $file_p
	echo "term"$'\t'"a"$'\t'"b"$'\t'"c"$'\t'"d"$'\t'"p_value"$'\t'"description"$'\t'"genes" > "${file_p%%.*}"_reactome_enrich.tsv
    grep -f $file_p $type | awk '{if ($2 != "") print $2}' | sort | uniq > _de_list
    awk '{if ($2 != "") print $2}' < $type | sort | uniq > _reactome
	cat $file_p $file_p $file_t | sort | uniq -u | cat - _reactome | sort | uniq -d > _nde_list
	cat _de_list | xargs -ipat grep pat $type | cut -f 1  | sort | uniq >  _reactome_list
	for reactome in `cat _reactome_list`; do grep $reactome $type | cut -f 2 | sort | uniq > "_"$reactome; done
	for f in `ls _REACT*`
	do
	    a=$(cat _de_list $f | sort | uniq -d | wc -l )
	    b=$(cat _de_list $f $f | sort | uniq -u | wc -l)
	    c=$(cat _nde_list $f | sort | uniq -d | wc -l)
	    d=$(cat _nde_list $f $f | sort | uniq -u | wc -l)
	    res=$(Rscript ../../../m-oryzae-polya/fisher_test.R $a $b $c $d )
	    echo -ne ${f/_/}$'\t'"$res"$'\t' >> "${file_p%%.*}"_reactome_enrich.tsv
	    grep -m 1 "${f/_/}" $type | cut -f 3 | tr '\n' '\t' >> "${file_p%%.*}"_reactome_enrich.tsv
	    cat _de_list $f | sort | uniq -d | tr '\n' ',' | sed 's/,$/\n/' >>  "${file_p%%.*}"_reactome_enrich.tsv
	done
    Rscript ../../../m-oryzae-polya/FDR.R "${file_p%%.*}"_reactome_enrich.tsv
	rm _reactome _de_list _nde_list _reactome_list _REACT*
}

# to crate kegg graphs
for f in `ls *up_kegg_enrich.tsv`; do n=`echo $f | sed 's/expr_up_kegg_enrich.tsv/plot/'` ; awk  -F "\t" '{if($7<0.05)print $1"\t"$2"\t"$7"\t"$8"\tUP"}' < $f > _$n; done
for f in `ls *down_kegg_enrich.tsv`; do n=`echo $f | sed 's/expr_down_kegg_enrich.tsv/plot/'` ; awk  -F "\t" '{if($7<0.05)print $1"\t"$2"\t"$7"\t"$8"\tDOWN"}' < $f >> _$n; done


# extract and sort data above for xls
for f in `ls *up_kegg_enrich.tsv`; do echo $f ; awk -v n=$(sed -e 's/^_/diff_expr\//' -e 's/_kegg_enrich\.tsv/\.csv/'  <<< $f)  -F "\t" '{if($7<0.05){ print $1"\t"$7"\t"$8; split($9,k,","); for (a in k){ system(" grep "k[a]" "n" | cut -f 1,3 -d \",\" | tr \",\" \"\t\"  | tr \"\\n\" \"\\t\"   "); system(" grep "k[a]" _info | cut -f 2,3,4")  } }} ' < $f; done > _up
python -c "
import operator
curr = []
for line in open('_up', 'r'):
  if line[0] == '_' or line[0:3] == 'mgr': 
     for x in sorted(curr, key = operator.itemgetter(1), reverse=True)[0:5]:
        x[1] = str(x[1])
        print '\t'.join(x) 
     curr = [] 
     print
     print line
  else:
     items = line.strip().split('\t')
     items[1] = float(items[1])
     curr.append(items)
for x in sorted(curr, key = operator.itemgetter(1))[0:5]:
  x[1] = str(x[1])
  print '\t'.join(x) 
"


# calculate GO enrichments with gprofile
for f in `ls *up* `; do Rscript ../../../m-oryzae-polya/gprofiler.R moryzae $f "${f/up/back}" T; done
for f in `ls *down* `; do Rscript ../../../m-oryzae-polya/gprofiler.R moryzae $f "${f/down/back}" T; done
for f in `ls *short* `; do Rscript ../../../m-oryzae-polya/gprofiler.R moryzae $f "${f/short/back}" F; done
for f in `ls *long* `; do Rscript ../../../m-oryzae-polya/gprofiler.R moryzae $f "${f/long/back}" F; done
# yeast orthologs
for f in `ls *up* `; do Rscript ../../../../m-oryzae-polya/gprofiler.R scerevisiae $f "${f/up/back}" T; done
for f in `ls *down* `; do Rscript ../../../../m-oryzae-polya/gprofiler.R scerevisiae $f "${f/down/back}" T; done
for f in `ls *short* `; do Rscript ../../../../m-oryzae-polya/gprofiler.R scerevisiae $f "${f/short/back}" F; done
for f in `ls *long* `; do Rscript ../../../../m-oryzae-polya/gprofiler.R scerevisiae $f "${f/long/back}" F; done


# shuffle fasta
python -c "
import random
i = 0
for line in open('_all_s', 'r'):
 if line[0] == '>': continue
 l = list(line.strip())
 random.shuffle(l)
 i += 1
 print '>', i 
 print ''.join(l) 
"


# glam alignment
glam2 n WT-CM-X_ARICH_sgl_m.fam -n 40000 -w 6 -O glam2_ARICH_sgl
fimo -oc fimo_ARICH_sgl --norc --thresh 1e-2 glam2_ARICH_sgl/glam2.meme WT-CM-X_TOT_sgl_m.fam 

cut -f 9 fimo_out/fimo.txt | sort | uniq > _t
for i in `cat _t`; do echo -ne $i" "; grep $i Rozella_allomycis.introns.3prime.fa -c | awk -v num=`grep -c ">" Rozella_allomycis.introns.3prime.fa` '{print $1/num*100}'; done	
	
# RNA structure base probabilities
file=_k.fa
mkdir "_"$file".out"
cd "_"$file".out"
RNAfold -p -d2 --noLP < ../$file > /dev/null
for f in `ls *_dp.ps`; do grep -E "[0-9].[ul]box" $f | awk '{if (arr[$1] == 0) arr[$1] = $3 ; else {if ($3 > arr[$1]) arr[$1]=$3}; if (arr[$2] == 0) arr[$2] = $3; else {if ($3>arr[$2]) arr[$2]=$3}; count[$1]++; count[$2]++} END {for (k in arr) print k, arr[k]}' | sort -n > "_"$f; done
cat _*_dp.ps | awk '{arr[$1]+=$2; count[$1]++;}END{for (k in arr) print k, arr[k]/count[k]}' | sort -n | cut -f 2 -d " "
cd ..
rm -rf "_"$file".out"

	
# motif scan 
python -c "
import re, sys
from Bio import SeqIO
f = open(sys.argv[1], 'r')
s = sys.argv[2]
s = s.replace('R', '[GA]').replace('Y', '[TC]').replace('S', '[GC]').replace('W', '[TA]').replace('K', '[GT]').replace('M', '[AC]').replace('D', '[GTA]').replace('H', '[TAC]').replace('B', '[GTC]').replace('V', '[GAC]').replace('N', '[ATGC]')
print s
d = [0 for x in range(int(sys.argv[3]))]
c = 0.0
for record in SeqIO.parse(f, 'fasta'):
  for m in re.finditer(s, str(record.seq)):
    d[m.start(0)] += 1
  c += 1		
for v in d:
  print v / c * 100		
" 




# extract 3'UTR sequences
python -c "
import os,sys,re
gff_file = open(sys.argv[1], 'r')
polyA_file = open(sys.argv[2], 'r')

stops = {}
for line in gff_file:
    if line[0] == '#' or line[0] == '\n':
        continue
    items = line.split('\t')
    if items[2] == 'three_prime_UTR':
        chrx = items[0]
        for x in items[8].split(';'):
            if x.split('=')[0] == 'Parent':
                name = x.split('=')[1].strip()
        name = re.sub(r'T.', '', name)
        sense = items[6]
        if sense == '+':
            pos = int(items[3])
        else:
            pos = int(items[4])
        stops[name] = pos

for line in polyA_file:
    items = line.strip().split(' ')
    pos = int(items[1])
    chrx = items[2]
    sense = items[3]
    gene = items[4]
    if not stops.has_key(gene): continue
    if sense == '+':
        sense = '-'
        sense_ = '2'
        start = pos
        end = stops[gene]
    else:
        sense = '+'
        sense_ = '1'
        start = stops[gene]
        end = pos
    if start > end: continue
    os.system('fastacmd -s \"lcl|' + chrx + '\" -S ' + sense_ + ' -L ' +  str(start) + ',' + str(end) + ' -d Magnaporthe_oryzae.MG8.21.dna.toplevel.fa | sed -e \'s/ dna.*/@' + gene  + '@' + sense + '/\' -e \'s/lcl|//\'  ')

gff_file.close()
polyA_file.close()
" Magnaporthe_oryzae.MG8.21.gff3 2D4-CM-X.polyA_all_m > _2D4.fa


# extract 3'UTR or intra-APA sequences (for miRNA search)
grep three_prime_UTR Magnaporthe_oryzae.MG8.21.gff3 | sed  -e 's/Parent=//' -e 's/T[0-9]*;//' | cut -f 4,9| awk '{print $2,$1}' | sort > _g
sort -k 5,5 -k 2 WT-CM-X.polyA_all_m | awk '{arr[$5"@"$3"@"$4]=arr[$5"@"$3"@"$4]$2":"} END {for(x in arr) print x,arr[x]}  ' | sed -e 's/@/ /g' -e 's/+/2/' -e 's/-/1/'  -e 's/:$//' | sort > _t
join _g _t | awk '{split($5, arr, ":"); for (x in arr) if (($4 == 1 && arr[x] > $2  && arr[x+1] > $2 && arr[x+1] != "") || ($4 == 2 && arr[x] < $2  && arr[x+1] < $2 && arr[x+1] != "") ) system("echo -n "$1"; fastacmd -d Magnaporthe_oryzae.MG8.21.dna.toplevel.fa -s " "\"lcl|"$3"\" -S " "\""$4"\" -L "arr[x]","arr[x+1]" ")}' | awk -F ">" '{if ($2 != "") print ">"$2"@"$1; else print $0}' | sed 's/ .*@/@/' > _g.intra

# search for matching small rna in intra-APA
python -c "
import sys
threshold = 50
length = 18
length_max = 30
file = open('SRR099267.bed', 'r')
apa = open('../_WT-CM-X.intra', 'r')
record_chrx = ''
curr__stop = 0
record_start = 0
recording = False
hits = []
for line in file:
    (chrx, start, stop, val) = line.strip().split('\t')
    val = int(val)
    start = int(start)
    stop = int(stop)
    if recording:
        if chrx == record_chrx and start == curr_stop and val >= threshold:
            curr_stop = stop
        else:
            if (stop - record_start) >= length:
                
                if (stop - record_start) <= length_max:
                    #print record_chrx, record_start, stop, (stop - record_start)
                    hits.append((record_chrx, record_start, curr_stop, (curr_stop - record_start)))
            recording = False
    else:
        if val >= threshold:
            record_start = start
            record_chrx = chrx
            curr_stop = stop
            recording = True   
#exit(0)            
for line in apa:
    (gene, chrx, num, pos) = line.strip().split(' ')
    pos = pos.split(':')
    for (record_chrx, record_start, record_stop, length) in hits:
        if chrx == record_chrx:
            pos = [int(p) for p in pos]
            for i, p in enumerate(pos[:-1]):
                if (record_start >= pos[i] and record_start <= pos[i+1]) and (record_stop >= pos[i] and record_stop <= pos[i+1]):
                    print gene #, record_start, record_stop
                
" | awk '{system("fastacmd -d ../Magnaporthe_oryzae.MG8.18.dna.toplevel.fa -s " "\"lcl|"$1"\" -L "$2","$3" ")}' | fastx_collapser > _micro
#| awk '{system("fastacmd -d ../Magnaporthe_oryzae.MG8.18.dna.toplevel.fa -s " "\"lcl|"$1"\" -L "$2","$3" ")}' > SRR643875_.fasta
#blastn -task blastn-short -query SRR643875_.fasta -db _WT-CM-X.intra -outfmt 6 -max_target_seqs 1 | awk '{if ($4 >=20) print $0}'


# orphans (400 or 1000 nt) search against known db
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.21.dna.toplevel.fa WT-ALL-X.notpolyA_all_m_high  -400 0 print WT-ALL-X.notpolyA_all_m_high_400.fa
blastn -task dc-megablast -query WT-ALL-X.notpolyA_all_m_high_400.fa -db nt -remote -outfmt 5 | python ../../m-oryzae-polya/parse_blast_xml.py 
blastn -task dc-megablast -query WT-ALL-X.notpolyA_all_m_high_400.fa -db Rfam.fasta -outfmt 5 | python ../../m-oryzae-polya/parse_blast_xml.py

python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa WT-ALL-X.notpolyA_all_m_low  -218 0 print _t
formatdb -i _t -p F -o
blastn -query  small/CPA_sRNA.fa -db _t -outfmt "6 qseqid sseqid pident qlen length evalue" -max_target_seqs 1 | awk '{if ($3 >= 98 && $5/$4 >=0.8) print $0}' | cut -f 2 | sort | uniq | wc -l
# perl ../../m-oryzae-polya/rfam_scan.pl --nobig -blastdb Rfam.fasta Rfam.cm WT-ALL-X.notpolyA_all_m_400.fa 

# ophans overlapping annotated genes (antisense)
python -c "
table = {}
for line in open('Magnaporthe_oryzae.MG8.21.gff3', 'r'):
  (chrx, none_1, feat, start, end, none_2,none_3,none_4,infos) = line.strip().split('\t')
  start = int(start)
  end = int(end)
  if not table.has_key(chrx): table[chrx] = {}
  if feat in ('gene', 'protein_coding_gene', 'pseudogene', 'pseudogenic_tRNA', 'rRNA_gene', 'RNA', 'snoRNA_gene', 'snRNA_gene', 'tRNA_gene'):
    if infos.find('Parent') != -1: continue
    id_start = infos.index('ID=')
    id_end = infos.index(';', id_start)
    gene = infos[id_start + 3: id_end]
    table[chrx][gene] = (start, end)
for line in open('WT-ALL-X.notpolyA_all_m_low', 'r'):
 (val, pos, chrx, sense) = line.strip().split(' ')
 pos = int(pos)
 for gene, coord in table[chrx].items():
   if pos >= coord[0] and pos <= coord[1]:
     #print chrx + ':' + str(pos) + ':' + val + ':' + sense# + '\t' + gene
     print line.strip(), gene
"

# extract polyA trascript sequences with fastacmd
cut -f 1,2,3,4 WT_2D4_CM_polyA.diff | awk '{system("fastacmd -d Magnaporthe_oryzae.MG8.18.dna.toplevel.fa -s " "\"lcl|"$2"\" -L "$3","$4" -o blast/"$1".fasta ")}'
# call blast
for f in `ls blast/*fasta`;
do 
blastall -p blastx -m 8 -a 4 -e 1e-20 -i $f -d MG8_proteome.fasta | head -n 1 | cut -f 2 | awk -v f="${f%%.*}" '{system("fastacmd -d MG8_proteome.fasta -s  " "\""$1"\" > "f".hits.fasta  ")}'
done
# call interproscan
for f in `ls *.fasta`;
do
python /media/marco/Elements/m-oryzae-polya/iprscan_soappy.py --email=marco.marconi@gmail.com --title=marco --sequence=$f  --outfile=$f --outformat=out &
done









# let's dump some info

# sequencing infos
ls ../oryzae_reads/2013-06-10-C24C1ACXX/*_1_* | grep -e "2D4....." -e "WT....." -o
zgrep -c "^+$" ../oryzae_reads/2013-06-10-C24C1ACXX/*_1_*
for f in *filtered.bam; do sam-stats $f | grep "mapped reads"; done
for f in *filtered.bam; do sam-stats $f | grep "len mean"; done


# replicates correlation
for f in "WT-CM" "WT-MM" "WT--N" "WT--C" "2D4-CM" "2D4-MM" "2D4--N" "2D4--C"
do
   echo -ne $f
   Rscript ../../m-oryzae-polya/correlation.R $f
done

# number of expressed genes (by replicate)
for sample in *.expr
do
    echo -ne "${sample%%.*}"","
    cat $sample | awk '{if ($2 >= 10) print $1}' |  wc -l
done
# number of expressed genes (by condition) 1 read in at least 2 replicates
for sample in "2D4-CM" "2D4-MM" "2D4--N" "2D4--C" "WT-CM" "WT-MM" "WT--N" "WT--C" 
do
    echo -ne "${sample%%.*}"","
    cat $sample-1.expr $sample-2.expr $sample-3.expr | awk '{if ($2 >= 10) print $1}' | sort | uniq -c | awk '{if ($1 >= 2) print $0}' | wc -l
done


# number of genes with an recognizable generic polyA site
for f in *X.polyA_all_m 
do
    echo -ne "${f%%??.*}"","
    cut -f 5 -d " " $f | sort | uniq | wc -l
done
# number of genes with an recognizable single polyA site
for f in *X.polyA_sgl_m 
do
    echo -ne "${f%%??.*}"","
    cut -f 5 -d " " $f | sort | uniq | wc -l
done
# number of genes with an recognizable apa polyA site
for f in *X.polyA_apa_m 
do
    echo -ne "${f%%??.*}"","
    cut -f 5 -d " " $f | sort | uniq | wc -l
done

# always expressed genes (genes that always have a recognizable polyA) (WT)
cut -f 5 -d " " WT-CM-X.polyA_all_m | sort | uniq > _CM
cut -f 5 -d " " WT-MM-X.polyA_all_m | sort | uniq > _MM
cut -f 5 -d " " WT--N-X.polyA_all_m | sort | uniq > _N
cut -f 5 -d " " WT--C-X.polyA_all_m | sort | uniq > _C
cat _CM _MM _N _C | sort | uniq -c | awk '{if($1 == 4) print $2}' | wc -l
# never expressed genes (genes that never have a recognizable polyA)(WT)
cat _CM _MM _N _C | sort | uniq > _t 
cut -f 1 gene_summary.txt > _all
cat _t _t _all | sort | uniq -u | wc -l
rm _CM _MM _N _C
# always expressed genes (genes that always have a recognizable polyA) (2D4)
cut -f 5 -d " " 2D4-CM-X.polyA_all_m | sort | uniq > _CM
cut -f 5 -d " " 2D4-MM-X.polyA_all_m | sort | uniq > _MM
cut -f 5 -d " " 2D4--N-X.polyA_all_m | sort | uniq > _N
cut -f 5 -d " " 2D4--C-X.polyA_all_m | sort | uniq > _C
cat _CM _MM _N _C | sort | uniq -c | awk '{if($1 == 4) print $2}' | wc -l
# never expressed genes (genes that never have a recognizable polyA) (2D4)
cat _CM _MM _N _C | sort | uniq > _t 
cut -f 1 gene_summary.txt > _all
cat _t _t _all | sort | uniq -u | wc -l
rm _CM _MM _N _C



# distribution of APA in genes 
for f in *X.polyA_all_m
do
   echo -ne "${f%%??.*},"	
   python ../../m-oryzae-polya/polyA_distribution.py Magnaporthe_oryzae.MG8.21.gff3 $f 
done

# number of genes with APA in WT but not in 2D4, and viceversa
for f in "ALL" "CM" "MM" "-N" "-C" 
do
    echo -ne "$f,"
    cut WT-"$f"-X.polyA_apa_m  -d " " -f 5 | sort | uniq > _a
    cut 2D4-"$f"-X.polyA_apa_m  -d " " -f 5 | sort | uniq > _b
    a=$(cat _a _b _b | sort | uniq -u | wc -l)
    b=$(cat _a _a _b | sort | uniq -u | wc -l)
    echo $a","$b
done


# APA classification
for f in *X.polyA_apa_m
do
   echo -ne "${f%%??.*},"	
   cut -f 5 -d " " $f | sort | uniq -c | awk '{ arr[$1]++; count++; if ($1 > max) max = $1 } END { for (i=2;i<= max;i++) printf "%d,",arr[i] }'
   echo	
done


# number of cleveage sites per gene 
for f in *X.polyA_all_m
do
    echo -ne "${f%%.*}"","
    cut -f 5 -d " " $f | sort | uniq -c | awk '{count+=1;num+=$1} END {print num/count}'
done

# differential expressed genes number
for f in diff_expr/*down.csv
do
 echo -ne $(basename "${f%%.*}")"," | sed 's/_expr_down//'
 cat $f | wc -l
done
for f in diff_expr/*up.csv
do
 echo -ne $(basename "${f%%.*}")"," | sed 's/_expr_up//'
 cat $f | wc -l
done





# number of poly-a sites by replicate
for f in *[123].polyA_all_m
do
    echo -ne "${f%%.*}"","
    wc -l $f | awk '{print $1}'
done    

# number of poly-a sites comulative
for f in *X.polyA_all_m
do
    echo -ne "${f%%??.*}"","
    wc -l $f | awk '{print $1}'
done    

# localization of polyA sites 
for f in *X.polyA_all_m
do
   echo -ne "${f%%??.*},"	
   python ../../m-oryzae-polya/polyA_localization.py Magnaporthe_oryzae.MG8.21.gff3 $f 
done


# orphan polyA
for f in "WT-ALL" "WT-CM" "WT-MM" "WT--N" "WT--C" "2D4-ALL" "2D4-CM" "2D4-MM" "2D4--N" "2D4--C"
do
    echo -ne "${f%%.*}"","
    #normal=$(cat "${f%%.*}"-1.polyA "${f%%.*}"-2.polyA "${f%%.*}"-3.polyA |  sort -k 2,7 | uniq -f 1 -c | awk '{if ($1 >= 2) print 0,$3,$4,$5,$6,$7,$8}' | wc -l)
    #orphan=$(cat "${f%%.*}"-1.notpolyA "${f%%.*}"-2.notpolyA "${f%%.*}"-3.notpolyA |  sort -k 2,7 | uniq -f 1 -c | awk '{if ($1 >= 2) print 0,$3,$4,$5,$6,$7,$8}' | wc -l)
    #echo 'scale=3;'$orphan/\($normal+$orphan\)*100 | bc | sed 's/0*$/%,/' | tr -d "\n"
    #echo 'scale=3;'$normal/\($normal+$orphan\)*100 | bc | sed 's/0*$/%/'
    normal=$(cat "${f%%.*}"-X.polyA_all_m | wc -l )
    orphan=$(cat "${f%%.*}"-X.notpolyA_all_m_high | wc -l )
    echo $normal","$orphan
done

# orphans differentially expressed
for f in notdiff_polyA_high/WT*WT*.csv
do
	echo -ne $(basename "${f%%.*}")"," | sed 's/_notpolyA//'
 	cat $f | awk -F "," '{if ($7 < 0.05) print $1}' | wc -l
done

# differential polyA sites number
# P1
for f in diff_polyA/*down*polyA_all_m
do
 echo -ne $(basename "${f%%.*}")"," | sed 's/_down//'
 cat $f | wc -l
done
# P2
for f in diff_polyA/*up*polyA_all_m
do
 echo -ne $(basename "${f%%.*}")"," | sed 's/_up//'
 cat $f | wc -l
done
# differential polyA sites number (genes)
# G1
for f in diff_polyA/*down*polyA_all_m
do
 echo -ne $(basename "${f%%.*}")"," | sed 's/_down//'
 cut $f -d " " -f 5 | sort | uniq | wc -l
done
# G2
for f in diff_polyA/*up*polyA_all_m
do
 echo -ne $(basename "${f%%.*}")"," | sed 's/_up//'
 cut $f -d " " -f 5 | sort | uniq | wc -l
done
# G3
for f in diff_polyA/*polyA.csv
do
 echo -ne $(basename "${f%%_polyA.*}")","
 cat "${f%%_polyA.*}"_down.polyA_all_m  "${f%%_polyA.*}"_up.polyA_all_m | cut -f 5 -d " " | sort | uniq  | wc -l
done
# genes differentially expressed with differentially expressed polyA
for f in diff_polyA/*polyA.csv
do
 echo -ne $(basename "${f%%_polyA.*}")","
 a=diff_expr/$(basename "${f%%_polyA.*}")"_expr"
 cat "${f%%_polyA.*}"_up.polyA_all_m  "${f%%_polyA.*}"_down.polyA_all_m | cut -f 5 -d " " | sort | uniq  > _t
 up=`cat $a"_up.csv" | cut -f 1 -d "," | sort | uniq | xargs -ipat grep pat _t | wc -l ` 
 down=`cat $a"_down.csv"  | cut -f 1 -d "," | sort | uniq | xargs -ipat grep pat _t | wc -l `
 echo $up,$down
done
# G1 sgl and APA
for f in diff_polyA/*_polyA.csv
do
 echo -ne $(basename "${f%%_polyA.*}")","
 cut "${f%%_polyA.*}"_down.polyA_all_m  -d " " -f 5 | sort | uniq > _g1
 cat $(basename "${f%%_vs_*}")"-X.polyA_apa_m" | cut -d " " -f 5 | sort | uniq > _apa
 cat $(basename "${f%%_vs_*}")"-X.polyA_sgl_m" | cut -d " " -f 5 | sort | uniq > _sgl
 echo -ne `cat _g1 _sgl | sort | uniq -d | wc -l`","
 echo `cat _g1 _apa | sort | uniq -d | wc -l`
done
# G2 sgl and APA
for f in diff_polyA/*_polyA.csv
do
 echo -ne $(basename "${f%%_polyA.*}")","
 cut "${f%%_polyA.*}"_up.polyA_all_m -d " " -f 5 | sort | uniq > _g3
 cat $(basename "${f%%_vs_*}")"-X.polyA_apa_m" | cut -d " " -f 5 | sort | uniq > _apa
 cat $(basename "${f%%_vs_*}")"-X.polyA_sgl_m" | cut -d " " -f 5 | sort | uniq > _sgl
 echo -ne `cat _g3 _sgl | sort | uniq -d | wc -l`","
 echo `cat _g3 _apa | sort | uniq -d | wc -l`
done
# G3 sgl and APA
for f in diff_polyA/*_polyA.csv
do
 echo -ne $(basename "${f%%_polyA.*}")","
 cut "${f%%_polyA.*}"_down.polyA_all_m  "${f%%_polyA.*}"_up.polyA_all_m -d " " -f 5 | sort | uniq > _g3
 cat $(basename "${f%%_vs_*}")"-X.polyA_apa_m" | cut -d " " -f 5 | sort | uniq > _apa
 cat $(basename "${f%%_vs_*}")"-X.polyA_sgl_m" | cut -d " " -f 5 | sort | uniq > _sgl
 echo -ne `cat _g3 _sgl | sort | uniq -d | wc -l`","
 echo `cat _g3 _apa | sort | uniq -d | wc -l`
done


# polyA site usage change (only dependent)
for f in "CM" # "MM" 
do
    echo -ne WT-$f"_"vs_WT--C","
    cat diff_polyA/WT-$f"_"vs_WT--C_down.polyA_all_m  diff_polyA/WT-$f"_"vs_WT--C_up.polyA_all_m > _a
    python ../../m-oryzae-polya/polyA_localization.py Magnaporthe_oryzae.MG8.21.gff3 _a | grep "3'UTR" | awk '{print $5",E"$3"@"$2}' | grep -f - diff_polyA/WT-$f"_"vs_WT--C_polyA.csv > _g
    cat _g | python ../../m-oryzae-polya/polyA_usage_ratio.py # | awk '{x=($2-$4); print x < 0 ? -x : x , $5 }' 
done
for f in "CM" #"MM" "-N" "-C"
do
    echo -ne WT"-"$f"_"vs_2D4"-"$f","
    cat diff_polyA/WT"-"$f"_"vs_2D4"-"$f"_"down.polyA_all_m  diff_polyA/WT"-"$f"_"vs_2D4"-"$f"_"up.polyA_all_m > _a
    python ../../m-oryzae-polya/polyA_localization.py Magnaporthe_oryzae.MG8.21.gff3 _a | grep "3'UTR" | awk '{print $5",E"$3"@"$2}' | grep -f - diff_polyA/WT"-"$f"_"vs_2D4"-"$f"_"polyA.csv > _g
    cat _g | python ../../m-oryzae-polya/polyA_usage_ratio.py #| awk '{x=($2-$4); print x < 0 ? -x : x , $5 }' > _t
    #Rscript ../../m-oryzae-polya/plot_fold.R  _t
done
for f in "CM"  "MM" 
do
    echo -ne 2D4-$f"_"vs_2D4--C","
    cat diff_polyA/2D4-$f"_"vs_2D4--C_down.polyA_all_m  diff_polyA/2D4-$f"_"vs_2D4--C_up.polyA_all_m > _a
    python ../../m-oryzae-polya/polyA_localization.py Magnaporthe_oryzae.MG8.21.gff3 _a | grep "3'UTR" | awk '{print $5",E"$3"@"$2}' | grep -f - diff_polyA/2D4-$f"_"vs_2D4--C_polyA.csv > _g
    cat _g | python ../../m-oryzae-polya/polyA_usage_ratio.py | awk '{x=($2-$4); print x < 0 ? -x : x , $5 }' > _t
done


# 3' UTR length
for f in *X.polyA_all_m
do
    echo -ne "${f%%??.*}"","	
    python ../../m-oryzae-polya/3UTR_length.py Magnaporthe_oryzae.MG8.21.gff3 $f
done 

# 3' UTR length (cumulative only in affected genes)
for f in "CM" #"MM" "-N" "-C"
do
    echo "WT vs 2D4 "$f
    cat diff_polyA/WT"-"$f"_"vs"_"2D4"-"$f"_"down.polyA_all_m   | cut -f 5 -d " " | sort | uniq > _diff
    cat _diff | xargs -ipat grep pat WT-$f-X.polyA_all_m  > _t
    python ../../m-oryzae-polya/3UTR_length.py Magnaporthe_oryzae.MG8.18.gff3 _t
    cat _diff | xargs -ipat grep pat 2D4-$f-X.polyA_all_m > _t
    python ../../m-oryzae-polya/3UTR_length.py Magnaporthe_oryzae.MG8.18.gff3 _t    
done

for f in "CM" #"MM" "-N" "-C"
do
    echo "WT vs 2D4 "$f
    cat diff_polyA/WT"-"$f"_"vs"_"2D4"-"$f"_"up.polyA_all_m   | cut -f 5 -d " " | sort | uniq > _diff
    cat _diff | xargs -ipat grep pat WT-$f-X.polyA_all_m  > _t
    python ../../m-oryzae-polya/3UTR_length.py Magnaporthe_oryzae.MG8.18.gff3 _t
    cat _diff | xargs -ipat grep pat 2D4-$f-X.polyA_all_m > _t
    python ../../m-oryzae-polya/3UTR_length.py Magnaporthe_oryzae.MG8.18.gff3 _t    
done

for f in "CM" "MM" "-N" "-C"
do
    echo "WT vs 2D4 "$f
    cat diff_polyA/WT"-"$f"_"vs"_"2D4"-"$f"_"down.polyA_all_m diff_polyA/WT"-"$f"_"vs"_"2D4"-"$f"_"up.polyA_all_m   | cut -f 5 -d " " | sort | uniq > _diff
    cat _diff | xargs -ipat grep pat WT-$f-X.polyA_all_m  > _t
    python ../../m-oryzae-polya/3UTR_length.py Magnaporthe_oryzae.MG8.18.gff3 _t
    cat _diff | xargs -ipat grep pat 2D4-$f-X.polyA_all_m > _t
    python ../../m-oryzae-polya/3UTR_length.py Magnaporthe_oryzae.MG8.18.gff3 _t    
done





# APA average distance
for f in `ls WT-CM-X.polyA_apa_m`
do
echo -ne "${f%%.*}"" "	
python -c "
fp = open('$f', 'r')
(val, curr_pos, chrx, sense, curr_gene, start, end) = fp.readline().strip().split(' ')
curr_pos = int(curr_pos)
dist = 0
dists = {i: 0 for i in range(100000)}
count = 0
for line in fp:
 (val, pos, chrx, sense, gene, start, end) = line.strip().split(' ')
 pos = int(pos)
 if gene == curr_gene:
  dist += abs(pos - curr_pos)
  dists[abs(pos - curr_pos)] += 1
  count += 1
  curr_pos = pos
 else:
  curr_gene = gene
  curr_pos = pos      
print dist / float(count)
for k,v in dists.items():
 print k, v
 if k == 500: break
"
done











# DOMAINS

python ../m-oryzae-polya/fasta_squash.py HRP1_seqdump.fasta HRP1_seqdump_nr.fasta > table
muscle -in HRP1_seqdump_nr.fasta -phyi > HRP1_seqdump_nr.phy
./PhyML-3.1_linux64 -i HRP1_seqdump_nr.phy -d aa
python ../../m-oryzae-polya/fasta_split.py ../HRP1_seqdump_nr.fasta
for f in `ls *.fasta`; do python ../../m-oryzae-polya/iprscan_soappy.py --email=marco.marconi@gmail.com --title=marco --sequence=$f  --outfile="${f%%.*}" --outformat=out & done
rm ../HRP1_domains.txt; for f in `ls *txt`; do grep HMMPfam $f >> ../HRP1_domains.txt ; done
sed -f table HRP1_seqdump_nr.phy_phyml_tree.txt > HRP1.tree
grep -o -e "[A-Za-z].*_.*" HRP1_webtree.txt > HRP1_order.txt
python ../m-oryzae-polya/draw_domains.py HRP1_domains.txt HRP1_order.txt












# basurilla


fastq-mcf -o 2D4--C-1_1_trimmed.fastq -o 2D4--C-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4--C-1_1.fastq 2D4--C-1_2.fastq 
fastq-mcf -o 2D4-CM-1_1_trimmed.fastq -o 2D4-CM-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4-CM-1_1.fastq 2D4-CM-1_2.fastq 
fastq-mcf -o 2D4-CM-2_1_trimmed.fastq -o 2D4-CM-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4-CM-2_1.fastq 2D4-CM-2_2.fastq 
fastq-mcf -o 2D4-MM-1_1_trimmed.fastq -o 2D4-MM-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4-MM-1_1.fastq 2D4-MM-1_2.fastq 
fastq-mcf -o 2D4-MM-2_1_trimmed.fastq -o 2D4-MM-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4-MM-2_1.fastq 2D4-MM-2_2.fastq 
fastq-mcf -o 2D4--N-1_1_trimmed.fastq -o 2D4--N-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4--N-1_1.fastq 2D4--N-1_2.fastq 
fastq-mcf -o WT--C-1_1_trimmed.fastq -o WT--C-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT--C-1_1.fastq WT--C-1_2.fastq 
fastq-mcf -o WT-CM-1_1_trimmed.fastq -o WT-CM-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT-CM-1_1.fastq WT-CM-1_2.fastq
fastq-mcf -o WT-CM-2_1_trimmed.fastq -o WT-CM-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT-CM-2_1.fastq WT-CM-2_2.fastq 
fastq-mcf -o WT-MM-1_1_trimmed.fastq -o WT-MM-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT-MM-1_1.fastq WT-MM-1_2.fastq 
fastq-mcf -o WT-MM-2_1_trimmed.fastq -o WT-MM-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT-MM-2_1.fastq WT-MM-2_2.fastq 
fastq-mcf -o WT--N-1_1_trimmed.fastq -o WT--N-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT--N-1_1.fastq WT--N-1_2.fastq
fastq-mcf -o 2D4--C-2_1_trimmed.fastq -o 2D4--C-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4--C-2_1.fastq 2D4--C-2_2.fastq
fastq-mcf -o 2D4--C-3_1_trimmed.fastq -o 2D4--C-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4--C-3_1.fastq 2D4--C-3_2.fastq 
fastq-mcf -o 2D4-CM-3_1_trimmed.fastq -o 2D4-CM-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4-CM-3_1.fastq 2D4-CM-3_2.fastq 
fastq-mcf -o 2D4-MM-3_1_trimmed.fastq -o 2D4-MM-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4-MM-3_1.fastq 2D4-MM-3_2.fastq
fastq-mcf -o 2D4--N-2_1_trimmed.fastq -o 2D4--N-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4--N-2_1.fastq 2D4--N-2_2.fastq
fastq-mcf -o 2D4--N-3_1_trimmed.fastq -o 2D4--N-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4--N-3_1.fastq 2D4--N-3_2.fastq 
fastq-mcf -o WT--C-2_1_trimmed.fastq -o WT--C-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT--C-2_1.fastq WT--C-2_2.fastq
fastq-mcf -o WT--C-3_1_trimmed.fastq -o WT--C-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT--C-3_1.fastq WT--C-3_2.fastq
fastq-mcf -o WT-CM-3_1_trimmed.fastq -o WT-CM-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT-CM-3_1.fastq WT-CM-3_2.fastq 
fastq-mcf -o WT-MM-3_1_trimmed.fastq -o WT-MM-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT-MM-3_1.fastq WT-MM-3_2.fastq 
fastq-mcf -o WT--N-2_1_trimmed.fastq -o WT--N-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT--N-2_1.fastq WT--N-2_2.fastq
fastq-mcf -o WT--N-3_1_trimmed.fastq -o WT--N-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT--N-3_1.fastq WT--N-3_2.fastq


gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails 2D4--C-1_1_trimmed.fastq 2D4--C-1_2_trimmed.fastq > 2D4--C-1.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails 2D4-CM-1_1_trimmed.fastq 2D4-CM-1_2_trimmed.fastq > 2D4-CM-1.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails 2D4-CM-2_1_trimmed.fastq 2D4-CM-2_2_trimmed.fastq > 2D4-CM-2.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails 2D4-MM-1_1_trimmed.fastq 2D4-MM-1_2_trimmed.fastq > 2D4-MM-1.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails 2D4-MM-2_1_trimmed.fastq 2D4-MM-2_2_trimmed.fastq > 2D4-MM-2.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails 2D4--N-1_1_trimmed.fastq 2D4--N-1_2_trimmed.fastq > 2D4--N-1.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails WT--C-1_1_trimmed.fastq WT--C-1_2_trimmed.fastq > WT--C-1.sam 
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails WT-CM-1_1_trimmed.fastq WT-CM-1_2_trimmed.fastq > WT-CM-1.sam 
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails WT-CM-2_1_trimmed.fastq WT-CM-2_2_trimmed.fastq > WT-CM-2.sam 
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails WT-MM-1_1_trimmed.fastq WT-MM-1_2_trimmed.fastq > WT-MM-1.sam 
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails WT-MM-2_1_trimmed.fastq WT-MM-2_2_trimmed.fastq > WT-MM-2.sam 
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails WT--N-1_1_trimmed.fastq WT--N-1_2_trimmed.fastq > WT--N-1.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails 2D4--C-2_1_trimmed.fastq 2D4--C-2_2_trimmed.fastq > 2D4--C-2.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails 2D4--C-3_1_trimmed.fastq 2D4--C-3_2_trimmed.fastq > 2D4--C-3.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails 2D4-CM-3_1_trimmed.fastq 2D4-CM-3_2_trimmed.fastq > 2D4-CM-3.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails 2D4-MM-3_1_trimmed.fastq 2D4-MM-3_2_trimmed.fastq > 2D4-MM-3.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails 2D4--N-2_1_trimmed.fastq 2D4--N-2_2_trimmed.fastq > 2D4--N-2.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails 2D4--N-3_1_trimmed.fastq 2D4--N-3_2_trimmed.fastq > 2D4--N-3.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails WT--C-2_1_trimmed.fastq WT--C-2_2_trimmed.fastq > WT--C-2.sam 
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails WT--C-3_1_trimmed.fastq WT--C-3_2_trimmed.fastq > WT--C-3.sam 
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails WT-CM-3_1_trimmed.fastq WT-CM-3_2_trimmed.fastq > WT-CM-3.sam 
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails WT-MM-3_1_trimmed.fastq WT-MM-3_2_trimmed.fastq > WT-MM-3.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails WT--N-2_1_trimmed.fastq WT--N-2_2_trimmed.fastq > WT--N-2.sam
gsnap -B 5 -t 8 -A sam -d MG8_21 -D ./MG8_21/  --nofails WT--N-3_1_trimmed.fastq WT--N-3_2_trimmed.fastq > WT--N-3.sam



# P1 vs notP1
cat diff_polyA/WT-CM_vs_2D4-CM_down.polyA_all_m diff_polyA/WT-MM_vs_2D4-MM_down.polyA_all_m diff_polyA/WT--N_vs_2D4--N_down.polyA_all_m diff_polyA/WT--C_vs_2D4--C_down.polyA_all_m | sort -k 1,7 | uniq > _p1
cat WT-CM-X.polyA_all_m WT-MM-X.polyA_all_m WT--N-X.polyA_all_m WT--C-X.polyA_all_m | sort -k 1,7 | uniq > _t
cat _p1 _p1 _t | sort -k 1,7 | uniq -u | sort -R > _not_p1
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.21.dna.toplevel.fa _p1 -100 -30 print _s1
c=$(cat _p1 | wc -l)
head -n $c _not_p1 > _not_p1_1
tail -n $c _not_p1 > _not_p1_2
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.21.dna.toplevel.fa _not_p1_1 -100 -30  print _s2
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.21.dna.toplevel.fa _not_p1_2 -100 -30  print _s3

# P1 vs P2 vs P1_others
cat diff_polyA/WT-CM_vs_2D4-CM_down.polyA_all_m diff_polyA/WT-MM_vs_2D4-MM_down.polyA_all_m diff_polyA/WT--N_vs_2D4--N_down.polyA_all_m diff_polyA/WT--C_vs_2D4--C_down.polyA_all_m | sort -k 1,7 | uniq > _p1
cat diff_polyA/WT-CM_vs_2D4-CM_up.polyA_all_m diff_polyA/WT-MM_vs_2D4-MM_up.polyA_all_m diff_polyA/WT--N_vs_2D4--N_up.polyA_all_m diff_polyA/WT--C_vs_2D4--C_up.polyA_all_m | sort -k 1,7 | uniq > _p2
cat WT-CM-X.polyA_all_m WT-MM-X.polyA_all_m WT--N-X.polyA_all_m WT--C-X.polyA_all_m | sort -k 1,7 | uniq > _t
cut -f 5 -d " " _p1 | xargs -ipat grep pat _t > _p1_all
cat _p1 _p1 _p2 _p2 _p1_all | sort -k 1,7 | uniq -u > _p1_others
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.21.dna.toplevel.fa _p1 -100 100 print _s1
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.21.dna.toplevel.fa _p2 -100 100 print _s2
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.21.dna.toplevel.fa _p1_others -100 100 print _s1_others


for f in `ls diff_polyA/*.csv`
do
    n=${f/polyA.csv/}
    cat $n"down.polyA_all_m"  $n"up.polyA_all_m" > _a
     python ../../m-oryzae-polya/polyA_localization.py Magnaporthe_oryzae.MG8.21.gff3 _a | grep "3'UTR" | awk '{print $5",E"$3"@"$2}' | grep -f - $f > _g
    cat _g | python ../../m-oryzae-polya/polyA_usage_ratio.py | awk '{if ($3 < 0) print $1}' > "_"$(basename $n)"_polyA_short" 
    cat _g | python ../../m-oryzae-polya/polyA_usage_ratio.py | awk '{if ($3 > 0) print $1}' > "_"$(basename $n)"_polyA_long" 
    n=$(sed  -e 's/diff_polyA\///' -e 's/_polyA.*//' <<< $f)
    arrIN=(${n//_vs_/ })
    cut -f 5 -d " " ${arrIN[0]}"-X.polyA_all_m" ${arrIN[1]}"-X.polyA_all_m" | sort | uniq > "_"$(basename $n)"_polyA_back"

done


	
# extract up and down expr changed genes
for f in `ls diff_expr/*down.csv`
do
    cat $f | cut -f 1,3 -d "," | sed 's/,/\t/' | sort -k2 -n | cut -f 1 > "_"$(basename ${f/.csv/}) 
done
for f in `ls diff_expr/*up.csv`
do
    cat $f | cut -f 1,3 -d "," | sed 's/,/\t/' | sort -k2 -rn | cut -f 1 > "_"$(basename ${f/.csv/}) 
done
for f in `ls diff_expr/*_expr.csv`
do
    awk -F "," '{if($2!="0") print $1}' < $f > "_"$(basename ${f/.csv/})"_back"
done



# extract by distance
python -c "

p1 = open('_p1', 'r')
p2 = open('_p2', 'r')
t1 = {}
diffs = [0 for i in range(500)]
for l1 in p1:
    items = l1.strip().split(' ')
    gene = items[4]
    sense = items[3]
    pos = int(items[1])
    t1[gene] = pos
for l2 in p2:
    items = l2.strip().split(' ')
    gene = items[4]
    sense = items[3]
    pos = int(items[1])
    if t1.has_key(gene):
        if abs(pos - t1[gene]) > 150:
                print gene 


"

# extract conserved orphans from blast
python -c "
from Bio import SeqIO

g = {}
for line in open('_g', 'r'):
    items = line.strip().split('\t')
    orphan = items[0]
    chrx = items[1]
    start = items[8]
    end = items[9]
    g[orphan] = (chrx, int(start), int(end))
m = {}
for line in open('_m', 'r'):
    items = line.strip().split('\t')
    orphan = items[0]
    chrx = items[1]
    start = items[8]
    end = items[9]
    m[orphan] = (chrx, int(start), int(end))
import sys    
for line in open('_c', 'r'):
     l = line.strip()
     f = open(l + '.fa', 'w')
     record_ophans = SeqIO.index('../WT-ALL-X.notpolyA_all_m_low_200.fa', 'fasta')
     record_gae = SeqIO.index('../others/Gaeumannomyces_graminis.Gae_graminis_V2.21.dna.genome.fa', 'fasta')
     record_poae = SeqIO.index('../others/Magnaporthe_poae.Mag_poae_ATCC_64411_V1.21.dna.genome.fa', 'fasta') 
     f.write('> oryzae_' + l + '\n') 
     f.write( str(record_ophans[l].seq) + '\n') 
     f.write( '> gae_' + l+ '\n') 
     if g[l][1] < g[l][2]:
         f.write( str(record_gae[g[l][0]].seq[g[l][1]:g[l][2]])+ '\n') 
     else:
        g[l][1] > g[l][2]
        f.write( str(record_gae[g[l][0]].seq[g[l][2]:g[l][1]].reverse_complement())+ '\n') 
     f.write( '> poae_' + l+ '\n') 
     if m[l][1] < m[l][2]:
         f.write( str(record_poae[m[l][0]].seq[m[l][1]:m[l][2]])+ '\n') 
     else:
        m[l][1] > m[l][2]
        f.write( str(record_poae[m[l][0]].seq[m[l][2]:m[l][1]].reverse_complement())+ '\n') 
"



# genes 4 sets venn diagrams
cat _c _d | sort | uniq  > _t
cat _a _b | sort | uniq -d | cat _t _t - | sort | uniq -u > _ab
cat _b _d | sort | uniq  > _t
cat _a _c | sort | uniq -d | cat _t _t - | sort | uniq -u > _ac
cat _c _b | sort | uniq  > _t
cat _a _d | sort | uniq -d | cat _t _t - | sort | uniq -u > _ad
cat _a _d | sort | uniq  > _t
cat _b _c | sort | uniq -d | cat _t _t - | sort | uniq -u > _bc
cat _a _c | sort | uniq  > _t
cat _b _d | sort | uniq -d | cat _t _t - | sort | uniq -u > _bd
cat _a _b | sort | uniq  > _t
cat _c _d | sort | uniq -d | cat _t _t - | sort | uniq -u > _cd
cat _a | sort | uniq  > _t
cat _b _c _d | sort | uniq -c | awk '{if ($1 == 3) print $2}' | cat _t _t - | sort | uniq -u > _bcd
cat _b | sort | uniq  > _t
cat _a _c _d | sort | uniq -c | awk '{if ($1 == 3) print $2}' | cat _t _t - | sort | uniq -u > _acd
cat _c | sort | uniq  > _t
cat _a _b _d | sort | uniq -c | awk '{if ($1 == 3) print $2}' | cat _t _t - | sort | uniq -u > _abd
cat _d | sort | uniq  > _t
cat _a _b _c | sort | uniq -c | awk '{if ($1 == 3) print $2}' | cat _t _t - | sort | uniq -u > _abc
cat _a _b _c | sort | uniq  > _t
cat _d | sort | uniq  | cat _t _t - | sort | uniq -u > __d
cat _a _b _d | sort | uniq  > _t
cat _c | sort | uniq  | cat _t _t - | sort | uniq -u > __c
cat _a _d _c | sort | uniq  > _t
cat _b | sort | uniq  | cat _t _t - | sort | uniq -u > __b
cat _d _b _c | sort | uniq  > _t
cat _a | sort | uniq  | cat _t _t - | sort | uniq -u > __a
cat _a _b _c _d | sort | uniq -c | awk '{if ($1 == 4) print $2}' > _abcd

# biomart join
python -c "
import sys
arr = {}
for line in open('_e.csv', 'r'):
  data =  line.strip().split('\t')
  gene = data[0]
  all = data[1:]
  if not arr.has_key(gene):
    arr[gene] = [set((x,)) for x in all]
  else:
    for i, x in enumerate(arr[gene]):
      arr[gene][i].add(all[i])

for k,v in arr.items():
  sys.stdout.write(k + '\t')
  for x in v: 
    sys.stdout.write('; '.join(x))
    sys.stdout.write('\t')
  print
" > _f.csv

# extract start codon
python -c "
import re
table = {}
for line in open('Magnaporthe_oryzae.MG8.21.gff3', 'r'):
    if line[0] == '#' or line[0] == '\n':
        continue
    items = line.split('\t')
    sense = items[6]
    if items[2]  == 'protein_coding_gene':
        chrx = items[0]
        start = int(items[3])
        end = int(items[4])
        sense = items[6]
        for x in items[8].split(';'):
            if x.split('=')[0] == 'ID':
                gene = x.split('=')[1].strip()
        table[gene] = None
    elif items[2] == 'CDS':
        chrx = items[0]
        start = int(items[3])
        end = int(items[4])
        sense = items[6]
        for x in items[8].split(';'):
            if x.split('=')[0] == 'Parent':
                gene = x.split('=')[1].strip()
                gene= re.sub('T.', '', gene)
        if  sense == '+' and table.has_key(gene):
            if table[gene] == None:
                table[gene] = (chrx, start, sense)
        elif  sense == '-' and table.has_key(gene):
             table[gene] = (chrx, end, sense)
for gene, (chrx, pos, sense) in table.items():
 if  sense == '+': 
  start = pos
  end = pos + 3
 elif  sense == '-':
  start = pos - 3
  end = pos  
 print chrx, '\t', '.', '\t', 'start_codon', start, '\t', end, '\t', '.', '\t', sense, '\t', '.', '\t', 'ID='+gene+';Parent='+gene+';'
"
# orphan matching retrotransposons

python -c "
for line in open('WT-ALL-X.notpolyA_all_m_low', 'r'):
    val, pos, chr, sense = line.strip().split(' ')
    for line2 in open('_retro', 'r'):
        chr2, marco, mobile, start, end, p1, sense2, p2, info = line2.strip().split('\t')
        if chr == chr2:
            if sense == '-' and sense2 == '+':
                if int(pos) > int(start) and int(pos) < int(end) + 400:
                    print line.strip(), line2.strip()
                    break
            if sense == '+' and sense2 == '-':
                if int(pos) > int(start) - 400 and int(pos) < int(end):
                    print line.strip(), line2.strip()
                    break
" 

#kegg
cut -f 5 -d " " ../2D4-CM-X.polyA_all_m ../2D4--C-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _2D4-CM_vs_2D4--C_expr_down _back
cut -f 5 -d " " ../2D4-CM-X.polyA_all_m ../2D4--C-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _2D4-CM_vs_2D4--C_expr_up _back
cut -f 5 -d " " ../2D4-CM-X.polyA_all_m ../2D4-MM-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _2D4-CM_vs_2D4-MM_expr_down _back
cut -f 5 -d " " ../2D4-CM-X.polyA_all_m ../2D4-MM-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _2D4-CM_vs_2D4-MM_expr_up _back
cut -f 5 -d " " ../2D4-CM-X.polyA_all_m ../2D4--N-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _2D4-CM_vs_2D4--N_expr_down _back
cut -f 5 -d " " ../2D4-CM-X.polyA_all_m ../2D4--N-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _2D4-CM_vs_2D4--N_expr_up _back
cut -f 5 -d " " ../2D4-MM-X.polyA_all_m ../2D4--C-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _2D4-MM_vs_2D4--C_expr_down _back
cut -f 5 -d " " ../2D4-MM-X.polyA_all_m ../2D4--C-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _2D4-MM_vs_2D4--C_expr_up _back
cut -f 5 -d " " ../2D4-MM-X.polyA_all_m ../2D4--N-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _2D4-MM_vs_2D4--N_expr_down _back
cut -f 5 -d " " ../2D4-MM-X.polyA_all_m ../2D4--N-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _2D4-MM_vs_2D4--N_expr_up _back
cut -f 5 -d " " ../WT-CM-X.polyA_all_m ../2D4-CM-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-CM_vs_2D4-CM_expr_down _back
cut -f 5 -d " " ../WT-CM-X.polyA_all_m ../2D4-CM-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-CM_vs_2D4-CM_expr_up _back
cut -f 5 -d " " ../WT-CM-X.polyA_all_m ../WT--C-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-CM_vs_WT--C_expr_down _back
cut -f 5 -d " " ../WT-CM-X.polyA_all_m ../WT--C-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-CM_vs_WT--C_expr_up _back
cut -f 5 -d " " ../WT-CM-X.polyA_all_m ../WT-MM-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-CM_vs_WT-MM_expr_down _back
cut -f 5 -d " " ../WT-CM-X.polyA_all_m ../WT-MM-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-CM_vs_WT-MM_expr_up _back
cut -f 5 -d " " ../WT-CM-X.polyA_all_m ../WT--N-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-CM_vs_WT--N_expr_down _back
cut -f 5 -d " " ../WT-CM-X.polyA_all_m ../WT--N-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-CM_vs_WT--N_expr_up _back
cut -f 5 -d " " ../WT--C-X.polyA_all_m ../2D4--C-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT--C_vs_2D4--C_expr_down _back
cut -f 5 -d " " ../WT--C-X.polyA_all_m ../2D4--C-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT--C_vs_2D4--C_expr_up _back
cut -f 5 -d " " ../WT-CM-X.polyA_all_m ../2D4-MM-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-MM_vs_2D4-MM_expr_down _back
cut -f 5 -d " " ../WT-CM-X.polyA_all_m ../2D4-MM-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-MM_vs_2D4-MM_expr_up _back
cut -f 5 -d " " ../WT-MM-X.polyA_all_m ../WT--C-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-MM_vs_WT--C_expr_down _back
cut -f 5 -d " " ../WT-MM-X.polyA_all_m ../WT--C-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-MM_vs_WT--C_expr_up _back
cut -f 5 -d " " ../WT-MM-X.polyA_all_m ../WT--N-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-MM_vs_WT--N_expr_down _back
cut -f 5 -d " " ../WT-MM-X.polyA_all_m ../WT--N-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT-MM_vs_WT--N_expr_up _back
cut -f 5 -d " " ../2D4-CM-X.polyA_all_m ../WT--N-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT--N_vs_2D4--N_expr_down _back
cut -f 5 -d " " ../2D4-CM-X.polyA_all_m ../WT--N-X.polyA_all_m | sort | uniq > _back
kegg_enrich ../kegg.txt _WT--N_vs_2D4--N_expr_up _back





# compara create heatmap

 for f in `find . -iname "*.apa_orthologs_scerevisiae"`; do cat $f >> _total; done
 sort _total | uniq > __total
 for f in `find . -iname polyA.apa_orthologs_scerevisiae`; do grep -f __total $f | awk -v f=$f '{print $1,f}' | sed -e 's/\.\///' -e s'/\/.*//'; done > _ass
 awk -F " " '{arr[$1]=arr[$1]","$2} END {for (k in arr) print k,arr[k]}' < _ass | sed  -e 's/_NGS//g' -e 's/ ,/\t/' -e 's/,_test//'  > tree.txt
 python -c "
arr = {}
l = set()
for lin in open('tree.txt', 'r'):
 gene, orths = lin.strip().split('\t')
 arr[gene] = []
 for o in orths.split(','):
    l.add(o)
    arr[gene].append(o)
import sys
sys.stdout.write('gene,')
for s in l:
  sys.stdout.write(s + ',')
print
for gene, orths in arr.items():
  sys.stdout.write(gene + ',')
  for s in l:
     if s in orths:
       sys.stdout.write('1' + ',')
     else:
       sys.stdout.write('0' + ',')
  print
" | sed -e 's/ //' -e 's/,$//' > hist.txt

# compara retrive motif frequency
mkdir _TA
rm _TA/*
#for f in `ls -d */`; do python ../../m-oryzae-polya/polyA_nucleotide.py $f/*genome.fa $f/*polyA -100 100 print $f/"_"sequences; done
count=0; for f in `find . -iname "_sequences"`;  do python scan.py $f TGTA[ATC] > _TA/"${f##*/}"$count; count=$((count+1));  done
	paste _TA/* > _m
Rscript ../../m-oryzae-polya/average_table.R

# compara cutsite frequencies
for f in `ls -d */`; do python ../../m-oryzae-polya/polyA_nucleotide.py $f/*genome.fa $f/*polyA -1 0 print $f/"_"cutsite; done
	count=0; for f in `ls -d */`; do if [ -e $f/"_cutsite" ]; then total=`grep -c ">" $f/"_cutsite"`; echo $f > _TA/_sequence$count ; for m in `cat motifs`; do grep -v ">" $f/"_cutsite" | grep -c $m; done | awk -v tot=$total '{print $1/tot}' >> _TA/_sequence$count;echo $total;  count=$((count + 1)); fi; done
# compara get orthologs from biomart
source=moryzae
dest=hsapiens
genes=`tr '\n' ',' < polyA.apa | sed 's/,$//'`
query=`cat ../query.xml | sed -e 's/SOURCE/'$source'/' -e 's/DEST/'$dest'/' -e 's/VALUE/'$genes'/' | tr -d '\n'`
wget -O _results.txt --post-data "query=$query" "http://fungi.ensembl.org/biomart/martservice" 2> /dev/null > /dev/null
cut -f 2 _results.txt | grep . | sort | uniq > polyA.apa_orthologs_$dest
rm _results.txt



####spliceosome

# get fasta for spliceosome or others
for protein in `cat spliceosome.txt`; do
  rm fasta/$f".fa"
  for species in `cat ensembl_fungi_list.txt`; do
    gene=`grep $protein orthologs/$species".txt" -m 1 | cut -f 2 | grep . `
    if [ "$gene" != "" ]; then  
      query=`cat ../m-oryzae-polya/query_fasta.xml | sed -e 's/SOURCE/'$species'/' -e 's/GENE/'$gene'/' | tr -d '\n'`
      wget -O _results_$protein".txt" --post-data "query=$query" "http://fungi.ensembl.org/biomart/martservice" 2> /dev/null > /dev/null
      cat _results_$protein".txt" | sed 's/>.*/>'$species'/' >> fasta/$protein".fa"
    fi
  done
done

#remove duplicates
for f in `ls fasta/*fa`; 
do
python -c "
from Bio import SeqIO
import sys
s = {}
for seq_record in SeqIO.parse(sys.argv[1], 'fasta'):
    if s.has_key(seq_record.id):
      if len(seq_record.seq) > len(s[seq_record.id]):
         s[seq_record.id] = str(seq_record.seq)
    else:
      s[seq_record.id] = str(seq_record.seq)

for k in sorted(s.keys()):
  print '>', k
  print s[k]

" $f | fasta_formatter -w 80 > _t ; mv _t $f
done 

# create tree pdf
for f in `find . -iname "*tree.txt"`; do while read a b; do sed -i '0,/'$a'/s//'$b'/' $f; done  < labels.txt ; done
for f in `ls *_tree.txt`; do ~/Downloads/newick-utils-1.6/src/nw_display -i 'visibility:hidden' -b 'visibility:hidden' $f -s -S -w 900 | convert - ${f/txt/pdf}; done
rm all_trees.pdf; pdftk *pdf cat output all_trees.pdf

# get all logos
for f in `ls _*`; do ~/Downloads/weblogo/seqlogo -f $f -F PNG -o donor_logo/$f -c -a -n -Y ; done

# alternative splicing
sudo htseq-count -s reverse -r pos -f bam -t gene -i ID SRR081547.sorted.bam Gaeumannomyces_graminis.Gae_graminis_V2.24.gff3 > _gene_count
sudo htseq-count -s reverse -r pos -f bam -t intron -i ID SRR081547.sorted.bam _introns.gtf > _intron_count
awk -F $'\t' 'function abs(x){return ((x < 0.0) ? -x : x)}{if ($3 == "gene") print $9"\t"abs($4-$5)}' < Gaeumannomyces_graminis.Gae_graminis_V2.24.gff3   | sed 's/ID=gene:\(GGTG_.....\).*\t/\1\t/' > _gene_length
awk -F $'\t' 'function abs(x){return ((x < 0.0) ? -x : x)}{if ($3 == "intron") print $9"\t"abs($4-$5)}' < _introns.gtf  | sed 's/.*ID //' > _intron_length
sort -k 1,1 _gene_length -o _gene_length ; sort -k 1,1 _intron_length -o _intron_length
sed -i 's/gene://' _gene_count
sort -k 1,1 _intron_count -o _intron_count; join _intron_count _intron_length | awk '{print $1,$2,$3,$2/$3}' > _intron_norm
sort -k 1,1 _gene_count -o _gene_count; join _gene_count _gene_length | awk '{print $1,$2,$3,$2/$3}' > _gene_norm
rm _o; for f in `cut -f 1 _intron_count`; do grep -m 1 $f  _introns.gtf | sed -e 's/.*gene_id "//' -e 's/"; ID /\t/' >> _o; done; sort -k 1,1 _o -o _o
reset;join _o _gene_norm | sort -k 2,2 | join - _intron_norm -1 2 -2 1 | awk '{if ($3>0 && $5 > 0.05 && $8 > 0.05)print $1,$6/$3}' | scisort -k 2

# count codon/aa
rm _*
for f in `ls rhizopus_oryzae/Rhior3_cds.fa`;
do 

python -c "
import sys
from Bio import SeqIO, Seq
fasta_file = sys.argv[1]  
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
count = {}
aa = {}
for seq in fasta_sequences:
  for codon in [str(seq.seq)[i:i+3] for i in range(0,len(seq.seq),3)]:
    if len(codon) != 3:
      continue
    ok = True  
    for x in codon:
      if x not in ('A', 'T', 'C', 'G'):
        ok = False
    if not ok: continue    
    if not count.has_key(codon):
      count[codon] = 1
    else:
      count[codon] += 1
    a = Seq.translate(codon)
    if not aa.has_key(a):
      aa[a] = {}
      aa[a][codon] = 1
      aa[a]['total'] = 1
    else:
      if not aa[a].has_key(codon):
        aa[a][codon] = 1
        aa[a]['total'] = 1
      else:
        aa[a][codon] += 1 
        aa[a]['total'] += 1    
#for codon, val in count.items():
#  print codon, val     
for aa, codons in aa.items():
  for codon, val in codons.items():
     if codon != 'total': 
      print aa, codon, val / float(codons['total'])
  
" $f | sed 's/ /\t/' > "__"${f/\/*/.info}

done

# domains
for f in `find . -iname "*nw"` ; do g=${f/.\//}; nw_labels -I $f > ${g/ete*/IP}/${g/_ete*/.order}; done
for f in `ls -d ./*_IP`; do cd $f; grep Pfam *tsv -h | sed 's/:/_/g' > ${f/.\//}.Pfam ;  cd ..;done
for d in `ls -d *_IP/`;  do cd $d ; rm ${d/\/}.list; for f in `ls *fa`; do echo -ne ${f/.fa/}"\t" >> ${d/\/}.list;  grep -v ">" $f | tr -d '\n' | wc -c >> ${d/\/}.list; done ; cd ..; done
for f in `find . -iname "*Pfam"` ; do echo $f; if [ -s ${f/_IP.Pfam/.order} ]; then echo $f; python /media/marco/Elements/m-oryzae-polya/draw_domains.py $f ${f/_IP.Pfam/.order}  ${f/.Pfam/.list}  `sed -e 's/.*\///' -e 's/_IP.Pfam//' <<< $f`  5  ; fi; done

	
