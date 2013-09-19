# trim low quality reads, polyA tails and adapters

# fastx_clipper -v -l 17 -a TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -Q33 -i $f   -o temp  
# fastx_clipper -v -l 17  -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCG      -Q33  -i temp -o "${f%%.*}".trimmed.fastq
    
fastq-mcf -o file1_trimmed_1.fastq -o file2_trimmed_2.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa file1.fastq file2.fastq  

    
# build database

gmap_build -d MG8_18 -D ./MG8_18 $d

# align

gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18  a_1.fastq a_2.fastq > "${f%%.*}".sam

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
    python ../../m-oryzae-polya/filter.py "${f%%.*}".sorted.bam Magnaporthe_oryzae.MG8.18.dna.toplevel.fa 7 "${f%%.*}".filtered
done

# assign reads to transcripts

for f in `ls *.filtered.bam`
do
    python ../../m-oryzae-polya/assign.py Magnaporthe_oryzae.MG8.18.gff3 "${f%%.*}".filtered.bam "${f%%.*}".assign "${f%%.*}".notassign 7
done

# create bedgraphs  

for f in `ls *.assign`
do
    cat  "${f%%.*}".assign | cut -f 2,3,4,5 | awk '{  if ( $4 == "-" ) print $3,$1  }' | sort -n | uniq -c | awk '{ printf "%s\t%d\t%d\t%d\n",$3,$2,$2+1,$1 }' > "${f%%.*}"_plus.bedgraph 
    cat  "${f%%.*}".assign | cut -f 2,3,4,5 | awk '{  if ( $4 == "+" ) print $2,$1  }' | sort -n | uniq -c | awk '{ printf "%s\t%d\t%d\t%d\n",$3,$2,$2+1,$1 }' > "${f%%.*}"_minus.bedgraph
done

# create not bedgraphs  

for f in `ls *.notassign`
do
    cat  "${f%%.*}".notassign | cut -f 2,3,4,5 | awk '{  if ( $4 == "-" ) print $3,$1  }' | sort -n | uniq -c | awk '{ printf "%s\t%d\t%d\t%d\n",$3,$2,$2+1,$1 }' > "${f%%.*}"_not_plus.bedgraph 
    cat  "${f%%.*}".notassign | cut -f 2,3,4,5 | awk '{  if ( $4 == "+" ) print $2,$1  }' | sort -n | uniq -c | awk '{ printf "%s\t%d\t%d\t%d\n",$3,$2,$2+1,$1 }' > "${f%%.*}"_not_minus.bedgraph
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

cat Magnaporthe_oryzae.MG8.18.gff3 | awk '{if($3 == "gene") print $0}' | grep  "ID=.*;"  -o | sed -e 's/ID=//' -e 's/;//' > _t
for f in `ls *.assign`
do
    cut -f 6 $f | cat - _t | sort | uniq -c | awk '{print $2"\t"$1-1}' > "${f%%.*}".expr
done
rm _t 


# extract polyA most significant sites

for f in `ls *.polyA`; do python ../../m-oryzae-polya/polyA_extract.py Magnaporthe_oryzae.MG8.18.gff3 $f "${f%%.*}".expr 33 0.05 all > $f"_"all_m; done 

# extract notpolyA most significant sites

for f in `ls *.notpolyA`; do python ../../m-oryzae-polya/polyA_extract_not.py $f 33 100 > $f"_"all_m; done 


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
    cat $f"-"1.notpolyA_all_m $f"-"2.notpolyA_all_m $f"-"3.notpolyA_all_m |  sort -k 2,7 | uniq -f 1 -c | awk '{if ($1 >= 2) print 0,$3,$4,$5,$6,$7,$8}' > $f"-"X.notpolyA_all_m
done
cat WT*X*notpolyA_all_m |  sort -k 1,4 | uniq  > WT-ALL-X.notpolyA_all_m
cat 2D4*X*notpolyA_all_m | sort -k 1,4 | uniq  > 2D4-ALL-X.notpolyA_all_m


# gff of polyA
for f in `ls *.polyA_all_m`; do cut $f -d " " -f 1,2,3,4,5 | awk '{ printf "%s\tmarco\tpolyA_site\t%d\t%d\t.\t%s\t.\ttranscript=%s;value=%d\n", $3, $2, $2, $4, $5, $1  }' > "${f%%.*}"_polyA.gff ; done
for f in `ls *.notpolyA_all_m`; do cut $f -d " " -f 1,2,3,4,5 | awk '{ printf "%s\tmarco\tpolyA_site\t%d\t%d\t.\t%s\t.\ttranscript=%s;value=%d\n", $3, $2, $2, $4, $5, $1  }' > "${f%%.*}"_notpolyA.gff ; done


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
Rscript ../../m-oryzae-polya/diff_expr.R WT--N WT--C    
Rscript ../../m-oryzae-polya/diff_expr.R WT-CM 2D4-CM    
Rscript ../../m-oryzae-polya/diff_expr.R WT-MM 2D4-MM    
Rscript ../../m-oryzae-polya/diff_expr.R WT--N 2D4--N    
Rscript ../../m-oryzae-polya/diff_expr.R WT--C 2D4--C    
Rscript ../../m-oryzae-polya/diff_expr.R 2D4-CM 2D4-MM    
Rscript ../../m-oryzae-polya/diff_expr.R 2D4-CM 2D4--N     
Rscript ../../m-oryzae-polya/diff_expr.R 2D4-CM 2D4--C    
Rscript ../../m-oryzae-polya/diff_expr.R 2D4-MM 2D4--N    
Rscript ../../m-oryzae-polya/diff_expr.R 2D4-MM 2D4--C    
Rscript ../../m-oryzae-polya/diff_expr.R 2D4--N 2D4--C   
wait
for f in diff_expr/*expr.csv
do
 cat $f | awk -F "," '{if ($8 < 0.05 && $6 > 0) print $0}' > "${f%%.*}"_up.csv
 cat $f | awk -F "," '{if ($8 < 0.05 && $6 < 0) print $0}' > "${f%%.*}"_down.csv
done


# differential polyA (p-value < 0.1)

function diff {
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

diff "WT" "2D4" "CM" "CM"
diff "WT" "2D4" "MM" "MM"
diff "WT" "2D4" "-N" "-N"
diff "WT" "2D4" "-C" "-C"
diff "WT" "WT" "CM" "MM"
diff "WT" "WT" "CM" "-N"
diff "WT" "WT" "CM" "-C"
diff "WT" "WT" "MM" "-N"
diff "WT" "WT" "MM" "-C"
diff "WT" "WT" "-N" "-C"
diff "2D4" "2D4" "CM" "MM"
diff "2D4" "2D4" "CM" "-N"
diff "2D4" "2D4" "CM" "-C"
diff "2D4" "2D4" "MM" "-N"
diff "2D4" "2D4" "MM" "-C"
diff "2D4" "2D4" "-N" "-C"

# differential notpolyA

function notdiff {
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

notdiff "WT" "2D4" "CM" "CM"
notdiff "WT" "2D4" "MM" "MM"
notdiff "WT" "2D4" "-N" "-N"
notdiff "WT" "2D4" "-C" "-C"
notdiff "WT" "WT" "CM" "MM"
notdiff "WT" "WT" "CM" "-N"
notdiff "WT" "WT" "CM" "-C"
notdiff "WT" "WT" "MM" "-N"
notdiff "WT" "WT" "MM" "-C"
notdiff "WT" "WT" "-N" "-C"
notdiff "2D4" "2D4" "CM" "MM"
notdiff "2D4" "2D4" "CM" "-N"
notdiff "2D4" "2D4" "CM" "-C"
notdiff "2D4" "2D4" "MM" "-N"
notdiff "2D4" "2D4" "MM" "-C"
notdiff "2D4" "2D4" "-N" "-C"



# GO enrichment
function go_enrich {
	type=$1
	file="diff_expr/$2_vs_$3_expr.csv"
	echo "" > "${file%%.*}"_go_enrich.csv
	cat "${file%%.*}"_down.csv  "${file%%.*}"_up.csv | cut -f 1 -d "," | sort | uniq > _de_list
	grep -v -f _de_list gene_summary.txt | tail -n +2 | cut -f 1 > _nde_list
	cat _de_list | xargs -ipat grep pat $type | cut -f 2  | sort | uniq >  _go_list
	for go in `cat _go_list`; do grep $go $type | cut -f 1 | grep MGG | sort | uniq > "_"$go; done
	for f in `ls _GO*`
	do
	    a=$(cat $f _de_list | sort | uniq -d | wc -l )
	    b=$(cat $f _nde_list | sort | uniq -d | wc -l)
	    c=$(cat $f $f _de_list | sort | uniq -u | wc -l)
	    d=$(cat $f $f _nde_list | sort | uniq -u | wc -l)
	    echo ${f/_/}"," `Rscript ../../m-oryzae-polya/fisher_test.R $a $b $c $d` >> "${file%%.*}"_go_enrich.csv
	done
	#cat "${file%%.*}"_go_enrich.csv | awk -F "," '{if ($2 <= 0.05) print $0}' | sed 's/,/ /'
	#cat _t | xargs -ipat grep pat $type | cut -f 3  | sort | uniq
	rm _de_list _nde_list _go_list _GO*
}
go_enrich go_terms.csv WT-CM WT-MM
go_enrich go_terms.csv WT-CM WT--N
go_enrich go_terms.csv WT-CM WT--C
go_enrich go_terms.csv WT-MM WT--N
go_enrich go_terms.csv WT-MM WT--C
go_enrich go_terms.csv WT--N WT--C
go_enrich go_terms.csv WT-CM 2D4-CM



# extract regions and nucleotide composition
rm motifs/*fam
for f in `ls *X.polyA_sgl_m`
do
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -100 -36 print "${f%%.*}"_GURICH_sgl_m.fam &
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -35 -10 print "${f%%.*}"_ARICH_sgl_m.fam &
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -15 -3 print "${f%%.*}"_URICH_sgl_m.fam &
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -3 3 print "${f%%.*}"_CUTSITE_sgl_m.fam &
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  4 30 print "${f%%.*}"_LONGURICH_sgl_m.fam 
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -100 100 print "${f%%.*}"_TOT_sgl_m.fam
done
for f in `ls *X.polyA_apa_m`
do
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -100 -36 print "${f%%.*}"_GURICH_apa_m.fam &
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -35 -10 print "${f%%.*}"_ARICH_apa_m.fam &
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -15 -3 print "${f%%.*}"_URICH_apa_m.fam &
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -3 3 print "${f%%.*}"_CUTSITE_apa_m.fam &
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  4 30 print "${f%%.*}"_LONGURICH_apa_m.fam 
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -100 100 print "${f%%.*}"_TOT_apa_m.fam
done
for f in `ls *X.polyA_all_m`
do
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -100 -36 print "${f%%.*}"_GURICH_all_m.fam &
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -35 -10 print "${f%%.*}"_ARICH_all_m.fam &
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -15 -3 print "${f%%.*}"_URICH_all_m.fam &
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -3 3 print "${f%%.*}"_CUTSITE_all_m.fam &
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  4 30 print "${f%%.*}"_LONGURICH_all_m.fam 
    python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $f  -100 100 print "${f%%.*}"_TOT_all_m.fam
done
mv *.fam motifs

# glam alignment
glam2 n WT-CM-X_ARICH_sgl_m.fam -n 40000 -w 6 -O glam2_ARICH_sgl
fimo -oc fimo_ARICH_sgl --norc --thresh 1e-2 glam2_ARICH_sgl/glam2.meme WT-CM-X_TOT_sgl_m.fam 
for i in {1..10}
do
 echo "Motif "$i > "_motif"$i
 cut -f 1,3 fimo_ARICH_sgl/fimo.txt | awk -v motif=$i -v num=`grep -c ">" WT-CM-X_TOT_sgl_m.fam` '{ if ( $1 == motif ) arr[$2]++ } END { for (i=0; i<=200; i++) print arr[i]/num*100 }' >> "_motif"$i
done
paste _motif*
cut -f 9 fimo_ARICH_sgl/fimo.txt | sort | uniq > _t
for i in `cat _t`; do echo -ne $i" "; grep $i WT-CM-X_ARICH_sgl_m.fam -c | awk -v num=`grep -c ">" WT-CM-X_ARICH_sgl_m.fam` '{print $1/num*100}'; done	
	
	
# motif scan
function scan {
    f=$1
    m=$2
    echo "" > _count; for line in `cat $f`; do  echo $line | grep -v ">" | grep -b -o -e $m >> _count; done
    cat _count | awk -v num=`grep -c ">" $f` -F ":" '{ arr[$1]++ } END { for (i=0; i<=200; i++) print arr[i]/num*100 }'
}

 
# extract polyA trascript sequences with fastacmd

cut -f 1,2,3,4 WT_2D4_CM_polyA.diff | awk '{system("fastacmd -d Magnaporthe_oryzae.MG8.18.dna.toplevel.fa -s " "\"lcl|"$2"\" -L "$3","$4" -o blast/"$1".fasta ")}'


# call blast

for f in `ls blast/*fasta`;
do 
blastall -p blastx -m 8 -a 4 -e 1e-20 -i $f -d MG8_proteome.fasta | head -n 1 | cut -f 2 | awk -v f="${f%%.*}" '{system("fastacmd -d MG8_proteome.fasta -s  " "\""$1"\" > "f".hits.fasta  ")}'
done



# extract 3'UTR or intra-APA sequences (for miRNA search)
grep stop_codon Magnaporthe_oryzae.MG8.18.gff3 | sed  -e 's/ID=stop_codon://' -e 's/T.*//' | cut -f 4,9| awk '{print $2,$1}' | sort > _g
sort -k 5,5 -k 2 WT-CM-X.polyA_all_m | awk '{arr[$5"@"$3"@"$4]=arr[$5"@"$3"@"$4]$2":"} END {for(x in arr) print x,arr[x]}  ' | sed -e 's/@/ /g' -e 's/+/2/' -e 's/-/1/'  -e 's/:$//' | sort > _t
#join _g _t | awk '{split($5, arr, ":"); for (x in arr) if (($4 == 1 && arr[x] > $2  && arr[x+1] > $2 && arr[x+1] != "") || ($4 == 2 && arr[x] < $2  && arr[x+1] < $2 && arr[x+1] != "") ) system("fastacmd -d Magnaporthe_oryzae.MG8.18.dna.toplevel.fa -s " "\"lcl|"$3"\" -S " "\""$4"\" -L "arr[x]","arr[x+1]" ")}' > _WT-CM-X.intra
join _g _t | awk '{split($5, arr, ":"); for (x in arr) if ($4 == 1 && arr[x] > $2+3 && arr[x]-($2+3)>8) system("fastacmd -d Magnaporthe_oryzae.MG8.18.dna.toplevel.fa -s " "\"lcl|"$3"\" -S "$4" -L "$2+3","arr[x]" ")}' > _WT-CM-X.intra_sense
join _g _t | awk '{split($5, arr, ":"); for (x in arr) if ($4 == 2 && arr[x] < $2-3 && ($2+3)-arr[x]>8) system("fastacmd -d Magnaporthe_oryzae.MG8.18.dna.toplevel.fa -s " "\"lcl|"$3"\" -S "$4" -L "arr[x]","$2-3" ")}' > _WT-CM-X.intra_antisense


# orphans (400 or 1000 nt) search against known db
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa WT-ALL-X.notpolyA_all_m  -400 0 print WT-ALL-X.notpolyA_all_m_400.fa
blastn -task dc-megablast -query WT-ALL-X.notpolyA_all_m_400.fa -db nt -remote -outfmt 5 | python ../../m-oryzae-polya/parse_blast_xml.py 
blastn -task dc-megablast -query WT-ALL-X.notpolyA_all_m_400.fa -db Rfam.fasta -outfmt 5 | python ../../m-oryzae-polya/parse_blast_xml.py 
# blastn -task dc-megablast -query WT-ALL-X.notpolyA_all_m_400.fa -db fRNAdb_3.0.fasta -outfmt 5  | python ../../m-oryzae-polya/parse_blast_xml.py 
# blastn -task dc-megablast -query WT-ALL-X.notpolyA_all_m_400.fa -db hairpin.fa -outfmt 5 | python ../../m-oryzae-polya/parse_blast_xml.py 
# blastn -task dc-megablast -query WT-ALL-X.notpolyA_all_m_400.fa -db ncrna_NONCODE\[v3.0\].fasta -outfmt 5 | python ../../m-oryzae-polya/parse_blast_xml.py
# perl ../../m-oryzae-polya/rfam_scan.pl --nobig -blastdb Rfam.fasta Rfam.cm WT-ALL-X.notpolyA_all_m_400.fa 

# ophans overlapping annotated genes
python -c "
for line1 in open('WT-ALL-X.notpolyA_all_m', 'r'):
 (val, pos, chrx, sense) = line1.strip().split(' ')
 pos = int(pos)
 for line2 in open('Magnaporthe_oryzae.MG8.18.gff3', 'r'):
  (Chrx, none_1, feat, start, end, none_2,none_3,none_4,infos) = line2.strip().split('\t')
  start = int(start)
  end = int(end)
  if chrx == Chrx and pos >= start and pos <= end and feat == 'gene':
   id_start = infos.index('ID=')
   id_end = infos.index(';', id_start)
   print chrx + ':' + str(pos) + ':' + val + ':' + sense, '\t', infos[id_start + 3: id_end]
"
# orphans differentially expressed
for f in notdiff_polyA/WT*WT*.csv
do
	echo -ne $(basename "${f%%.*}")"," | sed 's/_notpolyA//'
 	cat $f | awk -F "," '{if ($8 < 0.05) print $1,$8}' 
done



# ncRNA search denovo
mugsy --prefix magna --directory /media/marco/Elements/3Tfill/oryzae_18/ncrna/magna/out_magna/ genomes/*
cat out_magna/magna.maf | python ../../../../m-oryzae-polya/maf_order.py Magnaporthe_oryzae > sorted.maf
rnazWindow.pl  --min-seqs=3 sorted.maf > windows.maf
RNAz --both-strands --no-shuffle --cutoff=0.9 windows.maf > rnaz.out
rnazOutputSort.pl rnaz.out | rnazCluster.pl > results.dat
#rnazIndex.pl --gff -w results.dat | sed -e 's/Magnaporthe_oryzae\.//' -e 's/\t\([+\-]\)\t/\t\1\t\.\t/' > results_w.gff
rnazIndex.pl --gff results.dat | sed -e 's/Magnaporthe_oryzae\.//' -e 's/\t.\t/\t.\t.\t/' > results_l.gff
cut -f 1,4,5,7 results_l.gff | awk '{if($4 == "+") S=1; else S=2; system("fastacmd -d ../../Magnaporthe_oryzae.MG8.18.dna.toplevel.fa -s " "\"lcl|"$1"\" -L "$2","$3" -S "S" ")}' > results.fa
blastn -task dc-megablast -query results.fa -db ../../Rfam.fasta -max_target_seqs 1 -outfmt 6 
blastn -task megablast -query results.fa -db ../../WT-CM-X.notpolyA_all_m_400.fa -max_target_seqs 1 -outfmt 6
rnazFilter.pl "z<-3" ../results.dat | grep -e "^locus" | cut -f 2,3,4 | sed 's/Magnaporthe_oryzae\.//' | awk '{system("fastacmd -d ../../../Magnaporthe_oryzae.MG8.18.dna.toplevel.fa -s " "\"lcl|"$1"\" -L "$2","$3" ")}' | RNAfold

# RL-SAGE blast with never expressed genes (sequences downloades from biomart)
blastn -task blastn -query never_expressed.fa -db OSJNGg.fa -outfmt 6 -max_target_seqs 1  | awk '{if ($3 == 100 && $4 == 21) print $0}' | cut -f 1 | xargs -ipat grep pat ../gene_summary.txt


# call interproscan
for f in `ls *.fa`;
do
python ../../../m-oryzae-polya/iprscan_soappy.py --email=marco.marconi@gmail.com --title=marco --sequence=$f  --outfile=$f --outformat=out &
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
   echo $f
   Rscript ../../m-oryzae-polya/correlation.R $f
done

# number of expressed genes (by replicate)
for sample in *.expr
do
    echo -ne "${sample%%.*}"","
    cat $sample | awk '{if ($2 > 0) print $1}' |  wc -l
done
# number of expressed genes (by condition) 1 read in at least 2 replicates
for sample in "2D4-CM" "2D4-MM" "2D4--N" "2D4--C" "WT-CM" "WT-MM" "WT--N" "WT--C" 
do
    echo -ne "${sample%%.*}"","
    cat $sample-1.expr $sample-2.expr $sample-3.expr | awk '{if ($2 > 0) print $1}' | sort | uniq -c | awk '{if ($1 >= 2) print $0}' | wc -l
done
# always expressed genes 
for sample in "WT-CM" "WT-MM" "WT--N" "WT--C" 
do
    cat $sample-1.expr $sample-2.expr $sample-3.expr | awk '{if ($2 > 0) print $1}' | sort | uniq -c | awk '{if ($1 >= 2) print $2}' > "_"$sample"_"t
done
cat _*_t | sort | uniq -c | awk '{if ($1 == 4) print $0}' > _all
cat _all | wc -l
# never expressed genes
cat gene_summary.txt _*_t | tail -n +2 | cut -f 1 | sort | uniq -u | wc -l
rm _*
# genes intersections (or disjunction, for ppt diagrams)
function intersection {
    a=$1
    b=$2
    c=$3
    d=$4
    cat $a-1.expr $a-2.expr $a-3.expr | awk '{if ($2 > 0) print $1}' | sort | uniq -c | awk '{if ($1 >= 2) print $2}' > _a
    cat $b-1.expr $b-2.expr $b-3.expr | awk '{if ($2 > 0) print $1}' | sort | uniq -c | awk '{if ($1 >= 2) print $2}' > _b   
    cat $c-1.expr $c-2.expr $c-3.expr | awk '{if ($2 > 0) print $1}' | sort | uniq -c | awk '{if ($1 >= 2) print $2}' > _c
    cat $d-1.expr $d-2.expr $d-3.expr | awk '{if ($2 > 0) print $1}' | sort | uniq -c | awk '{if ($1 >= 2) print $2}' > _d
    #cat _a | wc -l 
    #cat _b | wc -l
    #cat _c | wc -l
    #cat _a _b _c| sort | uniq -c | awk '{if ($1 == 3) print $0}' | wc -l
    #cat _a  _b | sort | uniq -d | wc -l 
}
intersection WT-CM WT-MM WT--N WT--C



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


# distribution of APA in genes (you have to change the script)
for f in *X.polyA_all_m
do
   echo -ne "${f%%??.*},"	
   python ../../m-oryzae-polya/polyA_distribution.py Magnaporthe_oryzae.MG8.18.gff3 $f 
done

# number of polyA sites per gene 
for f in *X.polyA_all_m
do
    echo -ne "${f%%.*}"","
    cut -f 5 -d " " $f | sort | uniq -c | awk '{count+=1;num+=$1} END {print num/count}'
done

# genes affected by rbp35 in wt (all)
for f in diff_polyA/WT*2D4*down.polyA_all_m
do
 echo -ne $(basename "${f%%.*}")"," | sed 's/_down//'
 cat $f | cut -f 5 -d " " | sort | uniq | wc -l
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

# distribution of polyA sites (you have to change the script)
for f in *X.polyA_all_m
do
   echo -ne "${f%%??.*},"	
   python ../../m-oryzae-polya/polyA_distribution.py Magnaporthe_oryzae.MG8.18.gff3 $f 
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
    orphan=$(cat "${f%%.*}"-X.notpolyA_all_m | wc -l )
    echo $normal","$orphan
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
 cat "${f%%_polyA.*}"_down.polyA_all_m  "${f%%_polyA.*}"_up.polyA_all_m | cut -f 5 -d " " | sort | uniq | wc -l
done
# genes differentially expressed with differentially expressed polyA
for f in diff_polyA/*polyA.csv
do
 echo -ne $(basename "${f%%_polyA.*}")","
 cat "${f%%_polyA.*}"_down.polyA_all_m  "${f%%_polyA.*}"_up.polyA_all_m | cut -f 5 -d " " | sort | uniq | xargs -ipat grep pat _$(basename "${f%%_polyA.*}") | wc -l
done
# G1 sgl and APA
for f in diff_polyA/*down*polyA_all_m
do
 echo -ne $(basename "${f%%.*}")"," | sed 's/_down//'
 cut $f -d " " -f 5 | sort | uniq > _g1
 cut $(basename "${f%%_vs_*}")"-X.polyA_apa_m" -d " " -f 5 | sort | uniq > _apa
 cut $(basename "${f%%_vs_*}")"-X.polyA_sgl_m" -d " " -f 5 | sort | uniq > _sgl
 echo -ne `cat _g1 _sgl | sort | uniq -d | wc -l`","
 echo `cat _g1 _apa | sort | uniq -d | wc -l`
done


# 3' UTR length
for f in *X.polyA_all_m
do
    echo -ne "${f%%??.*}"","
    python ../../m-oryzae-polya/3UTR_length.py Magnaporthe_oryzae.MG8.18.gff3 $f
done 

# 3' UTR length (cumulative only in affected genes)
for f in "CM" "MM" "-N" "-C"
do
    echo "WT vs 2D4 "$f
    cat diff_polyA/WT"-"$f"_"vs"_"2D4"-"$f"_"down.polyA_all_m   | cut -f 5 -d " " | sort | uniq > _diff
    cat WT-$f-X.polyA_all_m | grep -f _diff > _t
    python ../../m-oryzae-polya/3UTR_length.py Magnaporthe_oryzae.MG8.18.gff3 _t
    cat 2D4-$f-X.polyA_all_m | grep -f _diff > _t
    python ../../m-oryzae-polya/3UTR_length.py Magnaporthe_oryzae.MG8.18.gff3 _t    
done

# distance from annotated gene features
# ...


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


fastq-mcf -o 2D4--C-1_1_trimmed.fastq -o 2D4--C-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4--C-11.fastq 2D4--C-12.fastq 
fastq-mcf -o 2D4-CM-1_1_trimmed.fastq -o 2D4-CM-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4-CM-11.fastq 2D4-CM-12.fastq 
fastq-mcf -o 2D4-CM-2_1_trimmed.fastq -o 2D4-CM-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4-CM-21.fastq 2D4-CM-22.fastq 
fastq-mcf -o 2D4-MM-1_1_trimmed.fastq -o 2D4-MM-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4-MM-11.fastq 2D4-MM-12.fastq 
fastq-mcf -o 2D4-MM-2_1_trimmed.fastq -o 2D4-MM-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4-MM-21.fastq 2D4-MM-22.fastq 
fastq-mcf -o 2D4--N-1_1_trimmed.fastq -o 2D4--N-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4--N-11.fastq 2D4--N-12.fastq 
fastq-mcf -o WT--C-1_1_trimmed.fastq -o WT--C-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT--C-11.fastq WT--C-12.fastq 
fastq-mcf -o WT-CM-1_1_trimmed.fastq -o WT-CM-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT-CM-11.fastq WT-CM-12.fastq
fastq-mcf -o WT-CM-2_1_trimmed.fastq -o WT-CM-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT-CM-21.fastq WT-CM-22.fastq 
fastq-mcf -o WT-MM-1_1_trimmed.fastq -o WT-MM-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT-MM-11.fastq WT-MM-12.fastq 
fastq-mcf -o WT-MM-2_1_trimmed.fastq -o WT-MM-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT-MM-21.fastq WT-MM-22.fastq 
fastq-mcf -o WT--N-1_1_trimmed.fastq -o WT--N-1_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT--N-11.fastq WT--N-12.fastq
fastq-mcf -o 2D4--C-2_1_trimmed.fastq -o 2D4--C-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4--C-21.fastq 2D4--C-22.fastq
fastq-mcf -o 2D4--C-3_1_trimmed.fastq -o 2D4--C-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4--C-31.fastq 2D4--C-32.fastq 
fastq-mcf -o 2D4-CM-3_1_trimmed.fastq -o 2D4-CM-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4-CM-31.fastq 2D4-CM-32.fastq 
fastq-mcf -o 2D4-MM-3_1_trimmed.fastq -o 2D4-MM-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4-MM-31.fastq 2D4-MM-32.fastq
fastq-mcf -o 2D4--N-2_1_trimmed.fastq -o 2D4--N-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4--N-21.fastq 2D4--N-22.fastq
fastq-mcf -o 2D4--N-3_1_trimmed.fastq -o 2D4--N-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa 2D4--N-31.fastq 2D4--N-32.fastq 
fastq-mcf -o WT--C-2_1_trimmed.fastq -o WT--C-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT--C-21.fastq WT--C-22.fastq
fastq-mcf -o WT--C-3_1_trimmed.fastq -o WT--C-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT--C-31.fastq WT--C-32.fastq
fastq-mcf -o WT-CM-3_1_trimmed.fastq -o WT-CM-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT-CM-31.fastq WT-CM-32.fastq 
fastq-mcf -o WT-MM-3_1_trimmed.fastq -o WT-MM-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT-MM-31.fastq WT-MM-32.fastq 
fastq-mcf -o WT--N-2_1_trimmed.fastq -o WT--N-2_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT--N-21.fastq WT--N-22.fastq
fastq-mcf -o WT--N-3_1_trimmed.fastq -o WT--N-3_2_trimmed.fastq -0 -l 17 -u ../../m-oryzae-polya/adaptor.fa WT--N-31.fastq WT--N-32.fastq


gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  2D4--C-1_1_trimmed.fastq 2D4--C-1_2_trimmed.fastq > 2D4--C-1.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  2D4-CM-1_1_trimmed.fastq 2D4-CM-1_2_trimmed.fastq > 2D4-CM-1.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  2D4-CM-2_1_trimmed.fastq 2D4-CM-2_2_trimmed.fastq > 2D4-CM-2.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  2D4-MM-1_1_trimmed.fastq 2D4-MM-1_2_trimmed.fastq > 2D4-MM-1.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  2D4-MM-2_1_trimmed.fastq 2D4-MM-2_2_trimmed.fastq > 2D4-MM-2.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  2D4--N-1_1_trimmed.fastq 2D4--N-1_2_trimmed.fastq > 2D4--N-1.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  WT--C-1_1_trimmed.fastq WT--C-1_2_trimmed.fastq > WT--C-1.sam 
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  WT-CM-1_1_trimmed.fastq WT-CM-1_2_trimmed.fastq > WT-CM-1.sam 
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  WT-CM-2_1_trimmed.fastq WT-CM-2_2_trimmed.fastq > WT-CM-2.sam 
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  WT-MM-1_1_trimmed.fastq WT-MM-1_2_trimmed.fastq > WT-MM-1.sam 
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  WT-MM-2_1_trimmed.fastq WT-MM-2_2_trimmed.fastq > WT-MM-2.sam 
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  WT--N-1_1_trimmed.fastq WT--N-1_2_trimmed.fastq > WT--N-1.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  2D4--C-2_1_trimmed.fastq 2D4--C-2_2_trimmed.fastq > 2D4--C-2.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  2D4--C-3_1_trimmed.fastq 2D4--C-3_2_trimmed.fastq > 2D4--C-3.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  2D4-CM-3_1_trimmed.fastq 2D4-CM-3_2_trimmed.fastq > 2D4-CM-3.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  2D4-MM-3_1_trimmed.fastq 2D4-MM-3_2_trimmed.fastq > 2D4-MM-3.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  2D4--N-2_1_trimmed.fastq 2D4--N-2_2_trimmed.fastq > 2D4--N-2.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  2D4--N-3_1_trimmed.fastq 2D4--N-3_2_trimmed.fastq > 2D4--N-3.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  WT--C-2_1_trimmed.fastq WT--C-2_2_trimmed.fastq > WT--C-2.sam 
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  WT--C-3_1_trimmed.fastq WT--C-3_2_trimmed.fastq > WT--C-3.sam 
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  WT-CM-3_1_trimmed.fastq WT-CM-3_2_trimmed.fastq > WT-CM-3.sam 
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  WT-MM-3_1_trimmed.fastq WT-MM-3_2_trimmed.fastq > WT-MM-3.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  WT--N-2_1_trimmed.fastq WT--N-2_2_trimmed.fastq > WT--N-2.sam
gsnap -B 5 -t 8 -A sam -d MG8_18 -D ./MG8_18/  WT--N-3_1_trimmed.fastq WT--N-3_2_trimmed.fastq > WT--N-3.sam

RNAz --both-strands --no-shuffle --cutoff=0.5 maf_parse1.maf > maf_parse1.out &
RNAz --both-strands --no-shuffle --cutoff=0.5 maf_parse2.maf > maf_parse2.out &
RNAz --both-strands --no-shuffle --cutoff=0.5 maf_parse3.maf > maf_parse3.out &
RNAz --both-strands --no-shuffle --cutoff=0.5 maf_parse4.maf > maf_parse4.out &



# altered polyA  in WT vs non altered
cat diff_polyA/WT*2D4*down*polyA* | sort -k 1,7 | uniq > _g
cat _g _g WT-ALL-X.polyA_all_m | sort -k 1,7 | uniq -u | sort -R > _r
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa _g -100 -30 print _s1
c=$(cat _g | wc -l)
head -n $c _r > _r2
tail -n $c _r > _r3
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa _r2 -100 -30   print _s2
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa _r3 -100 -30   print _s3

# altered polyA  in WT  vs same genes in 2D4
cond="CM"
cut -f 5 -d " "  diff_polyA/WT-"$cond"_vs_2D4-"$cond"_down.polyA_all_m | sort | uniq  > _t1
grep -f _t1 2D4-"$cond"-X.polyA_all_m > _g
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa diff_polyA/WT-"$cond"_vs_2D4-"$cond"_down.polyA_all_m -100 100 print _seq1
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa _g -100 100 print _seq2
scan _seq1 "TGTA[TCA]"
scan _seq2 "TGTA[TCA]"

# altered polyA  in 2D4 vs same genes in WT
cond="CM"
cut -f 5 -d " "  diff_polyA/WT-"$cond"_vs_2D4-"$cond"_up.polyA_all_m | sort | uniq  > _t2
grep -f _t2 WT-"$cond"-X.polyA_all_m > _g
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa diff_polyA/WT-"$cond"_vs_2D4-"$cond"_up.polyA_all_m -100 100 print _seq1
python ../../m-oryzae-polya/polyA_nucleotide.py Magnaporthe_oryzae.MG8.18.dna.toplevel.fa _g -100 100 print _seq2
scan _seq1 "TGTA[TCA]"
scan _seq2 "TGTA[TCA]"

# genes with both group, G3 
cat _t1 _t2 | sort | uniq -d | wc -l

