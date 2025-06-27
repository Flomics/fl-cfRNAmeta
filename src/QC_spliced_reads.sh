#!/bin/bash
set -Eeuo pipefail

while getopts "b:g:p:" arg; do
  case $arg in
    b) bam=$OPTARG;;
    g) gtf=$OPTARG;;
    p) prefix=$OPTARG;;

  esac
done

: ${bam?"No bam file provided through -b"}
: ${gtf?"No gtf file provided through -g"}
: ${prefix?"No prefix provided through -p"}

cat $gtf | makeIntrons.pl - > introns.gtf

cat introns.gtf |  perl -lane '$F[11]=~s/"|;//g; $start=$F[3]-1; $end=$F[4]+1; print $F[11]."\t".$F[0]."_".$start."_".$end' | sort -k2,2 > gtf.introns.tsv

echo -e "uniqMappedReads\tsplicedReads\t%splicedReads" > ${prefix}.splicedReads.stats.tsv
# Convert bam to bed12:
#  keep only uniquely mapped reads with -q 255
#  (we want only one mapping per read, otherwise building introns won't work correctly)
samtools view -b -F4 -q 255 $bam | bedtools bamtobed -i stdin -bed12  > $prefix.Aligned.sortedByCoord.bed
# if you need the list of spliced reads, select those records in this BED that have more than one BED block (i.e. exon)
totalUniqMappedReads=$(cat $prefix.Aligned.sortedByCoord.bed | wc -l)
#extract introns from bed12 (need >20GB of RAM)
cat $prefix.Aligned.sortedByCoord.bed | bed12togff | sortgff | makeIntrons.pl - > $prefix.introns.gff

#extract transcript_id (field #1) and intron id (field #2). intron_id is formed of its <chr>_<start>_<end> coordinates, so it identifies introns uniquely
cat $prefix.introns.gff |  perl -lane '$F[11]=~s/"|;//g; $start=$F[3]-1; $end=$F[4]+1; print $F[11]."\t".$F[0]."_".$start."_".$end' | sort -k2,2 > $prefix.introns.tsv
# spliced reads:
splicedReads=$(cut -f1 $prefix.introns.tsv | sort|uniq|wc -l)

echo -e "$totalUniqMappedReads\t$splicedReads" | awk '{print $0"\t"$2/$1*100}' >> ${prefix}.splicedReads.stats.tsv

##################################
### SJ stats, annotated vs novel:
##################################
echo -e "totalSJs\tknownSJs\t%knownSJs" > ${prefix}.spliceJunctions.stats.tsv;
for ds in `tail -n+2 ${prefix}.splicedReads.stats.tsv|cut -f1`; do
echoerr $prefix
# non-redundant list of detected introns:
cut -f2 $prefix.introns.tsv|sort|uniq > tmp1
# non-redundant list of annotated introns:
cut -f2 gtf.introns.tsv |sort|uniq> tmp2
totalSJs=$(cat tmp1|wc -l)
#compare both lists:
comm -1 -2 tmp1 tmp2 | sed 's/_/\t/g' > $prefix.introns.known.txt
comm -2 -3 tmp1 tmp2 | sed 's/_/\t/g' > $prefix.introns.novel.txt
knownSJs=$(cat $prefix.introns.known.txt | wc -l)
echo -e "$totalSJs\t$knownSJs" | awk '{print $0"\t"$2/$1*100}' >> ${prefix}.spliceJunctions.stats.tsv
done
