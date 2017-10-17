#!/bin/bash -e

ml python/3.5
source samtools-1.4

PROJDIR=/tgac/workarea/group-tg/projects/Kronos_retilling/ucd_incoming
echo "POSTPROCESSING $1"

target=$(dirname $1)/$(basename $1 .bam).mapped_pairs.bam
if [[ ! -e $target ]]; then
 echo Making $target
 samtools view -bh -f 2 $1 > $target
 echo Making index
 samtools index $target
 echo Calculating MD5
 md5sum $target > $target.md5sum
fi

sample=$(basename $(dirname $1))
src_fq1=$(ls $PROJDIR/data/$sample/*_R1_trimmed.fq.gz)
src_fq2=$(ls $PROJDIR/data/$sample/*_R2_trimmed.fq.gz)

echo $sample
echo $src_fq1
echo $src_fq2

tgt_fq1=$(dirname $src_fq1)/$(basename $src_fq1 .fq.gz).unmapped_iwgsc10.fq.gz
tgt_fq2=$(dirname $src_fq2)/$(basename $src_fq2 .fq.gz).unmapped_iwgsc10.fq.gz

if [[ (! -e $tgt_fq1) || (! -e $tgt_fq2) ]]; then
 echo "Extracting unmapped reads"
 samtools view $target | python3 /tgac/workarea/group-tg/projects/Kronos_retilling/extract_unmapped.py $src_fq1 $src_fq2
 md5sum $tgt_fq1 > $tgt_fq1.md5sum
 md5sum $tgt_fq2 > $tgt_fq2.md5sum
fi

echo "FINISHED"
