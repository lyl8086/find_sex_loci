#!/bin/bash
bam=$1
ref=$2
out=$3
a=$4
b=$5
SAMTOOLS='samtools'
[ $# -lt 5 ] && echo $0 [bam_list] [ref] [out] [num_sex1] [num_sex2] && exit 1

$SAMTOOLS faidx $ref
cut -f1 $ref.fai | shuf > tmp.ref.lst
# split into 8 chunks.
split -d -n r/8 tmp.ref.lst ref_bed
:>depth.cmd
:>depth.cmd.completed
# save cmd into one file.
for f in ref_bed*; do
    let i=i+1
    sort $f -o $f
    echo "\
:>$i.depth.gz && \
for j in \`cat $f\`; do \
$SAMTOOLS depth -q 20 -Q 20 \
-G UNMAP,SECONDARY,QCFAIL,DUP \
--reference $ref -f $bam \
-d 0 -r \$j | \
python3 ttest.py - $a $b | \
gzip -c >>$i.depth.gz; done" >>depth.cmd
done

ParaFly -c depth.cmd -CPU 8 -vv

cat {1..8}.depth.gz >$out.gz
rm {1..8}.depth.gz
echo "===================================="
echo "Final result: $out.gz"