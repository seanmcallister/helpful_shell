#!/bin/bash
##Convert SAM files to sorted and indexed BAM files. Run from inside SAM file folder (otherwise the naming is going to get weird (and probably wrong).
fastafile=$1
samfile=$2

echo
echo "Converting SAM to BAM..."
echo
samtools view -b -T ${fastafile} ${samfile} > ${samfile}.bam
rm ${samfile}
echo
echo "Sorting BAM file..."
echo
samtools sort ${samfile}.bam -o ${samfile}_sorted.bam
rm ${samfile}.bam
echo
echo "Indexing BAM file..."
echo
samtools index ${samfile}_sorted.bam ${samfile}_sorted.bam.bai
echo
echo "Saving per contig IDX stats..."
echo
samtools idxstats ${samfile}_sorted.bam > ${samfile}_idxstats.txt
