#!/bin/bash

#Set inputs from STDIN
contigpath=$1 #CONTIG output from spades
read1path=$2 #reads.1.fq.gz
read2path=$3 #reads.2.fq.gz
readUpath=$4 #reads.U.fq.gz
outhandle=$5 #unique sample or subsample name
numthreads=$6 #Number of threads

#Copy files to scratch and create run folders
echo
echo "Creating folders and copying files to scratch..."
hostname
echo
mkdir /scratch/mcallis_contigqc
mkdir /scratch/mcallis_contigqc/${outhandle}
mkdir /scratch/mcallis_contigqc/${outhandle}/infiles
cp ${contigpath} /scratch/mcallis_contigqc/${outhandle}/infiles/${outhandle}_contigs_preqc.fasta
cp ${read1path} /scratch/mcallis_contigqc/${outhandle}/infiles/${outhandle}_reads.1.fq.gz
cp ${read2path} /scratch/mcallis_contigqc/${outhandle}/infiles/${outhandle}_reads.2.fq.gz
cp ${readUpath} /scratch/mcallis_contigqc/${outhandle}/infiles/${outhandle}_reads.U.fq.gz

##Simplify contig headers
echo
echo "Creating simplified contig headers..."
echo
cat /scratch/mcallis_contigqc/${outhandle}/infiles/${outhandle}_contigs_preqc.fasta | sed -E "s/^>NODE_([0-9]+)_length_[0-9]+_cov.+/>${outhandle}_NODE_\1/" > /scratch/mcallis_contigqc/${outhandle}/${outhandle}_contigs_simpname.fasta

##Create bowtie db from each assembly
echo
echo "Creating preqc bowtie db..."
echo
bowtie2-build --threads ${numthreads} /scratch/mcallis_contigqc/${outhandle}/${outhandle}_contigs_simpname.fasta /scratch/mcallis_contigqc/${outhandle}/${outhandle}_preqc_db

##Run reads against it's own assembly
echo
echo "Recruiting reads to preqc contigs..."
echo
bowtie2 --no-unal -p ${numthreads} -x /scratch/mcallis_contigqc/${outhandle}/${outhandle}_preqc_db -1 /scratch/mcallis_contigqc/${outhandle}/infiles/${outhandle}_reads.1.fq.gz -2 /scratch/mcallis_contigqc/${outhandle}/infiles/${outhandle}_reads.2.fq.gz -U /scratch/mcallis_contigqc/${outhandle}/infiles/${outhandle}_reads.U.fq.gz -S /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB.SAM 1> /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB_stdout.txt 2> /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB_stderr.txt

##View, sort, and index SAM files to BAM
echo
echo "SAM --> sorted BAM conversion"
echo
samtools view -b -T /scratch/mcallis_contigqc/${outhandle}/${outhandle}_contigs_simpname.fasta /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB.SAM > /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB.bam
rm /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB.SAM
samtools sort /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB.bam -o /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB_sorted.bam
rm /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB.bam
samtools index /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB_sorted.bam /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB_sorted.bam.bai
samtools idxstats /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB_sorted.bam > /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB_idxstats.txt

##Get contig seq lengths
echo
echo "Getting contig sequence lengths..."
echo
seq_lengths /scratch/mcallis_contigqc/${outhandle}/${outhandle}_contigs_simpname.fasta | sed 's/ flag.*len.*    /       /' > /scratch/mcallis_contigqc/${outhandle}/${outhandle}_final.contigs.lengths

##Calculate the non-zero coverage (genomeCoverageBed)
echo
echo "Calculating non-zero coverage..."
echo
genomeCoverageBed -d -split -ibam /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_preqcDB_sorted.bam -g /scratch/mcallis_contigqc/${outhandle}/${outhandle}_final.contigs.lengths | /home/moorer/projects/chimera/contig_cov_from_bed > /scratch/mcallis_contigqc/${outhandle}/${outhandle}_final_contigs.e2e_recruitment.non_zero_cov_perc.txt
#/scratch/mcallis_contigqc/${outhandle}/${outhandle}_final_contigs.e2e_recruitment.non_zero_cov_perc.txt
##Get contig headers with >90% coverage and >2kb (for counting)
echo
echo "Collect contig headers >90% e2e coverage and >2kb..."
echo; echo "2000 90"; awk -v len="2000" -v perc="90" '$2 >= len && $3 >= perc' /scratch/mcallis_contigqc/${outhandle}/${outhandle}_final_contigs.e2e_recruitment.non_zero_cov_perc.txt > /scratch/mcallis_contigqc/${outhandle}/${outhandle}_final_contigs.e2e_recruitment.non_zero_cov_perc.len_at_least_2000.perc_at_least_90.txt && wc -l /scratch/mcallis_contigqc/${outhandle}/${outhandle}_final_contigs.e2e_recruitment.non_zero_cov_perc.len_at_least_2000.perc_at_least_90.txt; echo; echo

##Get contig headers with >90% coverage
echo
echo "Getting contig headers >90% e2e coverage..."
echo
awk '$3 >= 90 { print $2,$1,$3 }' /scratch/mcallis_contigqc/${outhandle}/${outhandle}_final_contigs.e2e_recruitment.non_zero_cov_perc.txt | sort -nr | awk 'BEGIN { OFS="      " } { print $2,$1,$3 }' > /scratch/mcallis_contigqc/${outhandle}/${outhandle}_final_contigs.e2e_recruitment.non_zero_cov_perc.perc_at_least_90.txt

##Pull >90% coverage contigs
echo
echo "Pulling >90% coverage contigs..."
echo
cat /scratch/mcallis_contigqc/${outhandle}/${outhandle}_final_contigs.e2e_recruitment.non_zero_cov_perc.perc_at_least_90.txt | cut -d$'\t' -f1 > /scratch/mcallis_contigqc/${outhandle}/${outhandle}_final_contigs.e2e_recruitment.non_zero_cov_perc.perc_at_least_90_headers.txt
grep_ids /scratch/mcallis_contigqc/${outhandle}/${outhandle}_final_contigs.e2e_recruitment.non_zero_cov_perc.perc_at_least_90_headers.txt /scratch/mcallis_contigqc/${outhandle}/${outhandle}_contigs_simpname.fasta > /scratch/mcallis_contigqc/${outhandle}/${outhandle}_simpname_qccontigs_final.fasta

##Create bowtie db post-qc
echo
echo "Creating bowtie db for post qc contigs..."
echo
bowtie2-build --threads ${numthreads} /scratch/mcallis_contigqc/${outhandle}/${outhandle}_simpname_qccontigs_final.fasta /scratch/mcallis_contigqc/${outhandle}/${outhandle}_postqc_db

##Recruit subsampled reads and all reads against each qc'ed subassembly
echo
echo "Recruit reads to post qc contig database..."
echo
bowtie2 --no-unal -p ${numthreads} -x /scratch/mcallis_contigqc/${outhandle}/${outhandle}_postqc_db -1 /scratch/mcallis_contigqc/${outhandle}/infiles/${outhandle}_reads.1.fq.gz -2 /scratch/mcallis_contigqc/${outhandle}/infiles/${outhandle}_reads.2.fq.gz -U /scratch/mcallis_contigqc/${outhandle}/infiles/${outhandle}_reads.U.fq.gz -S /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB.SAM 1> /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB_stdout.txt 2> /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB_stderr.txt

##View, sort, and index SAM files to BAM
echo
echo "SAM --> sorted BAM conversion..."
echo
samtools view -b -T /scratch/mcallis_contigqc/${outhandle}/${outhandle}_simpname_qccontigs_final.fasta /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB.SAM > /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB.bam
rm /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB.SAM
samtools sort /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB.bam -o /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB_sorted.bam
rm /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB.bam
samtools index /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB_sorted.bam /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB_sorted.bam.bai
samtools idxstats /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB_sorted.bam > /scratch/mcallis_contigqc/${outhandle}/${outhandle}_v_postqcDB_idxstats.txt

echo
echo "DONE CONTIG QC"
echo
