#!/bin/bash
#################################################################################
#Given a folder with :
###1. Prepared bowtie databases for each sample. (i.e. unique_postqc_db.1.bt2)
###2. Fasta file w/ qc contigs for each sample. (i.e. unique_qccontigs.fasta)
###3. Paired and unpaired reads for each sample. (i.e. unique_reads.1/2/U.fq.gz)
###4. File with all the "unique" sample names, one per line.
#Run each pairwise bowtie between all samples.
#Convert SAM to sorted/indexed BAM files.
#################################################################################
source ~/.bash_profile
uniquenamefile=$1
numthreads=$2
self_inclusive_YN=$3
mkdir multisample_recruitment_out
mkdir temp
numlines=`cat ${uniquenamefile} | wc -l`
for i in $(seq 1 ${numlines})
	do for j in $(seq 1 ${numlines})
		do if [ "$i" != "$j" ]
			then
			 unique_namei=`cat ${uniquenamefile} | head -${i} | tail -1`
			 unique_namej=`cat ${uniquenamefile} | head -${j} | tail -1`
			 echo "Running bowtie2 for ${unique_namej} reads against ${unique_namei} contig database..."
			 bowtie2 --no-unal -p ${numthreads} -x ${unique_namei}_postqc_db -1 ${unique_namej}_reads.1.fq.gz -2 ${unique_namej}_reads.2.fq.gz -U ${unique_namej}_reads.U.fq.gz -S temp/${unique_namej}_v_${unique_namei}db.SAM 2> multisample_recruitment_out/${unique_namej}_v_${unique_namei}_stderr.txt
			 cp temp/${unique_namej}_v_${unique_namei}db.SAM multisample_recruitment_out/${unique_namej}_v_${unique_namei}db.SAM
			 cd multisample_recruitment_out
			 sam_to_bam.sh ../${unique_namei}_qccontigs.fasta ${unique_namej}_v_${unique_namei}db.SAM
			 cd ../
		fi
	done
	if [ $self_inclusive_YN = "Y" ]
		then
		 unique_namei=`cat ${uniquenamefile} | head -${i} | tail -1`
		 echo "Running bowtie2 for ${unique_namei} reads against its own contig database..."
		 bowtie2 --no-unal -p ${numthreads} -x ${unique_namei}_postqc_db -1 ${unique_namei}_reads.1.fq.gz -2 ${unique_namei}_reads.2.fq.gz -U ${unique_namei}_reads.U.fq.gz -S temp/${unique_namei}_v_${unique_namei}db.SAM 2> multisample_recruitment_out/${unique_namei}_v_${unique_namei}_stderr.txt
		 cp temp/${unique_namei}_v_${unique_namei}db.SAM multisample_recruitment_out/${unique_namei}_v_${unique_namei}db.SAM
		 cd multisample_recruitment_out
		 sam_to_bam.sh ../${unique_namei}_qccontigs.fasta ${unique_namei}_v_${unique_namei}db.SAM
		 cd ../
	fi
done
