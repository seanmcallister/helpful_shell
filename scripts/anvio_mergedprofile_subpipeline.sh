#!/bin/bash
bampath=$1
dastool_scaffolds2bin=$2
uniqsample=$3
renamedbinspath=$4
outpath=$5
mincontiglen=$6
numthreads=$7
bamsuffix=$8

##Goal to create an additional merged profiles from multisample mappings.
source ~/.bash_profile
source activate anvio3
sleep 10s
mkdir ${outpath}/multiprofile_run
cd ${bampath}
for f in *bam
	do base=$(basename $f ${bamsuffix})
	anvi-profile -i ${f} -c ${outpath}/${uniqsample}_contigs.db -o ${outpath}/multiprofile_run/${f} --sample-name S_${base} --min-contig-length ${mincontiglen} -T ${numthreads}
	echo "multiprofile_run/${f}/PROFILE.db" >> ${outpath}/multiprofile_run/sample_bam_names.txt
	sleep 5s
done
sleep 5s
cd ${outpath}
bamstring=`awk '$1=$1' RS="\n" OFS=" " ./multiprofile_run/sample_bam_names.txt`
anvi-merge profile_out/PROFILE.db ${bamstring} -o multiprofile_run/merged_profiles -c ${uniqsample}_contigs.db --enforce-hierarchical-clustering
sleep 10s

##IMPORT COLLECTIONS
anvi-import-collection ${dastool_scaffolds2bin} -p multiprofile_run/merged_profiles/PROFILE.db -c ${uniqsample}_contigs.db --contigs-mode --collection-name DASTool_bins
sleep 5s
anvi-import-collection ${renamedbinspath}/tabdelim_zeta_bins.txt -p multiprofile_run/merged_profiles/PROFILE.db -c ${uniqsample}_contigs.db --contigs-mode --collection-name Zeta_bins

