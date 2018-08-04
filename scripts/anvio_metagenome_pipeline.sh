#!/bin/bash
##Set from STDIN the locations of necessary files to copy
contigpath=$1 #476_BS1_simpname_qccontigs_final.fasta
bowtiebampath=$2 #476_BS1_v_postqcDB_sorted.bam
dastool_scaffolds2bin=$3 #476_BS1_DAStool_DASTool_scaffolds2bin.txt
renamedbinspath=$4 #folder with renamed bins
uniqsample=$5 #unique sample name
numthreads=$6
mincontiglen=$7 #Contig length for binning collection
outpath=$8

#Anvio RUN:
source ~/.bash_profile
mkdir ${outpath}
cd ${outpath}
source activate anvio3
sleep 10s
#Skipped the contig formatting step because they already have simplified headers. And I don't want to remap.
echo
echo "Making contigs database..."
echo
anvi-gen-contigs-database -f ${contigpath} -o ${uniqsample}_contigs.db -n "Contigs database for ${uniqsample}" #Creat contigs database
sleep 10s
echo
echo "Running SCG HMMs for completeness and contamination estimates..."
echo
anvi-run-hmms -c ${uniqsample}_contigs.db -T 1 #Do the default anvio HMM search for single copy genes (can also run a second time to identify custom HMMs)
#anvi-display-contigs-stats 476_BS1_anvio.db #this not on Biomix (generates an html)
#anvi-setup-ncbi-cogs #Can do this manually when desired. Sets up and updates the NCBI cog database (required an affirmative ENTER command)
sleep 10s
echo
echo "Running COG prediction..."
echo
anvi-run-ncbi-cogs -c ${uniqsample}_contigs.db -T ${numthreads}
sleep 10s
echo
echo "Exporting gene amino acid sequences.."
echo
anvi-get-aa-sequences-for-gene-calls -c ${uniqsample}_contigs.db -o ${uniqsample}_protein_sequences.fa #pull amino acid sequences for annotation
#interproscan doesn't seem to work straight off on Biomix; abandon annotation?
echo
echo "Exporting gene nucleotide sequences..."
echo
anvi-get-dna-sequences-for-gene-calls -c ${uniqsample}_contigs.db -o ${uniqsample}_DNA_gene_calls.fa #for taxonomy prep
sleep 10s
echo
echo "Running CENTRIFUGE taxonomy prediction..."
echo
centrifuge -f -x $CENTRIFUGE_BASE/p+h+v/p+h+v ${uniqsample}_DNA_gene_calls.fa -S ${uniqsample}_centrifuge_hits.tsv #taxonomy calls
anvi-import-taxonomy -c ${uniqsample}_contigs.db -i centrifuge_report.tsv ${uniqsample}_centrifuge_hits.tsv -p centrifuge #import centrifuge taxonomy calls
sleep 10s
echo
echo "Profiling contigs database with bam..."
echo
anvi-profile -i ${bowtiebampath} -c ${uniqsample}_contigs.db -o profile_out --sample-name S_${uniqsample} --min-contig-length ${mincontiglen} --cluster-contigs -T ${numthreads}
sleep 10s

#Also I wanted to try what Rika did (mapping other similar samples to one library and merging that
#anvi-profile -i ../../../qc_contigs/479_BS3/multisample_recruitment/483_BS63_v_postqc479BS3db.SAM.bam_sorted.bam -c ../479_BS3_contigs.db -o ./483_BS63 --sample-name S_483_BS63_map --min-contig-length 2000 -T 24
#anvi-merge ../profile_out/PROFILE.db 476_BS1/PROFILE.db 479_BS4/PROFILE.db 481_BS4/PROFILE.db 482_BS8/PROFILE.db 483_BS63/PROFILE.db -o merged_profiles -c ../479_BS3_contigs.db

echo
echo "Importing DASTool Bins and Zeta Bins Collections..."
anvi-import-collection ${dastool_scaffolds2bin} -p profile_out/PROFILE.db -c ${uniqsample}_contigs.db --contigs-mode --collection-name DASTool_bins
sleep 10s
#cd ${renamedbinspath}
#for f in *Zeta*; do cat $f | grep ">" | sed -E "s/>/${f}\t/" | sed -E 's/\.fa//' | awk 'BEGIN{FS="\t";OFS="\t";} {print $2,$1}' >> tabdelim_zeta_bins.txt; done # create Zeta bin only collection
#cd ${outpath}
anvi-import-collection ${renamedbinspath}/tabdelim_zeta_bins.txt -p profile_out/PROFILE.db -c ${uniqsample}_contigs.db --contigs-mode --collection-name Zeta_bins
sleep 10s
#TO Do Next: Transfer files to mac and run these
#anvi-interactive -p profile_out/PROFILE.db -c 476_BS1_anvio.db

#anvi-summarize -p profile_out/PROFILE.db -c 476_BS1_anvio.db -o Bin-Summary -C DASTool_bins
#anvi-summarize -p profile_out/PROFILE.db -c 476_BS1_anvio.db -o Bin-Summary -C Zeta_bins

#anvi-refine -p profile_out/PROFILE.db -c 476_BS1_anvio.db -C DASTool_bins -b BIN-OF-INTEREST

echo
echo "Done"
echo
