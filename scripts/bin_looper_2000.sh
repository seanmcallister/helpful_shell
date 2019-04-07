#!/bin/bash
##Set from STDIN the locations of necessary files to copy
contigpath=$1 #qc'ed contig fasta file with simple headers
bowtiebampath=$2 #sorted bam file with read recruitment to qc'ed contigs
bowtiebambaipath=$3 #bam.bai file to go with the sorted bam
read1path=$4 #reads.1.fq
read2path=$5 #reads.2.fq
readUpath=$6 #reads.U.fq
outhandle=$7 #Prefix for outhandles
numthreads=$8 #Number of threads
readLENGTH=$9 #Read length in bp

##Copy necessary files to scratch and create folders
echo
echo "Creating folders and copying files to scratch..."
echo
mkdir /scratch/mcallis_binlooper
mkdir /scratch/mcallis_binlooper/${outhandle}
mkdir /scratch/mcallis_binlooper/${outhandle}/infiles
mkdir /scratch/mcallis_binlooper/${outhandle}/binsanity
mkdir /scratch/mcallis_binlooper/${outhandle}/concoct
mkdir /scratch/mcallis_binlooper/${outhandle}/DASTool
mkdir /scratch/mcallis_binlooper/${outhandle}/maxbin
mkdir /scratch/mcallis_binlooper/${outhandle}/metabat
mkdir /scratch/mcallis_binlooper/${outhandle}/metabat/superspecific
mkdir /scratch/mcallis_binlooper/${outhandle}/metabat/verysensitive
cp ${contigpath} /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta
cp ${bowtiebampath} /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_sorted.bam
cp ${bowtiebambaipath} /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_sorted.bam.bai
cp ${read1path} /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_reads.1.fq.gz
cp ${read2path} /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_reads.2.fq.gz
cp ${readUpath} /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_reads.U.fq.gz
echo
echo "Unzipping read files"
echo
unpigz -p ${numthreads} /scratch/mcallis_binlooper/${outhandle}/infiles/*.fq.gz

##Run metabat superspecific and very sensitive
echo
echo "Running metabat superspecific..."
echo
cd /scratch/mcallis_binlooper/${outhandle}/metabat/superspecific
runMetaBat.sh --superspecific -t ${numthreads} --unbinned --keep -B 30 /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_sorted.bam
mkdir /scratch/mcallis_binlooper/${outhandle}/metabat/superspecific/bins
mv /scratch/mcallis_binlooper/${outhandle}/metabat/superspecific/*.fa /scratch/mcallis_binlooper/${outhandle}/metabat/superspecific/bins/
mv /scratch/mcallis_binlooper/${outhandle}/metabat/superspecific/bins/*unbinned.fa /scratch/mcallis_binlooper/${outhandle}/metabat/superspecific/
echo
echo "Running metabat verysensitive..."
echo
cd /scratch/mcallis_binlooper/${outhandle}/metabat/verysensitive
runMetaBat.sh --verysensitive -t ${numthreads} --unbinned --keep -B 30 /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_sorted.bam
mkdir /scratch/mcallis_binlooper/${outhandle}/metabat/verysensitive/bins
mv /scratch/mcallis_binlooper/${outhandle}/metabat/verysensitive/*.fa /scratch/mcallis_binlooper/${outhandle}/metabat/verysensitive/bins/
mv /scratch/mcallis_binlooper/${outhandle}/metabat/verysensitive/bins/*unbinned.fa /scratch/mcallis_binlooper/${outhandle}/metabat/verysensitive/

##Run maxbin
echo
echo "Running maxbin..."
echo
cd /scratch/mcallis_binlooper/${outhandle}/maxbin
/home/mcallis/software/idba-1.1.1/bin/fq2fa --filter --merge /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_reads.1.fq /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_reads.2.fq /scratch/mcallis_binlooper/${outhandle}/maxbin/reads.1_and_2_interleaved.fa
/home/mcallis/software/idba-1.1.1/bin/fq2fa --filter /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_reads.U.fq /scratch/mcallis_binlooper/${outhandle}/maxbin/reads.unpaired.fa
~/software/MaxBin-2.2/run_MaxBin.pl -contig /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta -reads /scratch/mcallis_binlooper/${outhandle}/maxbin/reads.1_and_2_interleaved.fa -reads2 /scratch/mcallis_binlooper/${outhandle}/maxbin/reads.unpaired.fa -min_contig_length 2000 -thread ${numthreads} -out ${outhandle}_maxbin_out
mkdir /scratch/mcallis_binlooper/${outhandle}/maxbin/bins
mv /scratch/mcallis_binlooper/${outhandle}/maxbin/*.fasta /scratch/mcallis_binlooper/${outhandle}/maxbin/bins/

##Run concoct
echo
echo "Running concoct..."
echo
cd /scratch/mcallis_binlooper/${outhandle}/concoct
cat /scratch/mcallis_binlooper/${outhandle}/metabat/superspecific/${outhandle}_qccontig.fasta.depth.txt | cut -d$'\t' -f1,3 > /scratch/mcallis_binlooper/${outhandle}/concoct/${outhandle}_contigdepth.simple.txt
sleep 10s
concoct --coverage_file /scratch/mcallis_binlooper/${outhandle}/concoct/${outhandle}_contigdepth.simple.txt --composition_file /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta -l 2000 -r ${readLENGTH} -b /scratch/mcallis_binlooper/${outhandle}/concoct/${outhandle}.concoct.out -s 0
mkdir /scratch/mcallis_binlooper/${outhandle}/concoct/variance
mv /scratch/mcallis_binlooper/${outhandle}/concoct/*variances* /scratch/mcallis_binlooper/${outhandle}/concoct/variance/
concoct_pull_contigs.pl -f /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta -c /scratch/mcallis_binlooper/${outhandle}/concoct/${outhandle}.concoct.out_clustering_gt2000.csv -o ${outhandle}
mkdir /scratch/mcallis_binlooper/${outhandle}/concoct/bins
mv /scratch/mcallis_binlooper/${outhandle}/concoct/*.fasta /scratch/mcallis_binlooper/${outhandle}/concoct/bins/
sleep 10s
echo
echo "CONCOCT finished..."
echo

##Run Binsanity
echo
echo "Starting Binsanity binning..."
echo
cd /scratch/mcallis_binlooper/${outhandle}/binsanity
grep ">" /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta | cut -d ">" -f2 > /scratch/mcallis_binlooper/${outhandle}/binsanity/${outhandle}_contigids.txt
Binsanity-profile -i /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta -s /scratch/mcallis_binlooper/${outhandle}/infiles/ --ids /scratch/mcallis_binlooper/${outhandle}/binsanity/${outhandle}_contigids.txt -c ${outhandle}_binsanity --transform scale
cp /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta /scratch/mcallis_binlooper/${outhandle}/binsanity/
Binsanity-wf -c ${outhandle}_binsanity.cov.x100.lognorm -f ./ -l ${outhandle}_qccontig.fasta -x 2000 --threads 2
##NEED TO ADD IN AN IF STATEMENT HERE TO FIX THE BINSANITY OUTPUT IF THE low_completion.fna refinement step fails (which it does if there aren't any contigs in that file). Will require making final binsanity folders and copying the final complete and refined bins as well as results files to the proper folder.
cp -r /scratch/mcallis_binlooper/${outhandle}/binsanity/BinSanity-Final-bins /scratch/mcallis_binlooper/${outhandle}/binsanity/BinSanity-Final-bins_renamed
cd /scratch/mcallis_binlooper/${outhandle}/binsanity/BinSanity-Final-bins_renamed
binsanity_cleanup.sh ${outhandle}_qccontig-bin_


##Convert outputs for DASTool
echo
echo "Converting for DASTool input..."
echo
mkdir /scratch/mcallis_binlooper/${outhandle}/binsanity/DASTool_in
mkdir /scratch/mcallis_binlooper/${outhandle}/concoct/DASTool_in
mkdir /scratch/mcallis_binlooper/${outhandle}/maxbin/DASTool_in
mkdir /scratch/mcallis_binlooper/${outhandle}/metabat/superspecific/DASTool_in
mkdir /scratch/mcallis_binlooper/${outhandle}/metabat/verysensitive/DASTool_in
cd /scratch/mcallis_binlooper/${outhandle}/binsanity/DASTool_in
DASTool_converter.pl -f /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta -m binsanity -b /scratch/mcallis_binlooper/${outhandle}/binsanity/BinSanity-Final-bins_renamed -p ${outhandle}_qccontig-bin_renamed_ -s .fna -o ${outhandle}
cd /scratch/mcallis_binlooper/${outhandle}/concoct/DASTool_in
DASTool_converter.pl -f /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta -m concoct -b /scratch/mcallis_binlooper/${outhandle}/concoct/bins -p ${outhandle}_bin_ -s .fasta -o ${outhandle}
cd /scratch/mcallis_binlooper/${outhandle}/maxbin/DASTool_in
DASTool_converter.pl -f /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta -m maxbin -b /scratch/mcallis_binlooper/${outhandle}/maxbin/bins -p ${outhandle}_maxbin_out. -s .fasta -o ${outhandle}
cd /scratch/mcallis_binlooper/${outhandle}/metabat/superspecific/DASTool_in
DASTool_converter.pl -f /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta -m metabat_super -b /scratch/mcallis_binlooper/${outhandle}/metabat/superspecific/bins -p ${outhandle}_qccontig.fasta.metabat-bins-_--superspecific_-t_${numthreads}_--unbinned_--keep_-B_30. -s .fa -o ${outhandle}
cd /scratch/mcallis_binlooper/${outhandle}/metabat/verysensitive/DASTool_in
DASTool_converter.pl -f /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta -m metabat_very -b /scratch/mcallis_binlooper/${outhandle}/metabat/verysensitive/bins -p ${outhandle}_qccontig.fasta.metabat-bins-_--verysensitive_-t_${numthreads}_--unbinned_--keep_-B_30. -s .fa -o ${outhandle}


##Run DASTool
echo
echo "Running DASTool..."
echo
cd /scratch/mcallis_binlooper/${outhandle}/DASTool
~/software/DAS_Tool/DAS_Tool.sh -i /scratch/mcallis_binlooper/${outhandle}/binsanity/DASTool_in/${outhandle}_binsanity_scaffold2bin.tsv,/scratch/mcallis_binlooper/${outhandle}/concoct/DASTool_in/${outhandle}_concoct_scaffold2bin.tsv,/scratch/mcallis_binlooper/${outhandle}/maxbin/DASTool_in/${outhandle}_maxbin_scaffold2bin.tsv,/scratch/mcallis_binlooper/${outhandle}/metabat/superspecific/DASTool_in/${outhandle}_metabat_super_scaffold2bin.tsv,/scratch/mcallis_binlooper/${outhandle}/metabat/verysensitive/DASTool_in/${outhandle}_metabat_very_scaffold2bin.tsv -c /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta -o ${outhandle}_DAStool -l binsanity,concoct,maxbin,metabat_super,metabat_very --write_bins 1 -t ${numthreads}
cp -r /scratch/mcallis_binlooper/${outhandle}/DASTool/${outhandle}_DAStool_DASTool_bins /scratch/mcallis_binlooper/${outhandle}/DASTool/bins

##Run Checkm DASTool
echo
echo "Running checkm on DASTool output..."
echo
cd /scratch/mcallis_binlooper/${outhandle}/DASTool
checkm lineage_wf -x fa -t 2 --pplacer_threads 2 ./bins ./checkm_out 1> checkm_out_stdout.txt 2> checkm_out_stderr.txt

##Run Checkm 16S finder
echo
echo "Running checkm 16Sfinder on DASTool output..."
echo
cd /scratch/mcallis_binlooper/${outhandle}/DASTool
checkm ssu_finder -x fa -t 10 /scratch/mcallis_binlooper/${outhandle}/infiles/${outhandle}_qccontig.fasta bins checkm_16Sfinder

##Run phylosift. Not in this round.
#for f in 667_BS4_metabat_very_scaffold2bin.031.contigs.fa; do phylosift all --output=${f}_out --updated --threads=10 $f; done

##Run bin_coverage_individualassembly.pl
echo
echo "Calculating bin coverage for abundance estimates..."
echo
cd /scratch/mcallis_binlooper/${outhandle}/DASTool
bin_coverage_individualassembly.pl -b /scratch/mcallis_binlooper/${outhandle}/DASTool/bins -c /scratch/mcallis_binlooper/${outhandle}/metabat/superspecific/${outhandle}_qccontig.fasta.depth.txt -o ${outhandle}_coverage -x fa

##Run prodigal for ORF calling and create ORF sets per DASTool bin
echo
echo "Predicting ORFs with prodigal..."
echo
mkdir /scratch/mcallis_binlooper/${outhandle}/DASTool/ORFs
cd /scratch/mcallis_binlooper/${outhandle}/DASTool/bins
for f in *.fa
	do prodigal -a /scratch/mcallis_binlooper/${outhandle}/DASTool/ORFs/${f}_ORFs.faa -d /scratch/mcallis_binlooper/${outhandle}/DASTool/ORFs/${f}_ORFs.fna -i $f -o /scratch/mcallis_binlooper/${outhandle}/DASTool/ORFs/${f}_out.gbk -s /scratch/mcallis_binlooper/${outhandle}/DASTool/ORFs/${f}_ORFs_scores.txt
done

echo
echo "BIN LOOPAGE COMPLETED"
echo
