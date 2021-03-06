#!/bin/bash
#SBATCH --job-name=RaxML_zetas
#SBATCH --ntasks=22
#SBATCH --nodelist=biomix16
#SBATCH --mem=32000
raxml-pthreads -x 2421 -f a -m PROTCATJTT -n 02_ryan_tree_out.txt -p 748 -s /home/mcallis/jobs/trees/02_ryan_tree_R2.phy -T 22 -w /home/mcallis/jobs/trees/RAXML_WORKING2 -N 300

#-x passes the random interger for fast bootstrapping; if you want to run regular bootstraps, pass -b
#-m PROTCATJTT. CAT is good for large datasets. JTT has specific models for membrane proteins.
#-f a = rapid bootstrap AND search for best-scoring ML tree

#Change:
#-n = name of output file
#change the working directory (-d in fourth line)
#-o = name of outgroup
#-s=name of input file
#-T=number of processors used; change number to match in second line
#-w=working directory where everything will be stored
