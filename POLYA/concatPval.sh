#!/usr/bin/env bash

set -o errexit
set -o pipefail
set -o nounset
#set -o xtrace

echo $PWD 
exp=( DSS OS DFI PFI )
	
touch ALL_PVAL.csv
for i in "${exp[@]}"
do
	echo $i ;
	survivalType="$i"
	
	#Rscript /home/jean-philippe.villemin/code/RNA-SEQ/Rscript/SurvivalOctobre.R -e $i -m $PWD/MatricePsi_ALL.POLYA.bed_diff_NewgroupsLumBLumABasal.tsv -s /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/TCGA_CDR.csv

	#1 : 			#high is blue(A), low is red(B)
	#Rscript /home/jean-philippe.villemin/code/RNA-SEQ/Rscript/SurvivalV4.R -g "CD44_chr11:35204511-35204640" -e $i -c 1 -m $PWD/MatricePsi_ALL.POLYA.bed_diff_NewgroupsLumBLumABasal.tsv -s /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/TCGA_CDR.csv
	
done

## MEGA BUG
#/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  { if( $4 < 0.05 &&  sqrt(($6-$5)*($6-$5)) >= sqrt(0.1*0.1) && $7 >= 40 && $8 >= 40 ) print $_;}' /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/ALL_PVAL.csv > /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/ALL_PVAL.FILTRED.csv

