
Survival Analysis
=============

https://github.com/LucoLab/Survival

## Overview

You need two files : 

1. TCGA_CDR.csv contains the patients survival information.
2. matrice.tsv contains the psi for TCGA patients.

There is two scripts : 

_**SurvivalOctobre.R**_ create a file for each end points.

_**SurvivalV4**_ plot survival curve for a specific exon.


## prerequisite

_NB_ : matrice.tsv is created on lakitu server using a script that parse all psi files (one for each patient) and filter for a list of exon you want to check.

So here, I give this precomputed matrice for polyA Exons and Ribo0 Exons to check survival in Luminal,Basal TCGA patients.


```shell
#/home/luco/localLib/anaconda3/bin/python3 /home/luco/code/python/prepareDataForHeatmap_withReads.py  -l /home/luco/PROJECT/CCLE/CLEAN/TCGA_Annotation/NewgroupsLumBLumABasal.tsv -d /home/luco/PROJECT/TCGA/TCGA_whippet.10.4/output/ -e /home/luco/PROJECT/EMT_workspace/POLYA/CE.whippet.Garfield.10.4.bed -t /home/luco/PROJECT/CCLE/ALL.POLYA.bed -a CE -f diff
```

## How it works

In each dir you will find two bash script (concatPval.sh). You need to modify the path pointing to R script.

They will run survival analysis using SurvivalOctobre. If you need to plot a specific exon, you need to comment SurvivalV4.R and set exon id.


```shell
	
	Rscript SurvivalOctobre.R -e $i -m $PWD/MatricePsi_ALL.POLYA.bed_diff_NewgroupsLumBLumABasal.tsv -s /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/TCGA_CDR.csv

	#1 : 			#high is blue(A), low is red(B)
	Rscript SurvivalV4.R -g "CD44_chr11:35204511-35204640" -e $i -c 1 -m $PWD/MatricePsi_ALL.POLYA.bed_diff_NewgroupsLumBLumABasal.tsv -s /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/TCGA_CDR.csv
```	

