#################################################################
#
# date: Ocotober 22, 2016
# platform: Ubuntu 16.04
# R.version : 3.2.2
# author: Villemin Jean-Philippe
# team: Epigenetic Component of Alternative Splicing - IGH
#
# 
# Usage : 
# Rscript /home/jean-philippe.villemin/code/RNA-SEQ/Rscript/SurvivalV2.R -e OS -m /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/outputBasal.tsv -s /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/TCGA_CDR.csv

# Note : Update of V2 because the matrice used with psi doesn't have the same format.
#
#################################################################


####################################################################################
#####################################  Package Loading  ############################
####################################################################################
library(optparse)
library(data.table)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
require(gridExtra)
library(survminer)
library(survival)

####################################################################################
######################### Parameters  ##############################################
####################################################################################

option_list = list(
  make_option(c("-m", "--matricePsi"), type="character", default=NULL, help="Absolute File Input Path for PSI", metavar="character"),
  make_option(c("-s", "--survival"), type="character", default=NULL, help="Absolute File Input Path for Survival dataset", metavar="character"),
  make_option(c("-e", "--endPoint"), type="character", default=NULL, help="OS,DSS,DSI,PSI", metavar="character"),
  make_option(c("-g", "--gene"), type="character", default=NULL, help="geneToPrint", metavar="character"),
  make_option(c("-c", "--color"), type="numeric", default=1, help="Color Trick to set if low is BasalA(blue) or BasalB(red)", metavar="character")
  
)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options
warnings()

time =""

if (opt$endPoint =="OS") {
	
	time = "OS.time"
}
if (opt$endPoint =="DSS") {
	
	time = "DSS.time"
}
if (opt$endPoint =="DFI") {
	
	time = "DFI.time"
}
if (opt$endPoint =="PFI") {
	
	time = "PFI.time"
}

print(time)
print(opt$endPoint)

####################################################################################
#####################################  MAIN  ############################
####################################################################################


survival <- read.table(opt$survival,sep=";", header=TRUE,stringsAsFactors=FALSE, comment.char = "@", na.strings = "#N/A", quote="" )
survival_subset = subset(survival, select = c("bcr_patient_barcode",time,opt$endPoint))
survival_subset_unique_id <- survival_subset[!duplicated(survival_subset$bcr_patient_barcode), ]

# Create Matrice PSI using MolSubtype
type_cancer <- scan(opt$matricePsi, skip = 3 ,sep="\t", nlines = 1, what = character())
id_patient <- scan(opt$matricePsi, skip = 0 ,sep="\t", nlines = 1, what = character())

id_patient <- id_patient[-1] #delete column 1

id_patient<-gsub("Patient: ","",id_patient)
id_patient<-gsub("_.*","",id_patient)
id_patient <- as.data.frame(as.list(id_patient))
id_patient_transposed <-  transpose(id_patient)

type_cancer <- type_cancer[-1] #delete column 1
type_cancer <- gsub("Group2: ","",type_cancer)
type_cancer <- as.data.frame(as.list(type_cancer))
type_cancer_transposed <- transpose(type_cancer)


psi <- read.table(opt$matricePsi,skip = 32 ,sep="\t", header=FALSE,stringsAsFactors=FALSE, comment.char = "@",  quote="" )

endPoint <- opt$endPoint

df <- data.frame(NAME=character(),PVAL=double(),LOW_CUTOFF=double(),HIGH_CUTOFF=double(),LOWCOUNT=integer(),HIGHCOUNT=integer(),HR=double(),PVALHR=double(),CILOW=double(),CIHIGH=double())
for( i in 1:nrow(psi) ){

	if( psi[i,1]==opt$gene) {
		print(i)
		print(psi[i,1])
		
		print("MATCH")
		namefile = psi[i,1]
		
		namefile = paste(opt$endPo,namefile,sep="_")
		print(namefile)
		# Select a specific exon skipping
		res <- subset(psi ,  V1==psi[i,1] )
		
		
		res_psi = subset(res, select = -c(V1) )
		
		res_psi_transposed <- transpose(res_psi)
		
		psi_matrice_clean <- cbind(type_cancer_transposed,res_psi_transposed,id_patient_transposed )
		
		colnames(psi_matrice_clean)[1] = "type_cancer"
		colnames(psi_matrice_clean)[2] = "res_psi"
		colnames(psi_matrice_clean)[3] = "bcr_patient_barcode"
	
		psi_matrice_clean_unique_id <- psi_matrice_clean[!duplicated(psi_matrice_clean$bcr_patient_barcode), ]
	
		final <- merge(x = psi_matrice_clean_unique_id, y = survival_subset_unique_id, by = "bcr_patient_barcode", all.x = TRUE)
	  
		# Remove NA
		final<- final[complete.cases(final[,opt$endPoint]),]
		
		#print(final)
		first_quartile <- quantile(final$res_psi, .30,na.rm=TRUE) #get 1st quartile tertile
		third_quartile <- quantile(final$res_psi, .70,na.rm=TRUE) #get 3rd quartile .70
		
	    # Remove tha
		if (first_quartile==third_quartile){next}
		
		final$group[final$res_psi<=first_quartile]<-"low"
		final$group[final$res_psi>=third_quartile]<-"high"

		
		final_group_annotated <- final[!is.na(final$group),]
	     
		print(endPoint)
		print(time)
	
		print("low")
		
		print(dim(final_group_annotated[final_group_annotated$group=="low",]))
		low_count <- nrow(final_group_annotated[final_group_annotated$group=="low",])
		
		write.csv(final_group_annotated[final_group_annotated$group=="low",],row.names=FALSE,file=paste(namefile,"low.csv",sep="_"))
		
		
		print("high")
		print(dim(final_group_annotated[final_group_annotated$group=="high",]))
		
		high_count <- nrow(final_group_annotated[final_group_annotated$group=="high",])
		
		write.csv(final_group_annotated[final_group_annotated$group=="high",],row.names=FALSE,file=paste(namefile,"high.csv",sep="_"))
		
	    #"#00BFC4", "BASALA"
	    #"BASALB", "#F8766D"
	
		print("total")
		dim(final_group_annotated)
		
	 	print(length(final_group_annotated[,time]))
		if (length(final_group_annotated[,time])<40){next}
	
		surv_object <- Surv(time = as.numeric(final_group_annotated[,time]), event = final_group_annotated[,opt$endPoint] )
		#surv_object
		fit1        <- survfit(surv_object ~ group, data = final_group_annotated)
	
		test = surv_pvalue(fit1, final_group_annotated)
		print(test)
		print(test$pval)
	
		fit.coxph          <- coxph(surv_object ~ group, data = final_group_annotated)
		confidenceInterval <- summary(fit.coxph)$conf.int
		print(summary(fit.coxph))
	
		pvc             <- coef(summary(fit.coxph))[,5]
		hr              <- coef(summary(fit.coxph))[,2]
		low_confidence  <-confidenceInterval[,3]
		high_confidence <-confidenceInterval[,4]
	
		print(low_confidence)	
		print(high_confidence)
		print(hr)
		print(pvc)
		tobepasted= ""
		
		if (opt$color == 1) {
		# WHEN GOES DOWN IN BASAL
		tobepasted = c("HR = ",round(hr,digits=2),"[",round(low_confidence,digits=2),"-",round(high_confidence,digits=2),"]") }
		
		if (opt$color == 2) {
		# WHEN GOES UP IN BASAL
		tobepasted = c("HR = ",round(1/hr,digits=2),"[",round(1/high_confidence,digits=2),"-",round(1/low_confidence,digits=2),"]")
		} 
		
		title <- paste(tobepasted, collapse="",sep="")
		
		df    <- rbind(df, data.frame(NAME = namefile, PVAL = test$pval,LOW_CUTOFF=first_quartile,HIGH_CUTOFF=third_quartile,LOWCOUNT=low_count,HIGHCOUNT=high_count,HR=hr,PVALHR=pvc,CILOW=low_confidence,CIHIGH=high_confidence))
	
		# if 1 : 
		#low is Blue is basalA-Like  "#00BFC4" 
		#high is Red is basalB-like "#F8766D

		# if 2 : 
		#low is Blue is basalA-Like "#00BFC4" 
		#high is Red is basalB-like "#F8766D

		mypalette = ""
		if (opt$color == 1) {
			#high is blue(A), low is red(B)
			mypalette = c("#00BFC4", "#F8766D")
		}
		if (opt$color == 2) {
			#high is red (B), low is blue (A)
			mypalette = c("#F8766D", "#00BFC4")
		}
		
		tobepasted = c(opt$endPoint,"_",opt$color,"_",opt$gene,".png")
		png(paste(tobepasted,sep="/",collapse=""),width=10,height=10,units = 'cm', res = 300)
		ggsurv <- 	ggsurvplot(fit1, data = final_group_annotated, pval = TRUE , legend.title  = "Group :",title=title,break.x.by =365.25,xscale=365.25,xlim=c(0,2000),xlab="Years",palette=mypalette) 

		# use ggpar to change the font of one or a list of ggplots at once

		print({	ggpubr::ggpar(ggsurv,
							legend = "top",
							font.legend = list(size = 12, color = "black", face = "bold"),
							font.x = c(12 , "bold", "black"),           # font for x axises in the plot, the table and censor part
							font.y = c(12, "bold", "black"),       # font for y axises in the plot, the table and censor part
							font.tickslab = c(14, "plain", "black") 
					) })
		dev.off()
		
		# Switch 
		if (opt$color == 1) {
			final_group_annotated$group[final_group_annotated$group=="high"] <- "BasalA-like"
			final_group_annotated$group[final_group_annotated$group=="low"]  <- "BasalB-like"
		}
		if (opt$color == 2) {
			final_group_annotated$group[final_group_annotated$group=="high"] <- "BasalB-like"
			final_group_annotated$group[final_group_annotated$group=="low"]  <- "BasalA-like"
		}
		
		outputname <- paste0(namefile,"_PSI",sep="")
		png(paste0(outputname,".png",sep=""),width=4,height=10,units = 'cm', res = 300)#
		psiPlot   <- ggplot(final_group_annotated, aes(x=group, y=res_psi,color=group )) +
							stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE, label.y = 0.95,label.x = 1.5) + geom_boxplot(outlier.shape=NA) + #label.x = 2.45
							labs(x = "",fill = "",y = "") + ylim(c(0,1)) +
							geom_jitter( position=position_jitter(0.2))+ 
							scale_color_manual(values= c("BasalA-like" = "#00BFC4", "BasalB-like" = "#F8766D"))+ theme(legend.position="none",axis.text.x = element_blank()   ,axis.text.y = element_text(size=12))
		print({ psiPlot }) #+ coord_flip()
		dev.off()
		
		}
}
print("TO FILE")

tobepasted = c(opt$endPoint,"_PVAL_",opt$gene,".csv")
towrite =paste(tobepasted,sep="/",collapse="")
write.csv(df,row.names=TRUE,file=towrite)
