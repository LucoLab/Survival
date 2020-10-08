#################################################################
#
# date: October 4, 2020
# platform: Ubuntu 16.04
# R.version : 3.2.2
# author: Villemin Jean-Philippe
# team: Epigenetic Component of Alternative Splicing - IGH
#
# 
# Usage : 
# Rscript /home/jean-philippe.villemin/code/RNA-SEQ/Rscript/SurvivalOctobre.R -e OS -m /path2MatricePsiTCGA/outputBasal.tsv -s /path2Survival/TCGA_CDR.csv


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
  make_option(c("-e", "--endPoint"), type="character", default=NULL, help="OS,DSS,DSI,PSI", metavar="character")

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

psi 	 <- read.table(opt$matricePsi,skip = 32 ,sep="\t", header=FALSE,stringsAsFactors=FALSE, comment.char = "@",  quote="" )

endPoint <- opt$endPoint

df 		 <- data.frame(NAME=character(),PVAL=double(),LOW_CUTOFF=double(),HIGH_CUTOFF=double(),LOWCOUNT=integer(),HIGHCOUNT=integer(),HR=double(),PVALHR=double(),CILOW=double(),CIHIGH=double())

for( i in 1:nrow(psi) ){
	print(i)
	print(psi[i,1])

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
	
	print(head(final))
	
	first_quartile <- quantile(final$res_psi, .30,na.rm=TRUE) #get 1st quartile tertile
	third_quartile <- quantile(final$res_psi, .70,na.rm=TRUE) #get 3rd quartile .70
	
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
	file.remove(paste(namefile,"low.csv",sep="_"))
	
	
	print("high")
	print(dim(final_group_annotated[final_group_annotated$group=="high",]))
	
	high_count <- nrow(final_group_annotated[final_group_annotated$group=="high",])
	
	write.csv(final_group_annotated[final_group_annotated$group=="high",],row.names=FALSE,file=paste(namefile,"high.csv",sep="_"))
	file.remove(paste(namefile,"high.csv",sep="_"))
	

	print("final_group_annotated")
	outputname <- paste0(namefile,"_PSI",sep="")
	#png(paste0(outputname,".png",sep=""))#
	#print({ggplot(final_group_annotated, aes(x=final_group_annotated$group, y=final_group_annotated$res_psi,color=final_group_annotated$group )) +
	#stat_compare_means( label.x = 1.5, label.y = 1.05) + geom_boxplot(outlier.shape=NA) +
    #labs(x = "GROUP",fill = "GROUP",y = "PSI") +
	#geom_jitter(shape=16, position=position_jitter(0.2))+ 
	#scale_color_manual(name="GROUP",values=c( "#31a354", "#fc9272","#00FF00"))})
    #"#00BFC4", "BASALA"
    #"BASALB", "#F8766D"
	#dev.off()
	

	if (length(final_group_annotated[,time])<40){next}

	surv_object <- Surv(time = as.numeric(final_group_annotated[,time]), event = final_group_annotated[,opt$endPoint] )
	fit1 		<- survfit(surv_object ~ group, data = final_group_annotated)

	test=surv_pvalue(fit1, final_group_annotated)
	
	print(test)
	print(test$pval)

	fit.coxph 		   <- coxph(surv_object ~ group, data = final_group_annotated)
	confidenceInterval <- summary(fit.coxph)$conf.int
	
	print(summary(fit.coxph))

	pvc 			<- coef(summary(fit.coxph))[,5]
	hr 				<- coef(summary(fit.coxph))[,2]
	low_confidence  <-confidenceInterval[,3]
	high_confidence <-confidenceInterval[,4]

	print(low_confidence)	
	print(high_confidence)
	print(hr)
	print(pvc)
	
	df <- rbind(df, data.frame(NAME = namefile, PVAL = test$pval,LOW_CUTOFF=first_quartile,HIGH_CUTOFF=third_quartile,LOWCOUNT=low_count,HIGHCOUNT=high_count,HR=hr,PVALHR=pvc,CILOW=low_confidence,CIHIGH=high_confidence))

	#png(paste0(namefile,".png",sep=""),width=4,height=4,units = 'cm', res = 300)
	#print({ggsurvplot(fit1, data = final_group_annotated, pval = TRUE,title=namefile,break.x.by =365.25,xscale=365.25,xlim=c(0,2000),xlab="Years")})
	#dev.off()
	
}
print("TO FILE")
write.csv(df,row.names=FALSE,file=paste(opt$endPoint,"PVAL.csv",sep="_"),quote = FALSE)
