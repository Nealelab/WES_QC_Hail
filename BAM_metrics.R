
#! /usr/bin/Rscript

## ===== Rscript: BAM_metrics.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Script created: January 2017


## ----- DESCRIPTION:
## Pull out .bam level metrics using the predefined .bam directory path


## ----- REQUIREMENTS:

## R installed
## calling_metadata.txt file in Broad Platform VCF output
## BAM directory paths with default file structure in directory


## ----- Example usage:
## Rscript BAM_metrics.R [working directory] [input] [output]

## working directory: '/path/to/working/directory'
## input: 'VCF_filename.calling_metadata.txt'
## output: 'BAM_metrics.txt'

## ALS example
# Rscript BAM_metrics.R \
# /psych/genetics_data/sfarhan/ALSGEN_2017 \
# ALS_Case_Control_Exomes.calling_metadata.txt \
# ALS_hybrid_capture.txt & \

## UK_Ireland example
# Rscript BAM_metrics.R /psych/genetics_data/howrigan/projects/uk_ireland /seq/dax/C2003/Exome/v1/C2003.calling_metadata.txt uk_ireland_BAM_metrics.txt &





## read in command line arguments
args <- commandArgs(TRUE)
wdir <- args[1]
input <- args[2]
output <- args[3]

## set working directory
setwd(wdir)

## read in BAM paths

meta <- read.table(input,h=T,sep='\t')

## --- BAM level variables to add :

## from .hybrid_selection_metrics
BAIT_SET <- NA
TARGET_TERRITORY <- NA
PF_UQ_READS_ALIGNED <- NA
PF_UQ_BASES_ALIGNED <- NA
ON_TARGET_BASES <- NA
PCT_SELECTED_BASES <- NA
MEAN_TARGET_COVERAGE <- NA
ZERO_CVG_TARGETS_PCT <- NA
PCT_TARGET_BASES_2X <- NA
PCT_TARGET_BASES_10X <- NA
PCT_TARGET_BASES_20X <- NA
PCT_TARGET_BASES_30X <- NA

## from .selfSM
FREEMIX_CONTAMINATION <- NA

## from .alignment_summary_metrics
PCT_CHIMERAS <- NA
STRAND_BALANCE <- NA

## strip off the '.bam' trailer
tmp <- substr(as.character(meta$bam),start=1,stop=nchar(as.character(meta$bam))-4)


## LooP to get metrics:

for (i in 1:nrow(meta)) {

if( file.exists(paste0(tmp[i],'.hybrid_selection_metrics')) ) {

hybrid <- read.table(paste0(tmp[i],'.hybrid_selection_metrics'),h=T,skip=6,fill=T)
BAIT_SET[i] <- as.character(hybrid$BAIT_SET)[1]
TARGET_TERRITORY[i] <- hybrid$TARGET_TERRITORY[1]
PF_UQ_READS_ALIGNED[i] <- hybrid$PF_UQ_READS_ALIGNED[1]
PF_UQ_BASES_ALIGNED[i] <- hybrid$PF_UQ_BASES_ALIGNED[1]
ON_TARGET_BASES[i] <- hybrid$ON_TARGET_BASES[1]
PCT_SELECTED_BASES[i] <- hybrid$PCT_SELECTED_BASES[1]
MEAN_TARGET_COVERAGE[i] <- hybrid$MEAN_TARGET_COVERAGE[1]
ZERO_CVG_TARGETS_PCT[i] <- hybrid$ZERO_CVG_TARGETS_PCT[1]
PCT_TARGET_BASES_2X[i] <- hybrid$PCT_TARGET_BASES_2X[1]
PCT_TARGET_BASES_10X[i] <- hybrid$PCT_TARGET_BASES_10X[1]
PCT_TARGET_BASES_20X[i] <- hybrid$PCT_TARGET_BASES_20X[1]
PCT_TARGET_BASES_30X[i] <- hybrid$PCT_TARGET_BASES_30X[1]

selfSM <- read.table(paste(tmp[i],'.selfSM',sep=''),h=T,comment.char='',fill=T)
FREEMIX_CONTAMINATION[i] <- selfSM$FREEMIX

align <- read.table(paste(tmp[i],'.alignment_summary_metrics',sep=''),h=T,skip=6,fill=T)
PCT_CHIMERAS[i] <- align$PCT_CHIMERAS[1]
STRAND_BALANCE[i] <- mean(align$STRAND_BALANCE)

} else {

BAIT_SET[i] <- NA
TARGET_TERRITORY[i] <- NA
PF_UQ_READS_ALIGNED[i] <- NA
PF_UQ_BASES_ALIGNED[i] <- NA
ON_TARGET_BASES[i] <- NA
PCT_SELECTED_BASES[i] <- NA
MEAN_TARGET_COVERAGE[i] <- NA
ZERO_CVG_TARGETS_PCT[i] <- NA
PCT_TARGET_BASES_2X[i] <- NA
PCT_TARGET_BASES_10X[i] <- NA
PCT_TARGET_BASES_20X[i] <- NA
PCT_TARGET_BASES_30X[i] <- NA

FREEMIX_CONTAMINATION[i] <- NA

PCT_CHIMERAS[i] <- NA
STRAND_BALANCE[i] <- NA

	} ## end IF / ELSE

print(i)

} ## END of i LooP


## ---- Add these metrics to the file
meta$BAIT_SET <- BAIT_SET
meta$TARGET_TERRITORY <- TARGET_TERRITORY
meta$PF_UQ_READS_ALIGNED <- PF_UQ_READS_ALIGNED
meta$PF_UQ_BASES_ALIGNED <- PF_UQ_BASES_ALIGNED
meta$ON_TARGET_BASES <- ON_TARGET_BASES
meta$PCT_SELECTED_BASES <- PCT_SELECTED_BASES
meta$MEAN_TARGET_COVERAGE <- MEAN_TARGET_COVERAGE
meta$ZERO_CVG_TARGETS_PCT <- ZERO_CVG_TARGETS_PCT
meta$PCT_TARGET_BASES_2X <- PCT_TARGET_BASES_2X
meta$PCT_TARGET_BASES_10X <- PCT_TARGET_BASES_10X
meta$PCT_TARGET_BASES_20X <- PCT_TARGET_BASES_20X
meta$PCT_TARGET_BASES_30X <- PCT_TARGET_BASES_30X
meta$FREEMIX_CONTAMINATION <- FREEMIX_CONTAMINATION
meta$PCT_CHIMERAS <- PCT_CHIMERAS
meta$STRAND_BALANCE <- STRAND_BALANCE


## write to file
write.table(meta,output,col=T,row=F,quo=F,sep=',')


## ==== END of R script
