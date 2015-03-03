#!$Rbin

####
# consolidateAndAlign.R is a worker script that runs on the cluster nodes.
# It first consolidates the sorted reads from various demuxer, then aligns
# them to the template
# 
# Written by Jochen Weile <jochenweile@gmail.com> and Anjali Gopal <anjali.gopal91@gmail.com>

source("lib/liblogging.R")   #Logger
source("lib/cliargs.R")      #Command-line argument processing
source("lib/libyogitools.R") #Helper functions

# Working directory
dir.name <- getArg("dir",required=TRUE)
if (regexpr("/$",dir.name) < 1) {
	dir.name <- paste(dir.name,"/",sep="")
}
# Job ID
job.id <- getArg("id",required=TRUE)
# ORF Sequence DB
orf.db <- getArg("orfDB",required=TRUE)

# Turns on debug mode
debug.mode <- as.logical(getArg("debug",default=FALSE))

#Create logger
log.file <- paste(dir.name,"consolidator_",job.id,".log",sep="")
logger <- new.logger(log.file)

#####
# STEP 1: Consolidate results
#####
system(paste("cat ",dir.name,"R1_*.fastq>",dir.name,"R1.fastq&&rm R1_*.fastq",sep=""))
system(paste("cat ",dir.name,"R2_*.fastq>",dir.name,"R2.fastq&&rm R2_*.fastq",sep=""))
system(paste("cat ",dir.name,"BC_*.fastq>",dir.name,"BC.fastq&&rm BC_*.fastq",sep=""))

r1.file <- paste(dir.name,"R1.fastq",sep="")
r2.file <- paste(dir.name,"R2.fastq",sep="")
bc.file <- paste(dir.name,"BC.fastq",sep="")

#####
# STEP 2: Align to ORFs
#####
orf.sam <- bowtie(r1.file, orf.db, purge=FALSE, debug.mode=debug.mode)

#####
# STEP 3: Assemble the barcode
#####
