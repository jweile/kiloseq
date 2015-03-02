#!$Rbin

####
# demuxer.R is a worker script that runs on the cluster nodes
# It first aligns the tag reads to the tag library, then
# sorts the associated read pairs into appropriate bins
# It also tries to find the position of the barcode with in the 
# tag read and extracts the barcode sequence, adding to the bin.
# 
# Written by Jochen Weile <jochenweile@gmail.com> and Anjali Gopal <anjali.gopal91@gmail.com>

library("Biostrings")

source("lib/liblogging.R")   #Logger
source("lib/cliargs.R")      #Command-line argument processing
source("lib/libyogitools.R") #Helper functions

# R1 reads chunk file.
r1.file <- getArg("r1",required=TRUE)
# R2 reads chunk file.
r2.file <- getArg("r2",required=TRUE)
# Working directory
dir.name <- getArg("dir",required=TRUE)
# Job ID
job.id <- getArg("id",required=TRUE)
# Location of well tag DB
welltag.db <- getArg("welltags",default="res/welltags")
# Location of DNTAG DB
dntag.db <- getArg("dntags",default="res/dntags")

# Turns on debug mode
debug.mode <- as.logical(getArg("debug",default=FALSE))


#Create logger
log.file <- paste(dir.name,"demuxer_",job.id,".log",sep="")
logger <- new.logger(log.file)

#Load sequences
r1.seq <- read.DNAStringSet(r1.file)
r2.seq <- read.DNAStringSet(r2.file)

#Run Bowtie on R2 file against welltag DB and extract well information
logger$info("Aligning to well tags...")
welltag.sam <- bowtie(r2.file,welltag.db,debug.mode=debug.mode)
#Extract Well info
wells <- apply(
	extract.groups(welltag.sam$rname,"SET_ID=(\\w{1})\\|WELL=(\\w{1}\\d{2})"), 
	1, 
	function(groups) {
		if (!any(is.na(groups))) paste(groups,collapse="_") else NA
	}
)

#####
#TODO: Check that sam entries are in same order as fastq!
#####

#Run Bowtie on R2 file against DNTAG site to find Barcodes
logger$info("Aligning to DN tags...")
dntag.sam <- bowtie(r2.file,dntag.db,debug.mode=debug.mode)

#CIGAR: S=Soft clip, H=Hard clip, N=Intron skip, M=Match, D=Deletion, I=Insertion, P=Padded
cigar <- global.extract.groups(dntag.sam$cigar,"(\\d+)([SHNMDIP]{1})")
clip.width <- lapply(
	cigar,
	function(cigar) if (cigar[1,2]=="S") as.numeric(cigar[1,1]) else 0
)
match.width <- lapply(
	cigar, function(cigar) sum(as.numeric(cigar[1,!(cigar[,2] %in% c("D","S")) )))
)
#Extract barcode sequences
barcode.seq <- subseq(r2.seq,clip.width+1,clip.width+match.width)

####
# SORT READS INTO WELLS
#
logger$info("Sorting reads into wells...")
tapply(1:length(wells), wells, function(idx) {
	well <- wells[[ idx[[1]] ]]
	sub.dir <- paste(dir.name,well,sep="")
	if (!file.exists(sub.dir)) dir.create(sub.dir,showWarnings=FALSE)
	r1.file <- paste(sub.dir,"R1_",job.id,".fastq",sep="")
	r2.file <- paste(sub.dir,"R2_",job.id,".fastq",sep="")
	bc.file <- paste(sub.dir,"BC_",job.id,".txt")

	write.XStringSet(r1.seq[idx],r1.file,format="fastq")
	write.XStringSet(r2.seq[idx],r2.file,format="fastq")
})


###
# Delete input files when we're done, to save HD space.
#
if (!debug.mode) {
	logger$info("Cleaning up sequence fragments")
	file.remove(r1.file)
	file.remove(r2.file)
}
