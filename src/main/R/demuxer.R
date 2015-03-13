#!$Rbin

####
# demuxer.R is a worker script that runs on the cluster nodes
# It first aligns the tag reads to the tag library, then
# sorts the associated read pairs into appropriate bins
# It also tries to find the position of the barcode with in the 
# tag read and extracts the barcode sequence, adding to the bin.
# 
# Written by Jochen Weile <jochenweile@gmail.com> and Anjali Gopal <anjali.gopal91@gmail.com>

# library("Biostrings")

source("lib/liblogging.R")   #Logger
source("lib/cliargs.R")      #Command-line argument processing
source("lib/libyogitools.R") #Helper functions
source("lib/libyogiseq.R")   #FASTQ and bowtie

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
# Sequence snippet preceding the barcode in the tag read.
dntag.snippet <- getArg("snippet",default="TAGTGCGATTG")
# Whether or not to expect barcodes in the tag reads.
use.barcodes <- as.logical(getArg("useBarcodes",default=TRUE))
# Which read is considered the tag read? R1 or R2 ?
tag.orientation <- getArg("tagOrientation",default="R2")

# Turns on debug mode
debug.mode <- as.logical(getArg("debug",default=FALSE))

#Create logger
log.file <- paste(dir.name,"demuxer_",job.id,".log",sep="")
logger <- new.logger(log.file)

tag.read.file <- ifelse(tag.orientation=="R1",r1.file,r2.file)
orf.read.file <- ifelse(tag.orientation=="R1",r2.file,r1.file)

#Load sequences
read.fastq <- function(f) {
	tryCatch({
		con <- file(f, open="r")
		p <- new.fastq.parser(con)
		out <- list()
		while (length(s <- p$parse.next(1)) > 0) {
			out[[length(out)+1]] <- s[[1]]
		}
		out
	},
	error = function(ex) {
		logger$fatal(paste("Error reading file",f," :\n",ex))
		stop(ex)
	},
	finally = {
		if (exists("con") && isOpen(con)) {
			close(con)
		}
	})
}

seqnames <- function(seqs) sapply(seqs,function(s)s$getID())

logger$info("Loading sequence data...")
tag.read.seq <- read.fastq(tag.read.file)
orf.read.seq <- read.fastq(orf.read.file)

if (!all(seqnames(tag.read.seq)==seqnames(orf.read.seq))) {
	logger$fatal("R1 and R2 reads do not correspond to each other!")
	stop()
}

#####
# STEP 1: Run Bowtie on tag read file against welltag DB and extract well information
#####
logger$info("Aligning to well tags...")
#bowtie() function is defined in libyogitools.R
#TODO: adjust clip size automatically based on average read length!
welltag.sam <- bowtie(tag.read.file,welltag.db,clip3=47,debug.mode=debug.mode)
#Extract Well info
wells <- apply(
	extract.groups(welltag.sam$rname,"SET_ID=(\\w{1})\\|WELL=(\\w{1}\\d{2})"), 
	1, 
	function(groups) {
		# if (!any(is.na(groups))) paste(groups,collapse="_") else NA
		if (any(is.na(groups))) {
			"undetermined"
		} else if (!(groups[[1]] %in% LETTERS[1:4])) {
			"invalid"
		} else {
			paste(groups,collapse="_")
		}
	}
)

#sort wells according to order in tag read file
read.order <- sapply(seqnames(tag.read.seq), function(name) which(welltag.sam$cname == name))
wells <- wells[read.order]

#####
# STEP 2: Extract Barcodes
#####
if (use.barcodes) {
	seqs <- sapply(tag.read.seq,function(s)s$toString())
	barcode.seq <- mapply(function(m,s) {
		if (m > 0) subseq(s,m+nchar(dntag.snippet),length(s)) else NA
	},m=regexpr(dntag.snippet, seqs),s=tag.read.seq)
}


#Function for safely writing sequences
write.fastq <- function(f,seqs) {
	tryCatch({
		con <- file(f, open="w")
		writeFASTQ(con,seqs)
	},
	error = function(ex) {
		logger$fatal(paste("Error while writing file",f," :\n",ex))
		stop(ex)
	},
	finally = {
		if (exists("con") && isOpen(con)) {
			close(con)
		}
	})
}


#####
# STEP 3: SORT READS INTO WELLS
#####
logger$info("Sorting reads into wells...")
tapply(1:length(wells), wells, function(idx) {
	well <- wells[[ idx[[1]] ]]
	sub.dir <- paste(dir.name,well,"/",sep="")
	if (!file.exists(sub.dir)) dir.create(sub.dir,showWarnings=FALSE)

	orf.file <- paste(sub.dir,"OR_",job.id,".fastq",sep="")
	tag.file <- paste(sub.dir,"TR_",job.id,".fastq",sep="")
	write.fastq(orf.file,orf.read.seq[idx])
	write.fastq(tag.file,tag.read.seq[idx])

	if (use.barcodes) {
		bc.file <- paste(sub.dir,"BC_",job.id,".fastq",sep="")
		barcodes <- barcode.seq[idx]
		barcodes <- barcodes[!is.na(barcodes)]
		if (length(barcodes) > 0) {
			write.fastq(bc.file,barcodes)
		}
	}
})


#####
# CLEANUP: Delete input files when we're done, to save HD space.
####
if (!debug.mode) {
	logger$info("Cleaning up sequence fragments")
	file.remove(orf.read.file)
	file.remove(tag.read.file)
}