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
# Location of DNTAG DB
dntag.db <- getArg("dntags",default="res/dntags")

# Turns on debug mode
debug.mode <- as.logical(getArg("debug",default=FALSE))

#Create logger
log.file <- paste(dir.name,"demuxer_",job.id,".log",sep="")
logger <- new.logger(log.file)

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
r1.seq <- read.fastq(r1.file)
r2.seq <- read.fastq(r2.file)

if (!all(seqnames(r1.seq)==seqnames(r2.seq))) {
	logger$fatal("R1 and R2 reads do not correspond to each other!")
	stop()
}

#####
# STEP 1: Run Bowtie on R2 file against welltag DB and extract well information
#####
logger$info("Aligning to well tags...")
#bowtie() function is defined in libyogitools.R
welltag.sam <- bowtie(r2.file,welltag.db,clip3=47,debug.mode=debug.mode)
#Extract Well info
wells <- apply(
	extract.groups(welltag.sam$rname,"SET_ID=(\\w{1})\\|WELL=(\\w{1}\\d{2})"), 
	1, 
	function(groups) {
		if (!any(is.na(groups))) paste(groups,collapse="_") else NA
	}
)

read.order <- sapply(seqnames(r2.seq), function(name) which(welltag.sam$cname == name))
wells <- wells[read.order]

#####
# STEP 2: Run Bowtie on R2 file against DNTAG site to find Barcodes
#####
dntag.snippet <- "TAGTGCGATTG"
seqs <- sapply(r2.seq,function(s)s$toString())
barcode.seq <- mapply(function(m,s) {
	if (m > 0) subseq(s,m+11,length(s)) else NA
},m=regexpr(dntag.snippet, seqs),s=r2.seq)

# logger$info("Aligning to DN tags...")
# dntag.sam <- bowtie(r2.file,dntag.db,debug.mode=debug.mode)

# #CIGAR: S=Soft clip, H=Hard clip, N=Intron skip, M=Match, D=Deletion, I=Insertion, P=Padded
# cigar <- global.extract.groups(dntag.sam$cigar,"(\\d+)([SHNMDIP]{1})")
# #How many bases are clipped before the match?
# #The clipped bases should be the welltag and a bit of the loxP site.
# clip.width <- lapply(
# 	cigar,
# 	function(cigar) if (cigar[1,2]=="S") as.numeric(cigar[1,1]) else 0
# )
# #How long is the match in the read, including deletions?
# match.width <- lapply(
# 	cigar, function(cigar) sum(as.numeric(cigar[1,!(cigar[,2] %in% c("M","D")) )))
# )
# #the barcode starts at the end of the match
# bc.start <- clip.width+match.width+1
# #Extract barcode sequences. It's 25bp wide
# barcode.seq <- lapply(1:length(r2.seq), function(i) {
# 	bc.end <- bc.start[[i]]+25-1
# 	if (length(r2.seq[[i]] < bc.end)) bc.end <- length(r2.seq[[i]])
# 	subseq(r2.seq[[i]],bc.start[[i]],bc.end)
# })


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
	r1.file <- paste(sub.dir,"R1_",job.id,".fastq",sep="")
	r2.file <- paste(sub.dir,"R2_",job.id,".fastq",sep="")
	bc.file <- paste(sub.dir,"BC_",job.id,".fastq",sep="")

	barcodes <- barcode.seq[idx]
	barcodes <- barcodes[!is.na(barcodes)]

	write.fastq(r1.file,r1.seq[idx])
	write.fastq(r2.file,r2.seq[idx])
	write.fastq(bc.file,barcodes)
})


#####
# CLEANUP: Delete input files when we're done, to save HD space.
####
if (!debug.mode) {
	logger$info("Cleaning up sequence fragments")
	file.remove(r1.file)
	file.remove(r2.file)
}
