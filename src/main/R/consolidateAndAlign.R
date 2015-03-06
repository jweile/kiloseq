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
logger$info("Consolidating results...")
system(paste("cat ",dir.name,"R1_*.fastq>",dir.name,"R1.fastq&&rm ",dir.name,"R1_*.fastq",sep=""))
system(paste("cat ",dir.name,"R2_*.fastq>",dir.name,"R2.fastq&&rm ",dir.name,"R2_*.fastq",sep=""))
system(paste("cat ",dir.name,"BC_*.fastq>",dir.name,"BC.fastq&&rm ",dir.name,"BC_*.fastq",sep=""))

r1.file <- paste(dir.name,"R1.fastq",sep="")
r2.file <- paste(dir.name,"R2.fastq",sep="")
bc.file <- paste(dir.name,"BC.fastq",sep="")

#####
# STEP 2: Assemble the barcode
#####

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

logger$info("Reading barcodes...")
bcs <- read.fastq(bc.file)
bc.seqs <- sapply(bcs,function(x)x$toString())
#Compute frequences of each sequence
bc.freqs <- table(bc.seqs)
#Pick sequences with > 1% occurrence
top.freqs <- bc.freqs[bc.freqs > 0.01*length(bc.seqs)]

logger$info("Clustering barcodes...")
#Cluster based on edit distance
cm <- new.cluster.map(length(top.freqs))
for (i in 2:length(top.freqs)) {
	seq.i <- names(top.freqs)[i]
	for (j in 1:(i-1)) {
		seq.j <- names(top.freqs)[j]
		al <- new.alignment(seq.i,seq.j)
		if (al$getDistance() < 3) {
			cm$addLink(i,j)
		}
	}
}
top.clusters <- do.call(rbind,lapply(cm$getClusters(), function(cluster) {
	top <- names(top.freqs[cluster])[[which.max(top.freqs[cluster])]]
	list(seq=top,freq=sum(top.freqs[cluster])/length(bc.seqs))
}))
# top.clusters <- top.clusters[order(unlist(top.clusters[,2]),decreasing=TRUE),],
top.clusters <- data.frame(seq=unlist(top.clusters[,1]),freq=unlist(top.clusters[,2]),stringsAsFactors=FALSE)
top.clusters <- top.clusters[order(top.clusters$freq,decreasing=TRUE),]
write.table(
	top.clusters,
	paste(dir.name,"barcodes.csv",sep=""),
	sep=",",quote=FALSE,row.names=FALSE
)


#####
# STEP 3: Align to ORFs
#####
orf.sam <- bowtie(r1.file, orf.db, purge=FALSE, debug.mode=debug.mode)
