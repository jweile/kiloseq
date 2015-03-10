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
source("lib/libyogiseq.R")

#Function for safely loading sequences
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

# Working directory
dir.name <- getArg("dir",required=TRUE)
if (regexpr("/$",dir.name) < 1) {
	dir.name <- paste(dir.name,"/",sep="")
}
# Job ID
job.id <- getArg("id",required=TRUE)
# ORF Sequence DB
orf.db <- getArg("orfDB",required=TRUE)
orf.fa <- paste(orf.db,".fa",sep="")

# Turns on debug mode
debug.mode <- as.logical(getArg("debug",default=FALSE))

#Create logger
log.file <- paste(dir.name,job.id,".log",sep="")
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
# STEP 2: Assemble the barcodes
#####

logger$info("Reading barcodes...")
bcs <- read.fastq(bc.file)
bc.seqs <- sapply(bcs,function(x)x$toString())
bc.names <- sapply(bcs,function(x)x$getID())
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

#which reads belong to which cluster?
cluster.readnames <- lapply(cm$getClusters(), function(cluster) {
	seqs <- names(top.freqs)[cluster]
	bc.names[which(bc.seqs %in% seqs)]
})

#what is the barcode sequence of each cluster?
top.clusters <- do.call(rbind,lapply(cm$getClusters(), function(cluster) {
	top <- names(top.freqs[cluster])[[which.max(top.freqs[cluster])]]
	list(seq=top,freq=sum(top.freqs[cluster])/length(bc.seqs))
}))
top.clusters <- data.frame(seq=unlist(top.clusters[,1]),freq=unlist(top.clusters[,2]),stringsAsFactors=FALSE)

#sort by frequency
c.order <- order(top.clusters$freq,decreasing=TRUE)
top.clusters <- top.clusters[c.order,]
cluster.readnames <- cluster.readnames[c.order]

#write barcodes to file
write.table(
	top.clusters,
	paste(dir.name,"barcodes.csv",sep=""),
	sep=",",quote=FALSE,row.names=FALSE
)

####
# STEP 3: Segregate reads based on barcodes
####
logger$info("Loading sequence data...")
r1.reads <- read.fastq(r1.file)
names(r1.reads) <- sapply(r1.reads,function(x)x$getID())

logger$info("Segregating reads by barcode...")
segregated.r1.files <- sapply(1:nrow(top.clusters), function(i) {
	read.ids <- cluster.readnames[[i]]
	reads <- r1.reads[read.ids]
	sub.file <- paste(dir.name,"R1_BC",i,".fastq",sep="")
	write.fastq(sub.file,reads)
	sub.file
})

#####
# STEP 4: Align to ORFs and call variants
#####
logger$info("Aligning reads to ORF...")

ref.con <- file(orf.fa,open="r")
ref.seq <- readFASTA(ref.con)[[1]]
close(ref.con)

out <- do.call(rbind,lapply(segregated.r1.files, function(r1.file) {

	#Alignment
	sam.file <- bowtie(r1.file, orf.db, 
		purge=FALSE, parse=FALSE, header=TRUE, short=FALSE,
		debug.mode=debug.mode
	)
	sam <- read.delim(sam.file,stringsAsFactors=FALSE,comment.char="@",header=FALSE)
	colnames(sam) <- c(
		"cname","flag","rname","pos","mapq","cigar","mrnm","mpos",
		"isize","seq","qual","tags"
	)
	al.rate <- 1-(sum(sam$rname == "*")/nrow(sam))
	if (al.rate == 0) {
		return(list(call="No alignment!",al.rate=al.rate,dp5=0))
	}

	#Variant Caller
	vcf <- call.variants(sam.file,orf.fa)
	#Process results
	dp <- as.numeric(extract.groups(vcf$info,"DP=(\\d+)")[,1])
	depth <- sapply(1:length(ref.seq), function(i) {
		if (i %in% vcf$pos) {
			dp[[which(vcf$pos == i)[[1]]]]
		} else {
			0
		}
	})
	dp5 <- sum(depth > 5)/length(ref.seq)
	if (sum(depth < 5) > 10) {
		return(list(call="Low coverage!",al.rate=al.rate,dp5=dp5))
	}
	#TODO: Exctract significant SNPs and translate!
	return(list(call="Good!",al.rate=al.rate,dp5=dp5))
}))

top.clusters[,"al.rate"] <- unlist(out[,"al.rate"])
top.clusters[,"call"] <- unlist(out[,"call"])
top.clusters[,"dp5"] <- unlist(out[,"dp5"])

write.table(
	top.clusters,
	paste(dir.name,"calls.csv",sep=""),
	sep=",",quote=FALSE,row.names=FALSE
)


logger$info("Done!")
