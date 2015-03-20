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

# Whether or not to expect barcodes in the tag reads.
use.barcodes <- as.logical(getArg("useBarcodes",default=TRUE))

# Turns on debug mode
debug.mode <- as.logical(getArg("debug",default=FALSE))

#Create logger
log.file <- paste(dir.name,job.id,".log",sep="")
logger <- new.logger(log.file)


#####
# STEP 1: Consolidate results
#####

logger$info("Consolidating results...")
system(paste("cat ",dir.name,"OR_*.fastq>",dir.name,"OR.fastq&&rm ",dir.name,"OR_*.fastq",sep=""))
system(paste("cat ",dir.name,"TR_*.fastq>",dir.name,"TR.fastq&&rm ",dir.name,"TR_*.fastq",sep=""))
orf.read.file <- paste(dir.name,"OR.fastq",sep="")
tag.read.file <- paste(dir.name,"TR.fastq",sep="")

if (use.barcodes) {
	system(paste("cat ",dir.name,"BC_*.fastq>",dir.name,"BC.fastq&&rm ",dir.name,"BC_*.fastq",sep=""))
	bc.file <- paste(dir.name,"BC.fastq",sep="")
}

#####
# STEP 2: Assemble the barcodes
#####

if (use.barcodes) {

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
} else {
	#in case of no barcodes, just use a dummy table
	top.clusters <- data.frame(seq=NA,freq=NA,stringsAsFactors=FALSE)
}

####
# STEP 3: Segregate reads based on barcodes
####
if (use.barcodes) {

	logger$info("Loading sequence data...")
	orf.reads <- read.fastq(orf.read.file)
	names(orf.reads) <- sapply(orf.reads,function(x)x$getID())

	logger$info("Segregating reads by barcode...")
	segregated.or.files <- sapply(1:nrow(top.clusters), function(i) {
		read.ids <- cluster.readnames[[i]]
		reads <- orf.reads[read.ids]
		sub.file <- paste(dir.name,"OR_BC",i,".fastq",sep="")
		write.fastq(sub.file,reads)
		sub.file
	})
} else {
	#in case of no barcodes use a dummy file list
	segregated.or.files <- orf.read.file
}

#####
# STEP 4: Align to ORFs and call variants
#####
logger$info("Aligning reads to ORF...")

ref.con <- file(orf.fa,open="r")
ref.seq <- readFASTA(ref.con)[[1]]
close(ref.con)

all.calls <- do.call(rbind,lapply(1:nrow(top.clusters), function(cluster.idx) {

	orf.read.file <- segregated.or.files[[cluster.idx]]
	bc.info <- top.clusters[cluster.idx,]

	#Alignment
	sam.file <- bowtie(orf.read.file, orf.db, 
		purge=FALSE, parse=FALSE, header=TRUE, short=FALSE,
		debug.mode=debug.mode
	)
	#Look into SAM file to compute alignment effeciency
	sam <- read.delim(sam.file,stringsAsFactors=FALSE,comment.char="@",header=FALSE)
	colnames(sam) <- c(
		"cname","flag","rname","pos","mapq","cigar","mrnm","mpos",
		"isize","seq","qual","tags"
	)
	al.rate <- 1-(sum(sam$rname == "*")/nrow(sam))
	if (al.rate == 0) {
		return(cbind(bc.info,data.frame(
			al.rate=al.rate,
			dp5=0,
			call="No alignment!",
			share=NA,
			stringsAsFactors=FALSE
		)))
	}

	# Variant Caller
	variants <- call.variants(sam.file,orf.fa)
	dp5.abs <- sum(variants$depth < 5)
	dp5 <- 1- dp5.abs/length(ref.seq)
	if (dp5.abs > 10) {
		return(cbind(bc.info,data.frame(
			al.rate=al.rate,
			dp5=dp5,
			call="Low coverage!"
			,share=NA,
			stringsAsFactors=FALSE
		)))
	}

	#export variant calls
	vcalls <- variants$calls
	write.table(
		vcalls,
		paste(dir.name,"varcalls.csv",sep=""),
		sep=",",quote=FALSE,row.names=FALSE
	)

	#cluster calls based on frequency
	cm <- new.cluster.map(nrow(vcalls))
	for (i in 2:nrow(vcalls)) {
		for (j in 1:(i-1)) {
			fdist <- abs(vcalls$freq[[i]] - vcalls$freq[[j]])
			if (fdist < .05) {
				cm$addLink(i,j)
			}
		}
	}
	out <- do.call(rbind,lapply(cm$getClusters(), function(is) {
		freq <- mean(vcalls$freq[is])
		is <- is[order(vcalls$pos[is])]
		mutstr <- paste(sapply(is, function(i) {
			with(vcalls, paste(ref[[i]],pos[[i]],alt[[i]],sep=""))
		}),collapse=",")
		cbind(bc.info,data.frame(
			al.rate=al.rate,
			dp5=dp5,
			call=mutstr,
			share=freq,
			stringsAsFactors=FALSE
		))
	}))

	return(out[order(out$share,decreasing=TRUE),])

	# vcf <- call.variants(sam.file,orf.fa)
	# #Process results
	# dp <- as.numeric(extract.groups(vcf$info,"DP=(\\d+)")[,1])
	# depth <- sapply(1:length(ref.seq), function(i) {
	# 	if (i %in% vcf$pos) {
	# 		dp[[which(vcf$pos == i)[[1]]]]
	# 	} else {
	# 		0
	# 	}
	# })
	# dp5 <- sum(depth > 5)/length(ref.seq)
	# if (sum(depth < 5) > 10) {
	# 	return(list(call="Low coverage!",al.rate=al.rate,dp5=dp5))
	# }
	# #TODO: Exctract significant SNPs and translate!
	# return(list(call="Good!",al.rate=al.rate,dp5=dp5))
}))

# top.clusters[,"al.rate"] <- unlist(out[,"al.rate"])
# top.clusters[,"dp5"] <- unlist(out[,"dp5"])
# top.clusters[,"call"] <- unlist(out[,"call"])

write.table(
	all.calls,
	paste(dir.name,"calls.csv",sep=""),
	sep=",",quote=FALSE,row.names=FALSE
)


logger$info("Done!")
