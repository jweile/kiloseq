source("lib/liblogging.R")   #Logger
source("lib/cliargs.R")      #Command-line argument processing
source("lib/libyogitools.R") #Helper functions
source("lib/libyogiseq.R")



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
orf.bed <- paste(orf.db,".bed",sep="")

#Create logger
log.file <- paste(dir.name,job.id,".log",sep="")
logger <- new.logger(log.file)

orf.read.file <- paste(dir.name,"OR.fastq",sep="")

ref.con <- file(orf.fa,open="r")
ref.seq <- readFASTA(ref.con)[[1]]
close(ref.con)

orf.pos <- read.delim(orf.bed,stringsAsFactors=FALSE,header=FALSE)
colnames(orf.pos) <- c("chrom","from","to","name")
orf.pos <- orf.pos[orf.pos$name=="ORF",c("from","to")]


bc.info <- read.csv(paste(dir.name,"barcodes.csv",sep=""),stringsAsFactors=FALSE)
bc.info <- bc.info[bc.info$freq>.25,]


doAlign <- function() {
	#Alignment
	logger$info("Aligning...")
	sam.file <- bowtie(orf.read.file, orf.db, 
		purge=FALSE, parse=FALSE, header=TRUE, short=FALSE,
		debug.mode=debug.mode
	)
	#Look into SAM file to compute alignment efficiency
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
			stringsAsFactors=FALSE
		)))
	}

	logger$info("Calling variants...")
	# Variant Caller
	variants <- call.variants(sam.file,orf.fa)
	dp5.abs <- sum(variants$depth < 5)
	dp5 <- 1- dp5.abs/length(ref.seq)
	if (dp5.abs > 10) {
		return(cbind(bc.info,data.frame(
			al.rate=al.rate,
			dp5=dp5,
			call="Low coverage!",
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

	if (nrow(vcalls)==0) {
		return(cbind(bc.info,data.frame(
			al.rate=al.rate,
			dp5=dp5,
			call="WT",
			stringsAsFactors=FALSE
		)))
	}

	depth <- variants$depth
	vcalls <- vcalls[vcalls$freq > .25,]
	vcalls <- vcalls[order(vcalls$pos),]
	mutstr <- paste(sapply(1:nrow(vcalls),function(i) paste(vcalls[i,1],vcalls[i,2],vcalls[i,3],sep="")),collapse=",")
	# mutstr <- toupper(paste(apply(vcalls[,1:3],1,paste,collapse=""),collapse=","))
	if (any(vcalls$pos < orf.pos$from | vcalls$pos > orf.pos$to)) {
		mutstr <- paste("Corrupt:",mutstr)
	}
	if (regexpr("[\\*\\+-]",mutstr)>0) {
		mutstr <- paste("Frameshift:",mutstr)
	}
	return(cbind(bc.info,data.frame(
		al.rate=al.rate,
		dp5=dp5,
		call=mutstr,
		stringsAsFactors=FALSE
	)))
}

out <- doAlign()


write.table(
	out,
	paste(dir.name,"calls.csv",sep=""),
	sep=",",quote=TRUE,row.names=FALSE
)

logger$info("Done!")
