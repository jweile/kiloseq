
source("lib/liblogging.R")   #Logger
source("lib/cliargs.R")      #Command-line argument processing
source("lib/libyogitools.R") #Helper functions
source("lib/libyogiseq.R")

orf.db <- getArg("orfDB",required=TRUE)
orf.fa <- paste(orf.db,".fa",sep="")
orf.bed <- paste(orf.db,".bed",sep="")

con <- file(orf.fa,open="r")
ref.seq <- readFASTA(con)[[1]]
close(con)
ref.length <- length(ref.seq)

id <- getArg("id",required=TRUE)
logger <- new.logger(paste(out.dir,"/",id,".log",sep="")

dirsArg <- getArg("dirs",required=TRUE)

out.dir <- getArg("outDir",required=TRUE)

dirs <- strsplit(dirsArg,",")[[1]]
sam.files <- paste(dirs,"OR.sam",sep="/")

chunk.info <- do.call(rbind,strsplit(dir.chunk,"/"))
chunk.info <- chunk.info[,-1]
colnames(chunk.info) <- c("set","plate","well")

ds <- lapply(sam.files, function(sam.file) {
	pu <- sam2pileup(sam.file,orf.fa)
	bases <- c("A","C","G","T","*")
	freqs <- do.call(rbind,lapply(pu$pileup, function(pile.i) {
	    fpile <- pile.i[pile.i$p < .05,]
	    table(factor(fpile$base,levels=bases))
	}))
	d <- apply(freqs,1,sum) + sapply(pu$indel.track,length)
	names(d) <- names(pu$pileup)
	sapply(as.character(1:ref.length), function(pos) if (pos %in% names(d)) d[[pos]] else 0)
})


border.detect <- function(d) sapply(1:length(d), function(i){
	left <- if (i < 10) i-1 else 10
	right <- if (i > length(d)-10) length(d)-i else 10
	mean(d[(i-left):(i-1)]) - mean(d[(i+1):(i+right)])
})

dmat <- do.call(rbind,ds)
dnorm <- t(apply(dmat,1,function(d)d/sum(d)))
dmean <- apply(dnorm,2,mean)
ddev <- t(apply(dnorm,1,`-`,dmean))
edges <- t(apply(ddev,1,border.detect))

hits <- apply(edges,1,function(x)any(abs(na.omit(x)) > 0.0003))

out <- cbind(chunk.info,large.indel=hits)

write.table(out,paste(out.dir,"/",id,".csv",sep=""),sep=",",row.names=FALSE)


