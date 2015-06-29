source("lib/sge2.R")  
source("lib/liblogging.R") 
source("lib/cliargs.R")      #Command-line argument processing

library("hash")
library("parallel")

logger <- new.logger("findDeletions.log")

orf.db <- getArg("orfDB",required=TRUE)

ctfile <- getArg("cloneTable",required=TRUE)
base.dir <- getArg("baseDir",default="ks_joint")

mql <- as.integer(getArg("maxQueue",default=80))
ncpu <- as.integer(getArg("maxCPU",default=16))

out.dir <- getArg("outDir",default="findDeletions")

clone.table <- read.csv(ctfile,stringsAsFactors=FALSE)
dir.names <- with(clone.table,paste(base.dir,set,plate,well,sep="/"))

make.chunks <- function(x,size=100) split(x, ceiling(seq_along(x)/size))

dir.create(out.dir)

logger$info("Starting jobs")

sge <- new.sge(max.queue.length=mql, logger=logger, debug=FALSE)
dir.chunks <- make.chunks(dir.names)
for (i in 1:length(dir.chunks)) {

	dir.chunk <- dir.chunks[[i]]

	target.well.dirs <- paste(dir.chunk,collapse=",")

	id <- paste("findDelWorker",i,sep="")

	sge$enqueue(
		id=id,
		command="/software/R/bin/Rscript",
		arguments=list(
			"lib/findDelWorker.R",
			paste("dirs=",target.well.dirs,sep=""),
			paste("orfDB=",orf.db,sep=""),
			paste("outDir=",out.dir,sep=""),
			paste("id=",id,sep="")
		)
	)

}
sge$wait(verbose=TRUE)


border.detect <- function(d) sapply(1:length(d), function(i){
	left <- if (i < 10) i-1 else 10
	right <- if (i > length(d)-10) length(d)-i else 10
	mean(d[(i-left):(i-1)]) - mean(d[(i+1):(i+right)])
})


logger$info("Collecting results")

outfiles <- list.files(out.dir,pattern=".csv",full.names=TRUE)
results <- do.call(rbind,lapply(outfiles,read.csv,stringsAsFactors=FALSE))

result.idx <- hash(paste(results[,1],results[,2],results[,3],sep="_"),1:nrow(results))
result.order <- values(result.idx,paste(clone.table[,1],clone.table[,2],clone.table[,3],sep="_"))
results <- results[result.order,]

logger$info("Detecting deletions")

dmat <- as.matrix(results[,-(1:3)])
#normalize by number of reads per well
dnorm <- t(apply(dmat,1,function(d)d/sum(d)))
#compute positional average
dmean <- apply(dnorm,2,mean)
# dmed <- apply(dnorm,2,median)
#compute deviation from average
ddev <- t(apply(dnorm,1,`-`,dmean))
#detect sharp edges in deviations
# edges <- t(apply(ddev,1,border.detect))
edges <- do.call(rbind,mclapply(
	lapply(1:nrow(ddev),function(i)ddev[i,]),
	border.detect, 
	mc.cores=ncpu
))

hits <- apply(edges,1,function(x)any(abs(na.omit(x)) > 0.0003))

delRanges <- sapply(1:nrow(dmat),function(i) if(hits[[i]]) paste(which.max(edges[i,]),which.min(edges[i,]),sep="-") else NA)
clone.table <- cbind(clone.table,deletions=delRanges)
write.table(clone.table,paste(out.dir,"cloneTable.csv",sep="/"),sep=",",row.names=FALSE)

foo <- t(apply(edges[hits,],1,function(x)c(start=which.max(x),end=which.min(x))))
bar <- foo[foo[,1]<foo[,2],]

bar <- bar[order(bar[,1],bar[,2]),]

pdf(paste(out.dir,"deletionMap.pdf",sep="/"),8.5,11)
op <- par(mar=c(5,4,4,10)+.1,las=1,cex=.7)
plot(0,type="n",xlim=c(0,621),ylim=c(0,nrow(bar)),axes=FALSE,xlab="bp",ylab="deletions")
axis(1,at=seq(0,621,100),labels=seq(0,621,100))
# axis(4,at=bar.plate.starts,labels=names(bar.plate.starts))
arrows(bar[,1],1:nrow(bar),bar[,2],1:nrow(bar),length=0.05)
grid(NULL,NA)
abline(v=c(77,550),col="red")
par(op)
dev.off()
