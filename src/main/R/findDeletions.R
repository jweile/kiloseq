source("lib/sge2.R")  
source("lib/liblogging.R") 
source("lib/cliargs.R")      #Command-line argument processing


logger <- new.logger("findDeletions.log")

orf.db <- getArg("orfDB",required=TRUE)

ctfile <- getArg("cloneTable",required=TRUE)
base.dir <- getArg("baseDir",default="ks_joint")

mql <- as.integer(getArg("maxQueue",default=80))

out.dir <- getArg("outDir",default="findDeletions")

clone.table <- read.csv(ctfile,stringsAsFactors=FALSE)
dir.names <- with(clone.table,paste(base.dir,set,plate,well,sep="/"))

make.chunks <- function(x,size=100) split(x, ceiling(seq_along(x)/size))

dir.create(out.dir)

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


