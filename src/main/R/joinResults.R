
source("lib/libyogitools.R")
source("lib/sge2.R")  
source("lib/liblogging.R") 
source("lib/cliargs.R")      #Command-line argument processing


logger <- new.logger("joinResults.log")

orf.db <- getArg("orfDB",default="idxs/ube2i-ad")

a.root <- getArg("rootA",default="kiloseq_2015-04-24/ube2i-emlr-ad-10x_2015-04-28_15-11-12")
a.pattern <- getArg("patternA",default="Plate(\\d+)")
b.root <- getArg("rootB",default="kiloseq_2015-03-23/ube2i-emlr-ad_2015-03-24_17-12-23")
b.pattern <- getArg("patternB",default="Plate-(\\d+)")
target.root <- getArg("target",default="ks_joint/ube2i-emlr-ad")
mql <- as.integer(getArg("maxQueue",default=80))

subdirs <- function(d) {
	fs <- list.files(d,full.name=TRUE)
	fs[file.info(fs)$isdir]
}

a.plate.dirs <- subdirs(a.root)
a.plates <- as.integer(extract.groups(a.plate.dirs,a.pattern))
b.plate.dirs <- subdirs(b.root)
b.plates <- as.integer(extract.groups(b.plate.dirs,b.pattern))

plates <- to.df(do.call(rbind,lapply(sort(a.plates),function(i) {
	list(
		plate.num=i,
		a.dir=a.plate.dirs[[which(a.plates==i)]],
		b.dir=b.plate.dirs[[which(b.plates==i)]]
	)
})))

all.wells <- apply(
	expand.grid(LETTERS[1:4],LETTERS[1:8],1:12),1,
	function(x)paste(x[[1]],"_",x[[2]],sprintf("%02d",as.numeric(x[[3]])),sep="")
)

sge <- new.sge(max.queue.length=mql, logger=logger, debug=FALSE)

target.well.dirs <- character(0)

with(plates,lapply(1:nrow(plates),function(i) {
	a.well.dirs <- subdirs(a.dir[[i]])
	a.well.dirs <- a.well.dirs[regexpr("[A-D]_[A-H]\\d{2}",a.well.dirs)>0]
	a.well.ids <- substr(a.well.dirs,nchar(a.well.dirs)-4,nchar(a.well.dirs))
	b.well.dirs <- subdirs(b.dir[[i]])
	b.well.dirs <- b.well.dirs[regexpr("[A-D]_[A-H]\\d{2}",b.well.dirs)>0]
	b.well.ids <- substr(b.well.dirs,nchar(b.well.dirs)-4,nchar(b.well.dirs))
	target.plate <- paste(target.root,"/plate",sprintf("%02d",plate.num[[i]]),sep="")

	for (well in all.wells) {
		if (!(well %in% a.well.ids) || !(well %in% b.well.ids)) {
			next
		}
		a.well.dir <- a.well.dirs[which(a.well.ids==well)]
		b.well.dir <- b.well.dirs[which(b.well.ids==well)]
		target.well.dir <- paste(target.plate,"/",well,sep="")
		dir.create(target.well.dir,recursive=TRUE)
		target.well.dirs[length(target.well.dirs)+1] <<- target.well.dir
		system(paste("cat ",a.well.dir,"/OR.fastq ",b.well.dir,"/OR.fastq >",target.well.dir,"/OR.fastq",sep=""))
		a.bc <- paste(a.well.dir,"/barcodes.csv",sep="")
		b.bc <- paste(b.well.dir,"/barcodes.csv",sep="")
		if (file.exists(a.bc)) {
			file.copy(a.bc,target.well.dir)
		} else if (file.exists(b.bc)) {
			file.copy(b.bc,target.well.dir)
		}
		id <- paste("joinAndAlign",i,well,sep="_")
		sge$enqueue(
			id=id,
			command="/software/R/bin/Rscript",
			arguments=list(
				"lib/alignWorker.R",
				paste("dir=",target.well.dir,sep=""),
				paste("id=",id,sep=""),
				paste("orfDB=",orf.db,sep="")
			)
		)

	}
}))
sge$wait(verbose=TRUE)

all.calls <- do.call(rbind,lapply(target.well.dirs,function(target.well.dir) {
	well.info <- extract.groups(target.well.dir,"/([^/]+)/([^/]+)$")
	colnames(well.info) <- c("plate","well")

	calls.file <- paste(target.well.dir,"/calls.csv",sep="")
	if (file.exists(calls.file)) {
		calls <- read.csv(calls.file,stringsAsFactors=FALSE)
		return(cbind(well.info,calls))
	} else {
		return(NULL)
	}
}))

write.table(all.calls,paste(target.root,"/calls.csv",sep=""),sep=",",row.names=FALSE)

