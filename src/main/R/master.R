#!$Rbin

####
# master.R is the top-level controller script that chops large fastq
# files into smaller chunks on-the-fly and spawns worker jobs on the 
# SGE cluster to process them. Afterwards it collates the results for
# further processing.
# 
# Written by Jochen Weile <jochenweile@gmail.com> and Anjali Gopal <anjali.gopal91@gmail.com>

library("hash",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE)

source("lib/sge2.R")         #SUN Grid Engine
source("lib/cliargs.R")      #Command-line argument processing
source("lib/liblogging.R")   #Writing log files
source("lib/libyogitools.R") #Some handy tools

###
# Get command line arguments
#

# Files containing the R1 and R2 reads. Example file:
r1.files <- getArg("r1",required=TRUE)
r2.files <- getArg("r2",required=TRUE)
if (length(strsplit(r1.files,",")[[1]]) != length(strsplit(r2.files,","))) {
	stop("There needs to be an equal number of r1 and r2 files!")
}
rfile.table <- cbind(r1=strsplit(r1.files,",")[[1]],r2=strsplit(r2.files,",")[[1]])


# The session tag is used to name the ouptut directory together with a timestamp.
session.tag <- getArg("session",default="kiloseq")

# The chuck size determines how many reads are processed by each slave script.
chunk.size <- as.numeric(getArg("chunksize",default=20000))

# Turns on debug mode
debug.mode <- as.logical(getArg("debug",default=FALSE))

# This option determines how many jobs will be submitted to the cluster at maximum.
# When more jobs are available, they get buffered internally until enough room exists in the queue.
max.queue <- as.numeric(getArg("maxQueue",default=30))

# Location of well tag DB
welltag.db <- getArg("welltags",default="res/welltags")
# # Location of DNTAG DB
# dntag.db <- getArg("dntags",default="res/dntags")
# Location of ORF DB
orf.db <- getArg("orfDB",default="res/orfs")


###
# Create output directory
#
timestamp <- format(Sys.time(),format='%Y-%m-%d_%H-%M-%S')
out.dir <- paste(session.tag,"_",timestamp,"/", sep="")
dir.create(out.dir, mode="0755")

# Set up log file
logger <- new.logger(paste(out.dir,"master.log",sep=""))



###
# PHASE 1: CREATE CHUNKS OF READS ON-THE-FLY AND DEMUX THEM
#

###
# This function uses an already open file connection
# to create a new FASTQ chunk file. R
#
# incon = incoming file connection
# out.dir = output directory
# i = chunk number
# direction = R1 or R2
#
make.chunk <- function(incon, out.dir, i, direction) {

	outfile <- paste(out.dir,direction,"-",i,".fastq",sep="")
	outcon <- file(outfile,open="w")

	lines <- readLines(incon, chunk.size*4)
	writeLines(lines,outcon)

	close(outcon)

	list(file=outfile, last=(length(lines) < chunk.size*4))
}

###
# Function to find out whether file is GZip file
#
is.gz <- function(f) {
	substr(f,nchar(f)-2,nchar(f))==".gz" || 
		regexpr("gzip compressed data", system(paste("file",f),intern=TRUE) ) > 0
}

###
# work through read file pairs
#
# Construct SunGridEngine object with designated maximum queue size.
sge <- new.sge(max.queue.length=max.queue, logger=logger, debug=debug.mode)

logger$info("Demultiplexing...")

#Iterate over all SWIM wells, create chunks and submit jobs on the go.
#The return result directories.
result.dirs <- apply(rfile.table, 1, function(rfiles) {

	r1.file <- rfiles[["r1"]]
	r2.file <- rfiles[["r2"]]

	#TODO: name directory according to SWIM well
	dir.name <- paste(out.dir,gsub(".+/|\\.fastq(\\.gz)?$","",r1.file),"/",sep="")
	dir.create(dir.name)

	###
	# Open connections to sequencing results files
	#
	con.r1 <- file(r1.file,open="r")
	if (is.gz(r1.file)) {
		con.r1 <- gzcon(con.r1)
	}
	con.r2 <- file(r2.file,open="r")
	if (is.gz(r2.file)) {
		con.r2 <- gzcon(con.r2)
	}


	done <- FALSE
	i <- 0
	while (!done) {
		i <- i+1
		# Make chunks for R1 and R2
		r1.chunk <- make.chunk(con.r1,dir.name,i,"R1")
		r2.chunk <- make.chunk(con.r2,dir.name,i,"R2")

		#Create a job id
		job.id <- paste(session.tag,timestamp,i,sep="_")
		#Designate a log file.
		slave.log <- paste(dir.name,"slave_",i,".log",sep="")

		#Submit Slave job to SunGridEngine
		sge$enqueue(
			id=job.id,
			command="$Rbin",
			arguments=list(
				"lib/demuxer.R",
				paste("r1=",r1.chunk$file,sep=""),
				paste("r2=",r2.chunk$file,sep=""),
				paste("dir=",dir.name,sep=""),
				paste("id=",job.id,sep=""),
				paste("welltags=",welltag.db,sep=""),
				# paste("dntags=",dntag.db,sep=""),
				paste("debug=",debug.mode,sep="")
			)
		)

		#We're done if we run out of reads to process
		done <- r1.chunk$last || r2.chunk$last
	}
	#Close file connections
	close(con.r1)
	close(con.r2)

	dir.name
})
#Wait for the remaining jobs to finish
sge$wait(verbose=TRUE)


####
# PHASE 2: CONSOLIDATE JOB RESULTS AND ALIGN
#
logger$info("Consolidating and aligning...")

invisible(lapply(result.dirs, function(dir.name) {

	sub.dirs <- list.dirs(dir.name)

	lapply(sub.dirs, function(sub.dir) {

		swim.id <- gsub(".*/","",sub("/$","",dir.name))
		well.id <- gsub(".*/","",sub.dir)
		#ignore "undetermined" and "invalid" folders
		if (regexpr("[A-D]_\\w\\d{2}",well.id) < 0) {
			return(NULL)
		}
		job.id <- paste("consolidate",swim.id,well.id,sep="_")

		#Submit job to Sun Grid Engine
		sge$enqueue(
			id=job.id,
			command="$Rbin",
			arguments=list(
				"lib/consolidateAndAlign.R",
				paste("dir=",sub.dir,sep=""),
				paste("id=",job.id,sep=""),
				paste("orfDB=",orf.db,sep=""),
				paste("debug=",debug.mode,sep="")
			)
		)
	})
}))
#Wait for the remaining jobs to finish
sge$wait(verbose=TRUE)


logger$info("Done!")
