
source("../src/main/R/libyogitools.R") #Helper functions
source("../src/main/R/libyogiseq.R")	#FASTA parser

# con <- file("ube2i-ad.fa",open="r")
# ref.fa <- readFASTA(con)[[1]]
# close(con)

# ref.bed <- read.delim("ube2i-ad.bed",header=FALSE)
# colnames(ref.bed) <- c("ref","from","to","feature")

		

# muts <- strsplit(clones$aa.calls,",")


calc.reachables <- function(ref.bed,ref.fa) {
	tr <- init.translator()
	codon.starts <- seq(ref.bed$from,ref.bed$to,3)
	lapply(codon.starts, function(start) {
		codon <- subseq(ref.fa,start,start+2)$toString()
		from.aa <- tr$translate(codon)
		setdiff(unique(do.call(c,lapply(1:3, function(pos){
			from <- substr(codon,pos,pos)
			sapply(setdiff(c("A","C","G","T"),from),function(to) {
				to.codon <- codon
				substr(to.codon,pos,pos) <- to
				tr$translate(to.codon)
			})
		}))),from.aa)
	})
}

# Plots the mutation coverage for a given change matrix
plotMutCoverage <- function(.calls, ref.seq, ref.bed, main="") {

	# num.aa <- nchar(sequence)/3
	num.aa <- (ref.bed$to - ref.bed$from + 1)/3

	reachables <- calc.reachables(ref.bed,ref.seq)

	muts <- do.call(c,lapply(strsplit(.calls$aa.calls,","), function(ms) {
		if (regexpr("No|Low|WT|Corr|X-Over|Dupl|Frame|Stop",ms[[1]])>0) {
			return(NULL)
		}
		ms
	}))

	#initialize the change matrix
	change.matrix <- matrix(0,nrow=21,ncol=num.aa,
		dimnames=list(
			c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'),
			1:num.aa
		)
	)
	
	# Mark mutations in the matrix
	for (sac in muts) {

		if (sac == "silent") next

		pos <- as.numeric(substr(sac,2,nchar(sac)-1))
		aa <- substr(sac,nchar(sac),nchar(sac))
		from <- substr(sac,1,1)

		change.matrix[aa,pos] <- change.matrix[aa,pos] + 1
		change.matrix[from,pos] <- -1
	}

	coverage.all <- apply(change.matrix,2,function(x) sum(na.omit(x) > 0)) / 19
	coverage.lev <- sapply(1:num.aa, function(i) {
		sum(change.matrix[reachables[[i]],i] > 0)/length(reachables[[i]])
	})

	#define drawing layot, set drawing color to gray, adjust margins
	layout(cbind(1:3,c(5,5,4)), heights=c(1,.6,3),widths=c(9,.5))
	op <- par(fg="gray",mar=c(0,4.1,4.1,2.1),xaxs="i")	

	# draw a bar plot for coverage.all
	barplot(coverage.all,
		main=main, 
		xlab="Position",
		ylab="All", 
		ylim=c(0,1),
		border=NA,
		names.arg=NA,
		col="darkolivegreen3",
		axes=FALSE
	)
	axis(2,at=c(0,.5,1),labels=c(0,.5,1))
	grid(NA,4)

	op <- par(mar=c(0,4.1,1,2.1))	
	barplot(coverage.lev,
		xlab="Position",
		ylab="Access.", 
		ylim=c(0,1),
		border=NA,
		names.arg=NA,
		col="darkolivegreen4",
		axes=FALSE
	)
	axis(2,at=c(0,.5,1),labels=c(0,.5,1))
	grid(NA,4)


	# Compute a color gradient to represent the mutation counts
	maxVal <- max(apply(change.matrix,1,function(x) max(na.omit(x))))
	colors <- colorRampPalette(c("white", "orange"))(5)
	col.map <- function(x) {
		if (x < 0) "gray"
		else if (x < 4) colors[[x+1]]
		else if (x >= 4) colors[[5]]
	}

	### Draw the diagram
	# use horizontal axis labels
	op <- c(op,par(las=1))
	par(mar=c(5.1,4.1,0,2.1),xaxs="i")
	# create an empty plot
	plot(0,
		type='n',
		axes=FALSE,
		xlim=c(0,ncol(change.matrix)), 
		ylim=c(0,21),
		xlab="Position",
		ylab="Amino acid"
	)
	# iterate over each matrix entry and draw the contents on the plot
	for (x in 1:ncol(change.matrix)) {
		for (y in 1:21) {
			rect(x-1,22-y,x,21-y,col=col.map(change.matrix[y,x]), lty="blank")
			# if (change.matrix[y,x] > 0) {
			# 	#observed mutations are drawn in a color shade corresponding to their count
			# 	col <- colors[ceiling(4*change.matrix[y,x]/maxVal)+1]
			# 	rect(x-1,22-y,x,21-y,col=col, lty="blank")
			# } else if (change.matrix[y,x] < 0) {
			# 	rect(x-1,22-y,x,21-y,col="gray", lty="blank")
			# }
		}
	}
	# draw axes
	axis(1, at=c(1,seq(5,ncol(change.matrix),5))-.5, labels=c(1,seq(5,ncol(change.matrix),5)))
	axis(2, at=(1:21)-.5, labels=rev(rownames(change.matrix)))

	par(op)

	op <- par(mar=c(5.1,0,0,4.1),las=1)
	plot(0,type="n",ylim=c(0,6),xlim=c(0,1),axes=FALSE,ylab="",xlab="")
	# tops <- 1:4 * maxVal/4
	# bottoms <- round(tops - maxVal/4 + 1)
	# tops <- round(tops)
	rect(0,1:4,1,2:5,col=colors[2:5],lty="blank")
	rect(0,5,1,6,col="gray",lty="blank")
	# axis(4,at=0:5+.5,labels=c("0",paste(bottoms,tops,sep="-"),"wt"),tick=FALSE)
	axis(4,at=0:5+.5,labels=c(0:3,paste(4,maxVal,sep="-"),"wt"),tick=FALSE)
}


ref.seqs <- lapply(
	c("ube2i-ad.fa","ube2i-comp.fa","ube2i-db.fa"), 
	# rep("ube2i-pdonr.fa",6),
	function(fn) {
		con <- file(fn,open="r")
		ref.fa <- readFASTA(con)[[1]]
		close(con)
		ref.fa
	}
)

ref.beds <- do.call(rbind,lapply(
	c("ube2i-ad.bed","ube2i-comp.bed","ube2i-db.bed"), 
	# rep("ube2i-pdonr.bed",6),
	function(fn) {
		ref.bed <- read.delim(fn,header=FALSE)
		colnames(ref.bed) <- c("ref","from","to","feature")
		ref.bed
	}
))

clones <- lapply(c(
	"robot_ube2i-ad/clone_table.csv",
	"robot_ube2i-comp/clone_table.csv",
	"robot_ube2i-db/clone_table.csv"
	),read.csv,stringsAsFactors=FALSE
)

setnames <- c("UBE2I-AD","UBE2I-Complementation","UBE2I-DB")
nclones <- sapply(clones,nrow)

pdf("rearrayed_coverage_rescaled.pdf",11,5)
for (i in 1:3) {
	plotMutCoverage(clones[[i]],ref.seqs[[i]],ref.beds[i,],paste(setnames[[i]],":",nclones[[i]],"clones"))
}
dev.off()

filtered.clones <- lapply(clones, function(set) set[set$freq > .6,])
nclones <- sapply(filtered.clones,nrow)

pdf("rearrayed_coverage_safe.pdf",11,5)
for (i in 1:3) {
	plotMutCoverage(filtered.clones[[i]],ref.seqs[[i]],ref.beds[i,],paste(setnames[[i]],":",nclones[[i]],"clones"))
}
dev.off()
