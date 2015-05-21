
# calls <- rbind(
# 	cbind(set="sumo1gg-bplr-db",read.csv("sumo1gg-bplr-db/calls.csv",stringsAsFactors=FALSE)),
# 	cbind(set="ube2i-bplr-ad",read.csv("ube2i-bplr-ad/calls.csv",stringsAsFactors=FALSE)),
# 	cbind(set="ube2i-bplr-db",read.csv("ube2i-bplr-db/calls.csv",stringsAsFactors=FALSE)),
# 	cbind(set="ube2i-emlr-ad",read.csv("ube2i-emlr-ad/calls.csv",stringsAsFactors=FALSE)),
# 	cbind(set="ube2i-emlr-comp",read.csv("ube2i-emlr-comp/calls.csv",stringsAsFactors=FALSE))
# )
# calls$call[which(regexpr("NA",calls$call)>0)] <- "WT"
# write.table(calls,"all_calls.csv",sep=",",row.names=FALSE)


library("hash")
options(stringsAsFactors=FALSE)

# source("lib/libyogitools.R") #Helper functions
# source("lib/libyogiseq.R")	#FASTA parser

source("../src/main/R/libyogitools.R") #Helper functions
source("../src/main/R/libyogiseq.R")	#FASTA parser

calls <- read.csv("all_calls.csv")
top.calls <- calls[tapply(1:nrow(calls),paste(calls$plate,calls$well,sep="-"),min),]

set.names <- unique(top.calls$set)

ref.seqs <- lapply(
	c("sumo1gg-db.fa","ube2i-ad.fa","ube2i-db.fa","ube2i-ad.fa","ube2i-comp.fa"), 
	function(fn) {
		con <- file(fn,open="r")
		ref.fa <- readFASTA(con)[[1]]
		close(con)
		ref.fa
	}
)

ref.beds <- do.call(rbind,lapply(
	c("sumo1gg-db.bed","ube2i-ad.bed","ube2i-db.bed","ube2i-ad.bed","ube2i-comp.bed"), 
	function(fn) {
		ref.bed <- read.delim(fn,header=FALSE)
		colnames(ref.bed) <- c("ref","from","to","feature")
		ref.bed
	}
))



###
#DEAL WITH DUPLICATED BARCODES
#
bc.map <- hash()
for (i in 1:nrow(top.calls)) {
	bc <- top.calls$seq[[i]]
	if (is.na(bc) || bc=="") {
		next
	}
	bc.map[[bc]] <- c(bc.map[[bc]],i)
}
dup.bc <- keys(bc.map)[sapply(keys(bc.map),function(k)length(bc.map[[k]]))>1]
dup.bc.idxs <- values(bc.map,dup.bc)
#For each duplicated barcode identify the most useful clone and only retain
#that one
invisible(lapply(dup.bc.idxs,function(idxs) {
	if (length(idxs) > 6) {#BC crossover
		top.calls$call[idxs] <<- "Duplicated BC"
		return()
	} 
	bad <- regexpr("No|Low|WT|Corr|Cross",top.calls[idxs,"call"])>0
	good.idxs <- idxs[!bad]
	bad.idxs <- idxs[bad]
	if (length(good.idxs) == 0) {
		top.calls$call[idxs] <<- "Duplicated BC"
		return()
	} 
	if (length(good.idxs) > 1) {
		good.idxs <- good.idxs[order(nchar(top.calls[good.idxs,"call"]),decreasing=TRUE)]
		bad.idxs <- c(bad.idxs,good.idxs[2:length(good.idxs)])
	}
	top.calls$call[bad.idxs] <<- "Duplicated BC"
}))

###
#DEAL WITH PCR CROSSOVER
#
for (set.name in set.names) {

	set.idxs <- which(top.calls$set==set.name)

	muts <- global.extract.groups(top.calls$call[set.idxs],"([ACTG]\\d+[ACTG*])")
	muts <- lapply(muts,function(x)x[,1])
	mut.map <- hash()
	for (i in 1:length(muts)) {
		for (j in 1:length(muts[[i]])) {
			m <- muts[[i]][[j]]
			if (m=="") {
				next
			}
			mut.map[[m]] <- c(mut.map[[m]],set.idxs[[i]])
		}
	}
	ffs <- keys(mut.map)[sapply(keys(mut.map),function(k)length(mut.map[[k]])>20)]
	overlap.idxs <- Reduce(union,lapply(1:length(ffs), function(i) {
		do.call(c,lapply(1:length(ffs), function(j) {
			if (j > i) {
				intersection <- intersect(mut.map[[ffs[[i]]]],mut.map[[ffs[[j]]]])
				if (length(intersection) > 1) {
					intersection
				}
			} else NULL
		}))
	}))
	top.calls$call[overlap.idxs] <- paste("X-Over:",top.calls$call[overlap.idxs])

}


translate.calls <- function(.calls, ref.fa, ref.bed, aa.table="codontable.txt") {

	tr <- init.translator(aa.table)

	codon.starts <- seq(ref.bed$from,ref.bed$to,3)
	to.codon.idx <- function(i) {
		if (i < ref.bed$from || i > ref.bed$to) {
			c(codon.idx=NA,pos.in.codon=NA)
		} else {
			cidx <- min(which(c(codon.starts,ref.bed$to+1) > i))-1
			c(codon.idx=cidx,pos.in.codon=i-codon.starts[[cidx]]+1)
		}
	}
	muts <- strsplit(.calls$call,":|,")
	lapply(muts, function(ms) {
		if (length(ms)==0) {
			return(NA)
		}
		if (is.na(ms) || regexpr("No|Low|WT|Corr|X-Over|Dupl|Frame",ms[[1]])>0) {
			return(ms[[1]])
		}
		per.codon <- to.df(do.call(rbind,lapply(ms,function(m){
			from <- substr(m,1,1)
			pos <- as.numeric(substr(m,2,nchar(m)-1))
			to <- substr(m,nchar(m),nchar(m))
			codon.pos <- to.codon.idx(pos)
			list(cidx=codon.pos[[1]],cpos=codon.pos[[2]],to=to)
		})))
		if (any(is.na(per.codon$cidx))) {
			return("Corrupt!")
		}
		aa.muts <- tapply(1:nrow(per.codon),per.codon$cidx,function(rows) {
			cidx <- per.codon$cidx[[rows[[1]]]]
			cstart <- codon.starts[[cidx]]
			from.codon <- subseq(ref.fa,cstart,cstart+2)$toString()
			to.codon <- from.codon
			for (row in rows) {
				pos <- per.codon$cpos[[row]]
				substr(to.codon,pos,pos) <- per.codon$to[[row]]
			}
			from.aa <- tr$translate(from.codon)
			to.aa <- tr$translate(to.codon)
			if (from.aa==to.aa){
				return("silent")
			} else if (to.aa=="*"){
				return("Stop!")
			} else {
				paste(from.aa,cidx,to.aa,sep="")
			}
		})
		if (any(aa.muts=="Stop!")) {
			return("Stop!")
		} else {
			non.silent <- aa.muts[aa.muts!="silent"]
			if (length(non.silent)==0) {
				return("WT")
			}
			non.silent
		}
	})
}

top.calls$aa.calls <- NA
for (set.idx in 1:length(set.names)) {
	set.name <- set.names[[set.idx]]
	ref.seq <- ref.seqs[[set.idx]]
	ref.bed <- ref.beds[set.idx,]
	idxs <- which(top.calls$set==set.name)
	.calls <- top.calls[idxs,]
	tr <- translate.calls(.calls,ref.seq,ref.bed)
	tr.str <- sapply(tr,paste,collapse=",")
	top.calls[idxs,"aa.calls"] <- tr.str
}
write.table(top.calls,"processed_calls.csv",sep=",",row.names=FALSE)

###
# CALCULATE NUCLEOTIDE-BASED CENSUS OF MUTATIONS
###
pop.vs.snp <- function(.calls, ref.fa, ref.bed) {
	codon.starts <- seq(ref.bed$from,ref.bed$to,3)
	to.codon.idx <- function(i) {
		if (i < ref.bed$from || i > ref.bed$to) {
			c(codon.idx=NA,pos.in.codon=NA)
		} else {
			cidx <- min(which(c(codon.starts,ref.bed$to+1) > i))-1
			c(codon.idx=cidx,pos.in.codon=i-codon.starts[[cidx]]+1)
		}
	}
	muts <- strsplit(.calls$call,":|,")
	pvs <- to.df(do.call(rbind,lapply(muts, function(ms) {

		if (regexpr("No|Low|WT|Corr|X-Over|Dupl|Frame",ms[[1]])>0) {
			return(list(pop=NA,snp=NA,comment=ms[[1]]))
		}
		per.codon <- to.df(do.call(rbind,lapply(ms,function(m){
			from <- substr(m,1,1)
			pos <- as.numeric(substr(m,2,nchar(m)-1))
			to <- substr(m,nchar(m),nchar(m))
			codon.pos <- to.codon.idx(pos)
			list(cidx=codon.pos[[1]],cpos=codon.pos[[2]],to=to)
		})))
		c(table(factor(
			sapply(
				tapply(1:nrow(per.codon),per.codon$cidx,length),
				function(x) if(x>1)"pop"else"snp"
			),
			levels=c("pop","snp")
		)),list(comment=""))
	})))
}

calc.nc.census <- function(.calls,ref.fa,ref.bed) {
	muts <- pop.vs.snp(.calls, ref.fa, ref.bed)
	labels <- c("No alignment!","Low coverage!","X-Over","Duplicated BC","Corrupt","Frameshift","Stop","WT",1:14)
	counts <- c(rep(0,8),rep(list(list(pop=0,snp=0)),14))
	names(counts) <- labels
	for (i in 1:nrow(muts)) {
		if (is.na(muts$pop[[i]])) {
			x <- muts$comment[[i]]
			counts[[x]] <- counts[[x]]+1
		} else {
			n <- sum(muts[i,c("pop","snp")])
			if (n == 0) {
				counts[["WT"]] <- counts[["WT"]]+1
				next
			}
			x <- as.character(n)
			counts[[x]]$pop <- counts[[x]]$pop+muts[i,"pop"]/n
			counts[[x]]$snp <- counts[[x]]$snp+muts[i,"snp"]/n
		}
	}
	rowsum <- sum(unlist(counts))
	out <- do.call(cbind,lapply(counts,function(x) if(class(x)=="list") lapply(x,`/`,rowsum) else list(pop=0,snp=x/rowsum)))
	colnames(out) <- labels
	out
}

nc.census <- to.df(do.call(rbind,lapply(1:length(set.names), function(set.idx) {
	set.name <- set.names[[set.idx]]
	ref.seq <- ref.seqs[[set.idx]]
	ref.bed <- ref.beds[set.idx,]
	.calls <- top.calls[top.calls$set==set.name,]
	calc.nc.census(.calls,ref.seq,ref.bed)
})))
colnames(nc.census) <- c("No alignment!","Low coverage!","X-Over","Duplicated BC","Corrupt","Frameshift","Stop","WT",1:14)


draw.census <- function(census,filename,nc=FALSE) {

	pdf(filename,11,8.5)
	draw <- function(i) {
		ridx <- c(2*i-1,2*i)
		max.y <- max(apply(census[ridx,],2,sum))
		data <- as.matrix(census[ridx,])
		data.num <- data
		data.err <- data[2,]
		if (nc) {
			data.num <- data.num[,-7]
			data.err <- data.err[-7]
		}
		first <- if (nc) 8 else 9
		last <- ncol(data.num)
		succ.rate <- sum(data.num[,first:last])/(1-sum(data.err[1:3]))
		good.clones <- sum(data.num[,first:last]) * sum(top.calls$set==set.names[[i]])
		data.err[first:last] <- 0
		data.num[,1:(first-1)] <- 0
		xs <- barplot(
			data.num,
			ylim=c(0,if(max.y > .6) 1 else .6),
			col=c("darkolivegreen4","darkolivegreen3"),
			border=NA,
			main=set.names[[i]],
			ylab="rel. freq.",
			xlab=if(nc) "Codon changes" else "AA changes"
		)
		barplot(
			data.err,
			add=TRUE,
			border=NA,
			col=c(rep("gray",3),rep("firebrick3",ifelse(nc,4,5)),rep("darkolivegreen3",14))
		)
		if (i==2) {
			legend.labels <- if (nc) {
				c("SNV","POP","No data","Unusable")
			} else {
				c("SNP-accessible","SNP-inaccessible","No data","Unusable")
			}
			legend("topright",
				legend.labels,
				fill=c("darkolivegreen3","darkolivegreen4","gray","firebrick3"),
				border=NA
			)
		}
		y <- max(apply(census[ridx,first:last],2,sum))
		# last <- if (nc) 15 else 16
		clip(xs[1]-2,xs[last]+2,0,1.2)
		arrows(xs[first],y,xs[last],length=.02,angle=90,code=3)
		text(mean(xs[c(first,last)]),y,paste(good.clones,"=",signif(succ.rate*100,4),"%"),pos=3)
	}
	op <- par(mfrow=c(2,2),las=3,mar=c(8,4,4,2)+.1)
	for (i in 1:length(set.names)) {
		draw(i)
	}
	par(op)
	dev.off()
}

draw.census(nc.census,"nc-census_all.pdf",nc=TRUE)

####
## AMINO-ACID CENSUS:
####

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

acc.vs.inacc <- function(.calls, ref.fa, ref.bed) {
	reachables <- calc.reachables(ref.bed,ref.fa)
	# muts <- translate.calls(.calls,ref.fa,ref.bed)
	muts <- strsplit(.calls$aa.calls,",")
	to.df(do.call(rbind,lapply(muts, function(ms) {
		if (length(ms) == 1) {
			if (regexpr("No|Low|WT|Corr|X-Over|Dupl|Frame|Stop",ms[[1]])>0) {
			# if (ms %in% c("No alignment!","Low coverage!","Crossover!","WT","Frameshift!","Corrupt!","Stop!")) {
				return(list(pop=NA,snp=NA,comment=ms))
			} else if ("share" %in% colnames(.calls) && .calls[i,"share"] < .25) {
				return(list(pop=0,snp=0,comment="WT"))
			} 
		}
		counts <- table(factor(sapply(ms,function(m) {
			from <- substr(m,1,1)
			pos <- as.numeric(substr(m,2,nchar(m)-1))
			to <- substr(m,nchar(m),nchar(m))
			if (to %in% reachables[[pos]]) {
				"snp"
			} else {
				"pop"
			}
		}),levels=c("pop","snp")))
		c(counts,list(comment=""))
	})))	
}

calc.aa.census <- function(.calls,ref.fa,ref.bed) {
	muts <- acc.vs.inacc(.calls, ref.fa, ref.bed)
	labels <- c("No alignment!","Low coverage!","X-Over","Duplicated BC","Corrupt","Frameshift","Stop!","WT",1:14)
	counts <- c(rep(0,8),rep(list(list(pop=0,snp=0)),14))
	names(counts) <- labels
	for (i in 1:nrow(muts)) {
		if (is.na(muts$pop[[i]])) {
			x <- muts$comment[[i]]
			counts[[x]] <- counts[[x]]+1
			
		} else {
			n <- sum(muts[i,c("pop","snp")])
			x <- as.character(n)
			counts[[x]]$pop <- counts[[x]]$pop+muts[i,"pop"]/n
			counts[[x]]$snp <- counts[[x]]$snp+muts[i,"snp"]/n
		}
	}
	rowsum <- sum(unlist(counts))
	out <- do.call(cbind,lapply(counts,function(x) if(class(x)=="list") lapply(x,`/`,rowsum) else list(pop=0,snp=x/rowsum)))
	colnames(out) <- labels
	out
}


aa.census <- to.df(do.call(rbind,lapply(1:length(set.names), function(set.idx) {
	set.name <- set.names[[set.idx]]
	ref.seq <- ref.seqs[[set.idx]]
	ref.bed <- ref.beds[set.idx,]
	.calls <- top.calls[top.calls$set==set.name,]
	calc.aa.census(.calls,ref.seq,ref.bed)
})))
colnames(aa.census) <- c("No alignment!","Low coverage!","X-Over","Duplicated BC","Corrupt","Frameshift","Stop","WT",1:14)

draw.census(aa.census,"aa-census_all.pdf",nc=FALSE)



####
# DETERMINE POSITIONAL COVERAGE
####

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
			if (change.matrix[y,x] > 0) {
				#observed mutations are drawn in a color shade corresponding to their count
				col <- colors[ceiling(4*change.matrix[y,x]/maxVal)+1]
				rect(x-1,22-y,x,21-y,col=col, lty="blank")
			} else if (change.matrix[y,x] < 0) {
				rect(x-1,22-y,x,21-y,col="gray", lty="blank")
			}
		}
	}
	# draw axes
	axis(1, at=c(1,seq(5,ncol(change.matrix),5))-.5, labels=c(1,seq(5,ncol(change.matrix),5)))
	axis(2, at=(1:21)-.5, labels=rev(rownames(change.matrix)))

	par(op)

	op <- par(mar=c(5.1,0,0,4.1),las=1)
	plot(0,type="n",ylim=c(0,6),xlim=c(0,1),axes=FALSE,ylab="",xlab="")
	tops <- 1:4 * maxVal/4
	bottoms <- round(tops - maxVal/4 + 1)
	tops <- round(tops)
	rect(0,1:4,1,2:5,col=colors[2:5],lty="blank")
	rect(0,5,1,6,col="gray",lty="blank")
	axis(4,at=0:5+.5,labels=c("0",paste(bottoms,tops,sep="-"),"wt"),tick=FALSE)
}

pdf("positionalCoverage.pdf",width=11,height=4)
lapply(1:length(set.names), function(set.idx) {
	set.name <- set.names[[set.idx]]
	ref.seq <- ref.seqs[[set.idx]]
	ref.bed <- ref.beds[set.idx,]
	.calls <- top.calls[top.calls$set==set.name,]
	plotMutCoverage(.calls,ref.seq,ref.bed,set.name)
})
dev.off()

