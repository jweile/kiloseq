#!$Rbin

options(stringsAsFactors=FALSE)

source("lib/libyogitools.R") #Helper functions
source("lib/libyogiseq.R")	#FASTA parser

con <- file("ube2i-pdonr.fa",open="r")
ref.fa <- readFASTA(con)[[1]]
close(con)

ref.bed <- read.delim("ube2i-pdonr.bed",header=FALSE)
colnames(ref.bed) <- c("ref","from","to","feature")


calls <- read.csv("calls.csv")

joint.calls <- do.call(rbind,tapply(1:nrow(calls),paste(calls$set,calls$well,sep="-"),function(is) {
	common <- calls[min(is),1:6]
	if (any(regexpr("No alignment!",calls[is,"call"]) > 0)) {
		return(cbind(common,call="No alignment!"))
	}
	if (any(regexpr("Low coverage!",calls[is,"call"]) > 0)) {
		return(cbind(common,call="Low coverage!"))
	}
	if (all(regexpr("WT",calls[is,"call"]) > 0)) {
		return(cbind(common,call="WT"))
	}
	top.is <- is[which(calls[is,"share"] > .5)]
	if (length(top.is)==0) {
		return(cbind(common,call="WT"))
	}
	top.calls <- calls[top.is,"call"]
	if (any(regexpr("Frameshift",top.calls) > 0)) {
		return(cbind(common,call="Frameshift!"))
	}
	if (any(regexpr("Corrupt",top.calls) > 0)) {
		return(cbind(common,call="Corrupt!"))
	}
	muts <- unlist(global.extract.groups(top.calls,"([ACGT]\\d+[ACGT\\*-\\+])"))
	return(cbind(common,call=paste(muts,collapse=",")))

},simplify=FALSE))

top.calls <- calls[tapply(1:nrow(calls),paste(calls$set,calls$well,sep="-"),min),]

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
	muts <- strsplit(.calls$call,",")
	pvs <- to.df(do.call(rbind,lapply(muts, function(ms) {
		if (length(ms)==1) {
			if (ms %in% c("No alignment!","Low coverage!","WT","Frameshift!","Corrupt!")) {
				return(list(pop=NA,snp=NA,comment=ms))
			}
		}
		if (regexpr("Frameshift",ms[[1]])>0) {
			return(list(pop=NA,snp=NA,comment="Frameshift!"))
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
	# muts <- strsplit(.calls$call,",")
	# muts <- translate.calls(.calls,ref.fa,ref.bed)
	muts <- pop.vs.snp(.calls, ref.fa, ref.bed)
	subsets <- sapply(1:nrow(.calls),function(i) {
		if (regexpr("Polymerase",.calls[i,"set"]) > 0) {
			sub("_R1",substr(.calls[i,"well"],1,1),.calls[i,"set"])
		} else {
			sub("_R1","",.calls[i,"set"])
		}
	})
	labels <- c("No alignment!","Low coverage!","Corrupt!","Frameshift!","Stop!","WT",1:10)
	census <- do.call(rbind,tapply(1:length(subsets),subsets,function(is) {
		counts <- c(rep(0,6),rep(list(list(pop=0,snp=0)),10))
		names(counts) <- labels
		for (i in is) {
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
		counts
	}))
	out <- to.df(do.call(rbind,lapply(1:nrow(census),function(i){
		row <- census[i,]
		rowsum <- sum(unlist(row))
		do.call(cbind,lapply(row,function(x) if(class(x)=="list") lapply(x,`/`,rowsum) else list(pop=0,snp=x/rowsum)))
	})))
	colnames(out) <- labels
	out
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
	muts <- strsplit(.calls$call,",")
	lapply(muts, function(ms) {
		if (length(ms)==1) {
			if (ms %in% c("No alignment!","Low coverage!","WT","Frameshift!","Corrupt!")) {
				return(ms)
			}
		}
		if (regexpr("Frameshift",ms[[1]])>0) {
			return("Frameshift!")
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

# tr.joint <- translate.calls(joint.calls,ref.fa,ref.bed)
# tr.top <- translate.calls(top.calls,ref.fa,ref.bed)

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
	muts <- translate.calls(.calls,ref.fa,ref.bed)
	to.df(do.call(rbind,lapply(muts, function(ms) {
		if (length(ms) == 1) {
			if (ms %in% c("No alignment!","Low coverage!","WT","Frameshift!","Corrupt!","Stop!")) {
				return(list(pop=NA,snp=NA,comment=ms))
			} else if ("share" %in% colnames(.calls) && .calls[i,"share"] < .5) {
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
	# muts <- strsplit(.calls$call,",")
	# muts <- translate.calls(.calls,ref.fa,ref.bed)
	muts <- acc.vs.inacc(.calls, ref.fa, ref.bed)
	subsets <- sapply(1:nrow(.calls),function(i) {
		if (regexpr("Polymerase",.calls[i,"set"]) > 0) {
			sub("_R1",substr(.calls[i,"well"],1,1),.calls[i,"set"])
		} else {
			sub("_R1","",.calls[i,"set"])
		}
	})
	labels <- c("No alignment!","Low coverage!","Corrupt!","Frameshift!","Stop!","WT",1:10)
	census <- do.call(rbind,tapply(1:nrow(.calls),subsets,function(is) {
		counts <- c(rep(0,6),rep(list(list(pop=0,snp=0)),10))
		names(counts) <- labels
		for (i in is) {
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
		counts
	}))
	out <- to.df(do.call(rbind,lapply(1:nrow(census),function(i){
		row <- census[i,]
		rowsum <- sum(unlist(row))
		do.call(cbind,lapply(row,function(x) if(class(x)=="list") lapply(x,`/`,rowsum) else list(pop=0,snp=x/rowsum)))
	})))
	colnames(out) <- labels
	out
}

draw.census <- function(census,filename,nc=FALSE) {

	labels <- c(
		"Taq buffer stepwise","Taq buffer snapback",
		"H2O stepwise","H2O snapback",
		"Annealing buffer stepwise","Annealing buffer snapback",
		"Sulfolobus","Klenow", "Klenow Exo",
		"T7", "T4", "2 rounds","3 rounds", "Q5",
		"M13F"
	)

	pdf(filename,11,8.5)
	draw <- function(i) {
		ridx <- c(2*i-1,2*i)
		succ.rate <- sum(census[ridx,7:16])/(1-sum(census[ridx,1:2]))
		max.y <- max(apply(census[ridx,],2,sum))
		data <- as.matrix(census[ridx,])
		data.num <- data
		data.num[,1:6] <- 0
		data.err <- data[2,]
		data.err[7:16] <- 0
		if (nc) {
			data.num <- data.num[,-5]
			data.err <- data.err[-5]
		}
		xs <- barplot(
			data.num,
			ylim=c(0,if(max.y > .6) 1 else .6),
			col=c("darkolivegreen4","darkolivegreen3"),
			border=NA,
			main=labels[[i]],
			ylab="rel. freq.",
			xlab="#AA changes"
		)
		barplot(
			data.err,
			add=TRUE,
			border=NA,
			col=c(rep("gray",2),rep("firebrick3",4),rep("darkolivegreen3",10))
		)
		if (i==3) {
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
		y <- max(apply(census[ridx,7:16],2,sum))
		last <- if (nc) 15 else 16
		first <- if (nc) 6 else 7
		clip(xs[1]-2,xs[last]+2,0,1.2)
		arrows(xs[first],y,xs[last],length=.02,angle=90,code=3)
		text(mean(xs[c(first,last)]),y,paste(signif(succ.rate*100,4),"%"),pos=3)
	}
	op <- par(mfrow=c(2,3),las=3,mar=c(8,4,4,2))
	for (i in 1:6) {
		draw(i)
	}
	par(op)
	op <- par(mfrow=c(3,4),las=3,mar=c(8,4,4,2))
	for (i in 6+1:8) {
		draw(i)
	}
	draw(15)
	par(op)
	dev.off()
}

nc.census <- calc.nc.census(joint.calls,ref.fa,ref.bed)
# top.census <- calc.aa.census(top.calls,ref.fa,ref.bed)
aa.census <- calc.aa.census(joint.calls,ref.fa,ref.bed)

draw.census(nc.census,"census_nc.pdf",nc=TRUE)
draw.census(aa.census,"census_aa.pdf")





