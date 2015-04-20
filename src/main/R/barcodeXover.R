library("hash")

to.df <- function(x) {
	y <- lapply(1:ncol(x), function(col) {
		unlist(x[,col])
	})
	names(y) <- colnames(x)
	as.data.frame(y,stringsAsFactors=FALSE)
}

dirs <- list.dirs(recursive=FALSE)
dirs <- dirs[grep("ube2i",dirs)]
calls <- do.call(rbind,lapply(dirs, function(dir) {
	read.csv(paste(dir,"/calls.csv",sep=""),stringsAsFactors=FALSE)
}))

bc.idx <- hash()
for (i in 1:nrow(calls)) {
	bc <- calls[i,"seq"]
	if (nchar(bc) == 0) {
		next
	}
	bc.idx[[bc]] <- c(bc.idx[[bc]],i)
}

bc.reps <- sapply(keys(bc.idx),function(k) length(bc.idx[[k]]))


library("beeswarm")
example.bc <- names(which(bc.reps==10))
example.bc <- example.bc[nchar(example.bc)==25][1:10]
example.freqs <- lapply(example.bc,function(bc){
	calls[bc.idx[[bc]],"freq"]
})
example.plates <-lapply(example.bc,function(bc){
	as.numeric(factor(calls[bc.idx[[bc]],"set"]))
})
op <- par(las=1,mar=c(5,20,4,2)+.1)
beeswarm(example.freqs,
	dlab="Barcode Frequency",glab="",
	pch=16,pwcol=example.plates,
	labels=example.bc,
	dlim=c(0,1),
	vertical=FALSE
)
beeswarm(example.freqs,
	add=TRUE,
	vertical=FALSE
)
grid(NULL,NA)
abline(h=0:10+.5,col="gray")
par(op)


rep.stats <- table(bc.reps)
barplot(rep.stats,xlab="#wells for barcode",ylab="frequency",log="y")

nsets <- length(unique(calls$set))
pos.enr <- do.call(rbind,lapply(keys(bc.idx), function(k) {
	rows <- bc.idx[[k]]
	wells <- calls[rows,"well"]
	reps <- table(wells)
	candidate.wells <- names(reps[reps > 1])
	do.call(rbind,lapply(candidate.wells, function(well) {
		x <- reps[[well]]
		mat <- rbind(
			bc=c(well=x,nwell=length(rows)-x),
			nbc=c(well=nsets-x,nwell=384*nsets-(length(rows)-x))
		)
		p <- fisher.test(mat,alternative="greater")$p.value
		if (p < .05) {
			list(bc=k,well=well,rep=x,p=p)
		} else {
			NULL
		}
	}))
}))
pos.enr <- to.df(pos.enr)
pos.enr <- pos.enr[order(pos.enr$rep,decreasing=TRUE),]
write.table(pos.enr,"bc_well_enrichment.csv",sep=",",row.names=FALSE)

plate.counts <- table(calls$set)
plate.enr <- do.call(rbind,lapply(keys(bc.idx), function(k) {
	rows <- bc.idx[[k]]
	plates <- calls[rows,"set"]
	reps <- table(plates)
	candidate.plates <- names(reps[reps > 1])
	do.call(rbind,lapply(candidate.plates, function(plate) {
		x <- reps[[plate]]
		mat <- rbind(
			bc=c(plate=x,nplate=length(rows)-x),
			nbc=c(plate=plate.counts[[plate]]-x, nplate=sum(plate.counts)-plate.counts[[plate]]-(length(rows)-x))
		)
		p <- fisher.test(mat,alternative="greater")$p.value
		if (p < .05) {
			list(bc=k,plate=plate,rep=x,p=p)
		} else {
			NULL
		}
	}))
}))
plate.enr <- to.df(plate.enr)
plate.enr <- plate.enr[order(plate.enr$rep,decreasing=TRUE),]
write.table(plate.enr,"bc_plate_enrichment.csv",sep=",",row.names=FALSE)


out <- do.call(rbind,lapply(keys(bc.idx), function(k) {
	rows <- bc.idx[[k]]

	wells <- calls[rows,"well"]
	wreps <- table(wells)
	candidate.wells <- names(wreps[wreps > 1])
	well.mask <- wells %in% candidate.wells

	plates <- calls[rows,"set"]
	preps <- table(plates)
	candidate.plates <- names(preps[preps > 1])
	plate.mask <- plates %in% candidate.plates

	sub.plates <- calls[rows,"set"][!(wells %in% candidate.wells)]

	c(
		both.dups=sum(well.mask & plate.mask),
		well.dups=sum(well.mask-plate.mask > 0),
		plate.dups=sum(plate.mask-well.mask > 0),
		all.dups=length(rows)
	)
}))

bd <- tapply(out[,"both.dups"],out[,"all.dups"],sum)
wd <- tapply(out[,"well.dups"],out[,"all.dups"],sum)
pd <- tapply(out[,"plate.dups"],out[,"all.dups"],sum)
ad <- tapply(out[,"all.dups"],out[,"all.dups"],sum)
freq <- table(out[,"all.dups"])
rd <- ad-(wd+pd+bd)

cols <- c("gray30","steelblue3","darkolivegreen3","darkgoldenrod1")
# op <- par(mfrow=c(2,1))
layout(rbind(1,2))
op <- par(mar=c(0,4,4,2)+.1)
barplot(
	rbind(rest=freq*rd/ad,plate=freq*pd/ad,both=freq*bd/ad,well=freq*wd/ad)+.000001,
	xlab="#samples with same barcode",
	ylab="Absolute frequency",
	col=cols,border="gray",log="y",ylim=c(1,100000)
)
usr <- par("usr")
clip(usr[1],usr[2],.0001,1000000)
rect(usr[1],.0001,usr[2],1,col="white",border=NA)
legend("right",c("random","plate effect","plate/well effect","well effect"),fill=cols)
par(op)
op <- par(mar=c(5,4,1,2)+.1)
barplot(
	rbind(rest=rd/ad,plate=pd/ad,both=bd/ad,well=wd/ad),
	xlab="#samples with same barcode",
	ylab="Relative frequency",
	col=cols,border="gray"
)
par(op)
