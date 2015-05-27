options(stringsAsFactors=FALSE)
library("hash")

#Load calls
calls <- read.csv("processed_calls.csv")

#Filter down to only UBE2I-EMLR-AD and UBE2I-BPLR-AD sets
# calls <- calls[calls$set %in% c("ube2i-emlr-ad","ube2i-bplr-ad"),]
calls <- calls[calls$set %in% c("ube2i-emlr-comp","ube2i-bplr-comp"),]
# calls <- calls[calls$set == "ube2i-bplr-db",]
#Filter out unusable clones
usable.calls <- calls[regexpr("\\w\\d+\\w",calls$aa.calls)>0 & nchar(calls$seq)==25,]

#index mutations
mut.index <- hash()
for (i in 1:nrow(usable.calls)) {
	mut <- usable.calls$aa.calls[[i]]
	mut.index[[mut]] <- c(mut.index[[mut]],i)
}

#Filter down to up to 3 replicates per mutant
capped.calls <- usable.calls[unlist(lapply(values(mut.index),head,3)),]
#Restore proper sorting order
capped.calls <- capped.calls[order(with(capped.calls,paste(set,plate,well,sep="_"))),]

#number of complete destination plates (not including incomplete plate)
nplates <- nrow(capped.calls) %/% 384
#number of remaining clones that will go on the incomplete plate
nrest <- nrow(capped.calls) %% 384
#list of destination plates
dest.plate <- c(do.call(c,lapply(1:nplates,rep,384)),rep(nplates+1,nrest))

#table of rows and column coordinates on one plate
row.col <- expand.grid(col=1:24,row=1:16)
#list of corresponding well names
wells <- apply(row.col,1,function(x)paste(LETTERS[x[["row"]]],sprintf("%02d",x[["col"]]),sep=""))

#list of well names for each clone
dest.well <- c(rep(wells,nplates),wells[1:nrest])
#table of well coordinates for each clone
dest.well.coord <- rbind(do.call(rbind,replicate(nplates,row.col,simplify=FALSE)),row.col[1:nrest,])

#List of unique IDs for each clone
# clone.id <- paste("UBE2I-AD",sprintf("%05d",1:nrow(capped.calls)),sep="-")
# clone.id <- paste("UBE2I-DB",sprintf("%05d",1:nrow(capped.calls)),sep="-")
clone.id <- paste("UBE2I-HYC",sprintf("%05d",1:nrow(capped.calls)),sep="-")

capped.calls$id <- clone.id
capped.calls$dest.plate <- dest.plate
capped.calls$dest.well <- dest.well

#list of source plates
from.plate <- with(capped.calls,paste(set,substr(plate,6,7),sep="-"))
#list of source plate well coordinates
from.well.coord <- do.call(rbind,lapply(capped.calls$well,function(x) {
	q <- which(LETTERS==substr(x,1,1))
	r <- which(LETTERS==substr(x,3,3))
	c <- as.numeric(substr(x,4,nchar(x)))
	c(row=r*2 - (q<3),col=c*2 - q%%2)
}))

#compute source plate files for robot
source.tables <- tapply(1:nrow(capped.calls),from.plate,function(idxs) {
	data.frame(
		c=from.well.coord[idxs,"col"],
		r=from.well.coord[idxs,"row"],
		Gene=clone.id[idxs]
	)
})

#compute destination plate files for robot
dest.tables <- tapply(1:nrow(capped.calls),dest.plate,function(idxs) {
	data.frame(
		c=dest.well.coord[idxs,"col"],
		r=dest.well.coord[idxs,"row"],
		Gene=clone.id[idxs]
	)
})

sources.per.dest <- tapply(1:nrow(capped.calls),dest.plate,function(idxs) {
	unique(from.plate[idxs])
})

header.384 <- "ID,902
Name,384gl
Type,384
Category,Gel Plate with Lid"


# dir.name <- "robot/"
# dir.name <- "robot_ube2i-db/"
dir.name <- "robot_ube2i-comp/"
dir.create(dir.name)

write.table(capped.calls,paste(dir.name,"clone_table.csv",sep=""),sep=",",row.names=FALSE)

for (st in names(source.tables)) {
	outfile <- file(paste(dir.name,st,".csv",sep=""),open="a")
	writeLines(header.384,outfile)
	write.table(source.tables[[st]],outfile,sep=",",row.names=FALSE)
	close(outfile)
}

for (dt in names(dest.tables)) {
	outfile <- file(paste(dir.name,"dest_",sprintf("%02d",as.numeric(dt)),".csv",sep=""),open="a")
	writeLines(header.384,outfile)
	write.table(dest.tables[[dt]],outfile,sep=",",row.names=FALSE)
	close(outfile)
}

outfile <- file(paste(dir.name,"sourcesPerDest.txt",sep=""),open="w")
for (dest in names(sources.per.dest)) {
	cat("dest",dest,":",paste(sources.per.dest[[dest]],collapse=", "),"\n",file=outfile)
}
close(outfile)

