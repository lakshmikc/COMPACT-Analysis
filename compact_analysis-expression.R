input.data.dir <- '/Users/raj/Dropbox/research/liver-regeneration/timeseries'
results.save.dir <- '/Users/raj/Dropbox/research/liver-regeneration/timeseries'

input.data.file <- 'llmAvg-sub-from-batch-corrected-data.txt'


input.data.dir <- "/Users/raj/Dropbox/research/compact-analysis"
results.save.dir <- "/Users/raj/Dropbox/research/compact-analysis"


input.data.file <- "sample-data-1.txt"
comparative.condition <- "Diet"
pattern.condition     <- "TimePostPHx"
sample.type.row <- 2
comparative.condition.row <- 3
pattern.condition.row <- 4
data.column.start <- 9
data.row.start <- 6
header.row <- 1
discrete.levels <- c(cho=3,etoh=3,hf=3)


input.data.file <- "llmAvg-sub-from-batch-corrected-v2.txt"
comparative.condition <- "Diet"
pattern.condition     <- c("Treatment","TimePostPHx")
sample.type.row <- 2
comparative.condition.row <- 3
pattern.condition.row <- c(4,5)
data.column.start <- 22
data.row.start <- 7
data.columns <- c(26:29, 34:37, 42:44, 48:51, 55:57, 61:64, 68:71, 75:78, 81:84)
header.row <- 1
discrete.levels <- c(cho=3,etoh=3,hf=3)

compact.save.file <- "compact-matrix.txt"


#--------- Do Not Edit Below ---------

library(vegan)
library(gtools)
library(amap)

#setwd('/Users/raj/Documents/research/liver-regeneration/EtOH/timeseries')
#setwd('/Users/raj/Dropbox/research/liver-regeneration/timeseries')

setwd(results.save.dir)

#all.data <- read.table("llmAvg-sub-from-batch-corrected.txt",header=T,sep="\t",quote="\"",stringsAsFactors=F)
#data.column.offset <- 31

all.input.data <- read.input.data(input.data.file,input.data.dir,data.column.start,data.row.start,header.row,data.columns=data.columns)

sample.info <- read.sample.info(sample.info.data=all.input.data,sample.id.column=1,
								sample.type.column=sample.type.row,
								comparative.condition.column=comparative.condition.row,
								pattern.condition.column=pattern.condition.row)

average.per.group.per.gene <- calc.average.per.group.per.gene(input.data=all.input.data,sample.info=sample.info)

pattern.per.group.per.gene <- calc.pattern.discrete.per.group.per.gene(average.data=average.per.group.per.gene,
										sample.info=sample.info,threshold.value=log2(1.5))

patterns.counts.per.group <- calc.patterns.counts.per.group(pattern.data=pattern.per.group.per.gene,
												sample.info,discrete.levels)

compact.matrix <- calc.compact.matrix(pattern.data=pattern.per.group.per.gene,sample.info=sample.info,discrete.levels,symmetric=FALSE)

compact.matrix.symm <- calc.compact.matrix(pattern.data=pattern.per.group.per.gene,sample.info=sample.info,discrete.levels,symmetric=TRUE)

write.table(compact.matrix,file=compact.save.file,col.names=NA,sep="\t",quote=F)
write.table(compact.matrix.symm,file=compact.save.file,col.names=NA,sep="\t",quote=F)



# Optimization
test <- cbind(cho=patterns.counts.per.group$cho$patterns,etoh=patterns.counts.per.group$etoh$patterns)
test.res <- aggregate(x=test[,1],by=list(test[,1],test[,2]),FUN="length")
xtabs(x~.,data=test.res,sparse=T)
# Can also use table or ftable to generate compact matrix
table(patterns.counts.per.group$cho$patterns,patterns.counts.per.group$etoh$patterns)

test.res.2 <- ftable(patterns.counts.per.group$cho$patterns,patterns.counts.per.group$hf$patterns,patterns.counts.per.group$etoh$patterns,row.vars=c(1,2))


# Chi-squared test
# Can only do with built in function if no zero marginals, which might happen when forcing
#  to build a symmetric compact matrix.
# Some issue with the Chi-squared approximation being incorrect
#  if any of the random expected values are less than 5.
chisq.test(compact.matrix)$p.value



# Visualization
heatmap.2(compact.matrix.symm[[1]], Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=TRUE, scale="none", trace="none", breaks=c(0:50,800), col=c(colorpanel(n=50,low="white",high="red"),rep("#FF0000",1)), key=FALSE,colsep=1:ncol(compact.matrix.symm[[1]]),rowsep=1:nrow(compact.matrix.symm[[1]]),sepwidth=c(0.1,0.1))


color.lookup <- c(colorpanel(n=51,low="white",high="red"),rep("#FF0000",800))

library(rgl)
color.lookup <- c(colorpanel(n=51,low="white",high="red"),rep("#FF0000",2000))
cg1.indices <- 1:nrow(compact.matrix.symm[[1]])
cg2.indices <- 1:ncol(compact.matrix.symm[[1]])
z.data       <- compact.matrix.symm[[1]]
color.map    <- color.lookup[1+compact.matrix.symm[[1]]]

cg1.normals <- matrix(cg1.indices, nrow=length(cg1.indices), ncol=length(cg2.indices), byrow=F)
cg2.normals <- matrix(cg2.indices, nrow=length(cg1.indices), ncol=length(cg2.indices), byrow=T)

persp3d(cg1.indices, cg2.indices, z.data, col=color.map, normal_x=cg1.normals, normal_y=cg2.normals, normal_z=z.data, lit=F)


############
# Old code
############

all.data <- read.table(input.data.file.full.path,header=T,sep="\t",quote="\"",stringsAsFactors=F,row.names=1)

fc.thres <- 1.5

etoh.llm.1h.cols <- (1:4) + data.column.offset
etoh.phx.1h.cols <- (5:8) + data.column.offset
etoh.llm.6h.cols <- (9:12) + data.column.offset
etoh.phx.6h.cols <- (13:16) + data.column.offset
etoh.llm.24h.cols <- (17:20) + data.column.offset
etoh.phx.24h.cols <- (21:23) + data.column.offset

etoh.llm.cols <- c(etoh.llm.1h.cols,etoh.llm.6h.cols,etoh.llm.24h.cols)
etoh.phx.cols <- c(etoh.phx.1h.cols,etoh.phx.6h.cols,etoh.phx.24h.cols)
etoh.cols <- c(etoh.llm.1h.cols,etoh.phx.1h.cols,etoh.llm.6h.cols,etoh.phx.6h.cols,etoh.llm.24h.cols,etoh.phx.24h.cols)

hf.llm.1h.cols <- (24:26) + data.column.offset
hf.phx.1h.cols <- (27:30) + data.column.offset
hf.llm.6h.cols <- (31:33) + data.column.offset
hf.phx.6h.cols <- (34:36) + data.column.offset
hf.llm.24h.cols <- (37:39) + data.column.offset
hf.phx.24h.cols <- (40:43) + data.column.offset

hf.llm.cols <- c(hf.llm.1h.cols,hf.llm.6h.cols,hf.llm.24h.cols)
hf.phx.cols <- c(hf.phx.1h.cols,hf.phx.6h.cols,hf.phx.24h.cols)
hf.cols <- c(hf.llm.1h.cols,hf.phx.1h.cols,hf.llm.6h.cols,hf.phx.6h.cols,hf.llm.24h.cols,hf.phx.24h.cols)

cho.llm.1h.cols <- (44:46) + data.column.offset
cho.phx.1h.cols <- (47:50) + data.column.offset
cho.llm.6h.cols <- (51:53) + data.column.offset
cho.phx.6h.cols <- (54:57) + data.column.offset
cho.llm.24h.cols <- (58:59) + data.column.offset
cho.phx.24h.cols <- (60:63) + data.column.offset

cho.llm.cols <- c(cho.llm.1h.cols,cho.llm.6h.cols,cho.llm.24h.cols)
cho.phx.cols <- c(cho.phx.1h.cols,cho.phx.6h.cols,cho.phx.24h.cols)
cho.cols <- c(cho.llm.1h.cols,cho.phx.1h.cols,cho.llm.6h.cols,cho.phx.6h.cols,cho.llm.24h.cols,cho.phx.24h.cols)

data.cols <- c(cho.cols,hf.cols,etoh.cols)

median.data <- cbind(
					apply(all.data[,etoh.llm.1h.cols],1,median),
					apply(all.data[,etoh.phx.1h.cols],1,median),
					apply(all.data[,etoh.llm.6h.cols],1,median),
					apply(all.data[,etoh.phx.6h.cols],1,median),
					apply(all.data[,etoh.llm.24h.cols],1,median),
					apply(all.data[,etoh.phx.24h.cols],1,median),
					apply(all.data[,etoh.llm.1h.cols],1,median),
					apply(all.data[,hf.phx.1h.cols],1,median),
					apply(all.data[,hf.llm.6h.cols],1,median),
					apply(all.data[,hf.phx.6h.cols],1,median),
					apply(all.data[,hf.llm.24h.cols],1,median),
					apply(all.data[,hf.phx.24h.cols],1,median),
					apply(all.data[,cho.llm.1h.cols],1,median),
					apply(all.data[,cho.phx.1h.cols],1,median),
					apply(all.data[,cho.llm.6h.cols],1,median),
					apply(all.data[,cho.phx.6h.cols],1,median),
					apply(all.data[,cho.llm.24h.cols],1,median),
					apply(all.data[,cho.phx.24h.cols],1,median)
)

mean.data <- cbind(
					apply(all.data[,etoh.llm.1h.cols],1,mean),
					apply(all.data[,etoh.phx.1h.cols],1,mean),
					apply(all.data[,etoh.llm.6h.cols],1,mean),
					apply(all.data[,etoh.phx.6h.cols],1,mean),
					apply(all.data[,etoh.llm.24h.cols],1,mean),
					apply(all.data[,etoh.phx.24h.cols],1,mean),
					apply(all.data[,hf.llm.1h.cols],1,mean),
					apply(all.data[,hf.phx.1h.cols],1,mean),
					apply(all.data[,hf.llm.6h.cols],1,mean),
					apply(all.data[,hf.phx.6h.cols],1,mean),
					apply(all.data[,hf.llm.24h.cols],1,mean),
					apply(all.data[,hf.phx.24h.cols],1,mean),
					apply(all.data[,cho.llm.1h.cols],1,mean),
					apply(all.data[,cho.phx.1h.cols],1,mean),
					apply(all.data[,cho.llm.6h.cols],1,mean),
					apply(all.data[,cho.phx.6h.cols],1,mean),
					apply(all.data[,cho.llm.24h.cols],1,mean),
					apply(all.data[,cho.phx.24h.cols],1,mean)
)

mean.pattern <- t(apply(mean.data,1,function(c){(abs(c)>log2(fc.thres))*sign(c)}))
median.pattern <- t(apply(median.data,1,function(c){(abs(c)>log2(fc.thres))*sign(c)}))


median.data.phx <- cbind(
					apply(all.data[,etoh.phx.1h.cols],1,median),
					apply(all.data[,etoh.phx.6h.cols],1,median),
					apply(all.data[,etoh.phx.24h.cols],1,median),
					apply(all.data[,hf.phx.1h.cols],1,median),
					apply(all.data[,hf.phx.6h.cols],1,median),
					apply(all.data[,hf.phx.24h.cols],1,median),
					apply(all.data[,cho.phx.1h.cols],1,median),
					apply(all.data[,cho.phx.6h.cols],1,median),
					apply(all.data[,cho.phx.24h.cols],1,median)
				)

mean.data.phx <- cbind(
					apply(all.data[,etoh.phx.1h.cols],1,mean),
					apply(all.data[,etoh.phx.6h.cols],1,mean),
					apply(all.data[,etoh.phx.24h.cols],1,mean),
					apply(all.data[,hf.phx.1h.cols],1,mean),
					apply(all.data[,hf.phx.6h.cols],1,mean),
					apply(all.data[,hf.phx.24h.cols],1,mean),
					apply(all.data[,cho.phx.1h.cols],1,mean),
					apply(all.data[,cho.phx.6h.cols],1,mean),
					apply(all.data[,cho.phx.24h.cols],1,mean)
				)


mean.pattern.phx <- t(apply(mean.data.phx,1,function(c){(abs(c)>log2(fc.thres))*sign(c)}))
median.pattern.phx <- t(apply(median.data.phx,1,function(c){(abs(c)>log2(fc.thres))*sign(c)}))


patterns <- permutations(3,3,-1:1,repeats=T)
colnames(patterns) <- c("1h", "6h", "24h")
templates <- apply(patterns,1,paste,collapse="")

mean.pattern.decimal.etoh <- apply(mean.pattern.phx[,1:3]+1, 1, function(x){sum(x*3^(rev(seq_along(x))-1))})
mean.pattern.decimal.hf   <- apply(mean.pattern.phx[,4:6]+1, 1, function(x){sum(x*3^(rev(seq_along(x))-1))})
mean.pattern.decimal.cho  <- apply(mean.pattern.phx[,7:9]+1, 1, function(x){sum(x*3^(rev(seq_along(x))-1))})

# The three columns with counts correspond to the three diet groups
patterns.counts <- cbind(patterns, matrix(0,nrow=dim(patterns)[1],ncol=3))
colnames(patterns.counts) <- c(colnames(patterns),"cho.count","hf.count","etoh.count")
patterns.counts[1+as.numeric(names(table(mean.pattern.decimal.cho))),"cho.count"] <- table(mean.pattern.decimal.cho)
patterns.counts[1+as.numeric(names(table(mean.pattern.decimal.hf))),"hf.count"] <- table(mean.pattern.decimal.hf)
patterns.counts[1+as.numeric(names(table(mean.pattern.decimal.etoh))),"etoh.count"] <- table(mean.pattern.decimal.etoh)

patterns.counts.save.file <- "patterns.counts.txt"
#write.table(patterns.counts,file=patterns.counts.save.file,sep="\t",quote=F,row.names=F,col.names=T)

temp.data <- NULL
for (i in 1:27) {
	cho.matches <- (mean.pattern.decimal.cho == (i-1))
	hf.matches <- (mean.pattern.decimal.hf == (i-1))
	etoh.matches <- (mean.pattern.decimal.etoh == (i-1))
	temp.data <- rbind(temp.data)
}




mean.pattern.cho.hf.count <- matrix(NA,nrow=27,ncol=27)
rownames(mean.pattern.cho.hf.count) <- paste("cho",templates,sep='_')
colnames(mean.pattern.cho.hf.count) <- paste("hf",templates,sep='_')

for (i in 1:27) {
	cho.matches <- (mean.pattern.decimal.cho == (i-1))
	for (j in 1:27) {
		hf.matches <- (mean.pattern.decimal.hf == (j-1))
		mean.pattern.cho.hf.count[i,j] <- sum(cho.matches & hf.matches)
	}
}

cho.hf.save.file <- paste("mean.pattern.cho.hf.count-fc",fc.thres,".txt",sep="")
#write.table(mean.pattern.cho.hf.count,file=cho.hf.save.file,sep="\t",col.names=NA,quote=F)

mean.pattern.cho.etoh.count <- matrix(NA,nrow=27,ncol=27)
rownames(mean.pattern.cho.etoh.count) <- paste("cho",templates,sep='_')
colnames(mean.pattern.cho.etoh.count) <- paste("etoh",templates,sep='_')

for (i in 1:27) {
	cho.matches <- (mean.pattern.decimal.cho == (i-1))
	for (j in 1:27) {
		etoh.matches <- (mean.pattern.decimal.etoh == (j-1))
		mean.pattern.cho.etoh.count[i,j] <- sum(cho.matches & etoh.matches)
	}
}

cho.etoh.save.file <- paste("mean.pattern.cho.etoh.count-fc",fc.thres,".txt",sep="")
#write.table(mean.pattern.cho.etoh.count,file=cho.etoh.save.file,sep="\t",col.names=NA,quote=F)

mean.pattern.hf.etoh.count <- matrix(NA,nrow=27,ncol=27)
rownames(mean.pattern.hf.etoh.count) <- paste("hf",templates,sep='_')
colnames(mean.pattern.hf.etoh.count) <- paste("etoh",templates,sep='_')

for (i in 1:27) {
	hf.matches <- (mean.pattern.decimal.hf == (i-1))
	for (j in 1:27) {
		etoh.matches <- (mean.pattern.decimal.etoh == (j-1))
		mean.pattern.hf.etoh.count[i,j] <- sum(hf.matches & etoh.matches)
	}
}

hf.etoh.save.file <- paste("mean.pattern.hf.etoh.count-fc",fc.thres,".txt",sep="")
#write.table(mean.pattern.hf.etoh.count,file=hf.etoh.save.file,sep="\t",col.names=NA,quote=F)

mean.pattern.cho.hf.etoh.count <- array(NA,dim=c(27,27,27))
for (i in 1:27) {
	cho.matches <- (mean.pattern.decimal.cho == (i-1))
	for (j in 1:27) {
		hf.matches <- (mean.pattern.decimal.hf == (j-1))
		for (k in 1:27) {
			etoh.matches <- (mean.pattern.decimal.etoh == (k-1))
			mean.pattern.cho.hf.etoh.count[i,j,k] <- sum(cho.matches & etoh.matches & hf.matches)
		}
	}		
}

count.thres <- 1
pattern.genes.save.dir <- "/Users/raj/Dropbox/research/liver-regeneration/timeseries/genelists-cho-etoh-all"
for (i in 1:27) {
	cho.matches <- (mean.pattern.decimal.cho == (i-1))
	for (j in 1:27) {
		etoh.matches <- (mean.pattern.decimal.etoh == (j-1))
		pattern.genes <- rownames(all.data)[cho.matches & etoh.matches]
		if (length(pattern.genes) >= count.thres) {
			temp.prefix <- paste(pattern.genes.save.dir,"/",sprintf("%04d",length(pattern.genes)),sep="")
			pattern.genes.save.file <- paste(temp.prefix,"-c",templates[i],"-e",templates[j],".txt",sep="")
			write.table(pattern.genes,file=pattern.genes.save.file,quote=F,sep="\t",row.names=F,col.names=F)
		}
	}
}

count.thres <- 1
pattern.genes.save.dir <- "/Users/raj/Dropbox/research/liver-regeneration/timeseries/genelists-cho-hf-all"
for (i in 1:27) {
	cho.matches <- (mean.pattern.decimal.cho == (i-1))
	for (j in 1:27) {
		hf.matches <- (mean.pattern.decimal.hf == (j-1))
		pattern.genes <- rownames(all.data)[cho.matches & hf.matches]
		if (length(pattern.genes) >= count.thres) {
			temp.prefix <- paste(pattern.genes.save.dir,"/",sprintf("%04d",length(pattern.genes)),sep="")
			pattern.genes.save.file <- paste(temp.prefix,"-c",templates[i],"-h",templates[j],".txt",sep="")
			write.table(pattern.genes,file=pattern.genes.save.file,quote=F,sep="\t",row.names=F,col.names=F)
		}
	}
}

count.thres <- 1
pattern.genes.save.dir <- "/Users/raj/Dropbox/research/liver-regeneration/timeseries/genelists-cho"
for (i in 1:27) {
	cho.matches <- (mean.pattern.decimal.cho == (i-1))
	pattern.genes <- rownames(all.data)[cho.matches]
	if (length(pattern.genes) >= count.thres) {
		temp.prefix <- paste(pattern.genes.save.dir,"/",sprintf("%04d",length(pattern.genes)),sep="")
		pattern.genes.save.file <- paste(temp.prefix,"-c",templates[i],".txt",sep="")
		write.table(pattern.genes,file=pattern.genes.save.file,quote=F,sep="\t",row.names=F,col.names=F)
	}
}

count.thres <- 1
pattern.genes.save.dir <- "/Users/raj/Dropbox/research/liver-regeneration/timeseries/genelists-etoh"
for (i in 1:27) {
	etoh.matches <- (mean.pattern.decimal.etoh == (i-1))
	pattern.genes <- rownames(all.data)[etoh.matches]
	if (length(pattern.genes) >= count.thres) {
		temp.prefix <- paste(pattern.genes.save.dir,"/",sprintf("%04d",length(pattern.genes)),sep="")
		pattern.genes.save.file <- paste(temp.prefix,"-e",templates[i],".txt",sep="")
		write.table(pattern.genes,file=pattern.genes.save.file,quote=F,sep="\t",row.names=F,col.names=F)
	}
}

count.thres <- 1
pattern.genes.save.dir <- "/Users/raj/Dropbox/research/liver-regeneration/timeseries/genelists-hf"
for (i in 1:27) {
	hf.matches <- (mean.pattern.decimal.hf == (i-1))
	pattern.genes <- rownames(all.data)[hf.matches]
	if (length(pattern.genes) >= count.thres) {
		temp.prefix <- paste(pattern.genes.save.dir,"/",sprintf("%04d",length(pattern.genes)),sep="")
		pattern.genes.save.file <- paste(temp.prefix,"-h",templates[i],".txt",sep="")
		write.table(pattern.genes,file=pattern.genes.save.file,quote=F,sep="\t",row.names=F,col.names=F)
	}
}


count.thres <- 20
median.of.pattern.cho <- NULL
k <- 1
for (i in 1:27) {
	cho.matches <- (mean.pattern.decimal.cho == (i-1))
	pattern.genes <- rownames(all.data)[cho.matches]
	if (length(pattern.genes) >= count.thres) {
		median.of.pattern.cho <- rbind(median.of.pattern.cho,c(templates[i],length(pattern.genes),apply(all.data[cho.matches,],2,median)))
		k <- k+1
	}
}

colnames(median.of.pattern.cho)[1:2] <- c("Pattern","Size")
#write.table(median.of.pattern.cho,file="median_of_cho_patterns.txt",sep="\t",row.names=F,quote=F)


count.thres <- 1
all.data.ordered.cho.pattern <- NULL
for (i in 1:27) {
	cho.matches <- (mean.pattern.decimal.cho == (i-1))
	pattern.genes <- rownames(all.data)[cho.matches]
	if (length(pattern.genes) >= count.thres) {
		sort.order <- order(apply(abs(mean.data.phx[cho.matches,7:9]),1,sum),decreasing=T)
		all.data.ordered.cho.pattern <- rbind(all.data.ordered.cho.pattern,all.data[cho.matches,][sort.order,])
	}
}



count.thres <- 20
median.subpattern.phx <- NULL
k <- 1
for (i in 1:27) {
	cho.matches <- (mean.pattern.decimal.cho == (i-1))
	for (j in 1:27) {
		etoh.matches <- (mean.pattern.decimal.etoh == (j-1))
		pattern.genes <- rownames(all.data)[cho.matches & etoh.matches]
		if (length(pattern.genes) >= count.thres) {
			median.subpattern.phx <- rbind(median.subpattern.phx,apply(all.data[cho.matches&etoh.matches,],2,median))
			k <- k+1
		}
	}
}



count.thres <- 100
cho.etoh.indices <- which(mean.pattern.cho.etoh.count >= count.thres,arr.ind=T)
temp <- unique(c(cho.etoh.indices))
mean.pattern.cho.etoh.count.filtered <- mean.pattern.cho.etoh.count[temp,temp]
cho.etoh.filtered.save.file <- paste("mean.pattern.cho.etoh.count.filtered-fc",fc.thres,".txt",sep="")
#write.table(mean.pattern.cho.etoh.count.filtered,file=cho.etoh.filtered.save.file,sep="\t",col.names=NA,quote=F)


blank.cols <- matrix(NA,nrow=dim(all.data)[1],ncol=1)
data.4hmap <- cbind(all.data[,1],all.data[,etoh.cols],blank.cols,blank.cols,all.data[,hf.cols],blank.cols,blank.cols,all.data[,cho.cols])
#data.4hmap.phx <- cbind(all.data[,etoh.llm.cols], blank.cols, all.data[,etoh.phx.1h.cols],blank.cols,all.data[,etoh.phx.6h.cols],blank.cols,all.data[,etoh.phx.24h.cols],blank.cols,blank.cols,all.data[,hf.llm.cols], blank.cols,all.data[,hf.phx.1h.cols],blank.cols,all.data[,hf.phx.6h.cols],blank.cols,all.data[,hf.phx.24h.cols],blank.cols,blank.cols,all.data[,cho.llm.cols], blank.cols,all.data[,cho.phx.1h.cols],blank.cols,all.data[,cho.phx.6h.cols],blank.cols,all.data[,cho.phx.24h.cols])

data.4hmap.phx <- cbind(all.data[,etoh.phx.1h.cols],blank.cols,all.data[,etoh.phx.6h.cols],blank.cols,all.data[,etoh.phx.24h.cols],blank.cols,blank.cols,all.data[,cho.phx.1h.cols],blank.cols,all.data[,cho.phx.6h.cols],blank.cols,all.data[,cho.phx.24h.cols],blank.cols,blank.cols,all.data[,hf.phx.1h.cols],blank.cols,all.data[,hf.phx.6h.cols],blank.cols,all.data[,hf.phx.24h.cols])


#temp.data <- all.data
temp.data <- all.data.ordered.cho.pattern
#blank.cols <- cbind(blank.cols,blank.cols)
data.4hmap.phx.ordered.cho <- cbind(temp.data[,etoh.phx.1h.cols],blank.cols,temp.data[,etoh.phx.6h.cols],blank.cols,temp.data[,etoh.phx.24h.cols],blank.cols,blank.cols,temp.data[,cho.phx.1h.cols],blank.cols,temp.data[,cho.phx.6h.cols],blank.cols,temp.data[,cho.phx.24h.cols],blank.cols,blank.cols,temp.data[,hf.phx.1h.cols],blank.cols,temp.data[,hf.phx.6h.cols],blank.cols,temp.data[,hf.phx.24h.cols])
rm(temp.data)
#write.table(data.4hmap.phx.ordered.cho,file="all.data.ordered.pattern.cho.txt",sep="\t",quote=F,col.names=NA)

#temp.data <- all.data
temp.data <- all.data.ordered.cho.pattern
#blank.cols <- cbind(blank.cols,blank.cols)
data.4hmap.phx.ordered.cho <- cbind(temp.data[,cho.phx.1h.cols],blank.cols,temp.data[,cho.phx.6h.cols],blank.cols,temp.data[,cho.phx.24h.cols],blank.cols,blank.cols,temp.data[,hf.phx.1h.cols],blank.cols,temp.data[,hf.phx.6h.cols],blank.cols,temp.data[,hf.phx.24h.cols],blank.cols,blank.cols,temp.data[,etoh.phx.1h.cols],blank.cols,temp.data[,etoh.phx.6h.cols],blank.cols,temp.data[,etoh.phx.24h.cols])
rm(temp.data)
#write.table(data.4hmap.phx.ordered.cho,file="all.data.ordered.pattern.cho.hf.etoh.txt",sep="\t",quote=F,col.names=NA)



pdf("pattern-hmap-ordered-cho.pdf",width=8,height=8)
	heatmap.2(as.matrix(data.4hmap.phx.ordered.cho),Rowv=F,Colv=F,dendrogram="none",scale="none", col="greenred",trace="none",labRow="",breaks=seq(-1,1,0.01),cexCol=0.5,density.info="none",keysize=1)
dev.off()



key.pattern.counts <- cbind(cho.etoh.indices-1,mean.pattern.cho.etoh.count[cho.etoh.indices])

library(gplots)
count.thres <- 20
pdf("pattern-hmap.pdf",width=8,height=8)
for (i in 1:27) {
	for (j in 1:27) {
		temp.data <- data.4hmap.phx[which( mean.pattern.decimal.cho==(i-1) & mean.pattern.decimal.etoh==(j-1) ),]
		num.genes <- nrow(temp.data)
		if (num.genes >= count.thres) {
			if (num.genes==1) {temp.data <- rbind(temp.data,temp.data)}
			temp.title <- paste("etoh: ",templates[j],"  cho: ",templates[i],"  genes: ",num.genes,sep="")
			heatmap.2(as.matrix(temp.data),Rowv=T,Colv=F,dendrogram="none",scale="none", col="greenred",trace="none",labRow="",breaks=seq(-1.5,1.5,0.01),cexCol=0.5,main=temp.title)
		}
	}	
}
dev.off()

subset.data.save.file <- paste("cho",cho.i,".etoh",etoh.j,".txt",sep="")
#write.table(temp.data,file=subset.data.save.file,col.names=NA,sep="\t",quote=F)


# Minimum Spanning Trees
plot.minspantree <- function (data.matrix,node.colors="black",node.type="p",node.cex=2,plot.title="") {
	data.dist <- Dist(t(data.matrix),method="euclidean")
	data.spantree <- spantree (data.dist)
	plot(data.spantree,pch=20,axes=F,xlab="",ylab="",type=node.type,col=node.colors,cex=node.cex,main=plot.title)
}

cho.color <- "#0000DD"
etoh.color <- "#BB0000"
hf.color <- "#008800"

sample.colors <- NULL
sample.colors[cho.cols] <- cho.color
sample.colors[etoh.cols] <- etoh.color
sample.colors[hf.cols] <- hf.color

llm.cex <- 0.75
phx.1h.cex <- 1.5
phx.6h.cex <- 2.25
phx.24h.cex <- 3

sample.cex <- NULL
sample.cex[c(cho.llm.cols,hf.llm.cols,etoh.llm.cols)] <- llm.cex
sample.cex[c(cho.phx.1h.cols,hf.phx.1h.cols,etoh.phx.1h.cols)] <- phx.1h.cex
sample.cex[c(cho.phx.6h.cols,hf.phx.6h.cols,etoh.phx.6h.cols)] <- phx.6h.cex
sample.cex[c(cho.phx.24h.cols,hf.phx.24h.cols,etoh.phx.24h.cols)] <- phx.24h.cex

plot.data <- all.data
plot.title <- "All Data"

#plot.data <- median.subpattern.phx
#plot.title <- "median per pattern"

pdf("minimumSpanTree.pdf",width=8,height=6)

sample.cols <- c(cho.cols,hf.cols,etoh.cols)
data4tree <- plot.data[,sample.cols]
plot.minspantree(data.matrix=data4tree,node.colors=sample.colors[sample.cols],node.type="p",node.cex=sample.cex[sample.cols],plot.title=plot.title)
legend("topleft",bty="n",legend=c("CHO","HF","EtOH"),text.col=c(cho.color,hf.color,etoh.color))
legend("bottomleft",bty="n",legend=c(0,1,6,24),text.col="black",pt.cex=c(llm.cex,phx.1h.cex,phx.6h.cex,phx.24h.cex), pch=20,col="black",title="Hours Post PHx")

sample.cols <- c(cho.cols)
data4tree <- plot.data[,sample.cols]
plot.minspantree(data.matrix=data4tree,node.colors=sample.colors[sample.cols],node.type="p",node.cex=sample.cex[sample.cols],plot.title=plot.title)
legend("topleft",bty="n",legend=c("CHO"),text.col=c(cho.color))
legend("topright",bty="n",legend=c(0,1,6,24),text.col="black",pt.cex=c(llm.cex,phx.1h.cex,phx.6h.cex,phx.24h.cex), pch=20,col="black",title="Hours Post PHx")

sample.cols <- c(hf.cols)
data4tree <- plot.data[,sample.cols]
plot.minspantree(data.matrix=data4tree,node.colors=sample.colors[sample.cols],node.type="p",node.cex=sample.cex[sample.cols],plot.title=plot.title)
legend("topleft",bty="n",legend=c("HF"),text.col=c(hf.color))
legend("topright",bty="n",legend=c(0,1,6,24),text.col="black",pt.cex=c(llm.cex,phx.1h.cex,phx.6h.cex,phx.24h.cex), pch=20,col="black",title="Hours Post PHx")

sample.cols <- c(etoh.cols)
data4tree <- plot.data[,sample.cols]
plot.minspantree(data.matrix=data4tree,node.colors=sample.colors[sample.cols],node.type="p",node.cex=sample.cex[sample.cols],plot.title=plot.title)
legend("topleft",bty="n",legend=c("EtOH"),text.col=c(etoh.color))
legend("bottomright",bty="n",legend=c(0,1,6,24),text.col="black",pt.cex=c(llm.cex,phx.1h.cex,phx.6h.cex,phx.24h.cex), pch=20,col="black",title="Hours Post PHx")

sample.cols <- c(cho.cols,hf.cols)
data4tree <- plot.data[,sample.cols]
plot.minspantree(data.matrix=data4tree,node.colors=sample.colors[sample.cols],node.type="p",node.cex=sample.cex[sample.cols],plot.title=plot.title)
legend("topleft",bty="n",legend=c("CHO","HF"),text.col=c(cho.color,hf.color))
legend("bottomright",bty="n",legend=c(0,1,6,24),text.col="black",pt.cex=c(llm.cex,phx.1h.cex,phx.6h.cex,phx.24h.cex), pch=20,col="black",title="Hours Post PHx")

sample.cols <- c(cho.cols,etoh.cols)
data4tree <- plot.data[,sample.cols]
plot.minspantree(data.matrix=data4tree,node.colors=sample.colors[sample.cols],node.type="p",node.cex=sample.cex[sample.cols],plot.title=plot.title)
legend("topleft",bty="n",legend=c("CHO","EtOH"),text.col=c(cho.color,etoh.color))
legend("topright",bty="n",legend=c(0,1,6,24),text.col="black",pt.cex=c(llm.cex,phx.1h.cex,phx.6h.cex,phx.24h.cex), pch=20,col="black",title="Hours Post PHx")

sample.cols <- c(hf.cols,etoh.cols)
data4tree <- plot.data[,sample.cols]
plot.minspantree(data.matrix=data4tree,node.colors=sample.colors[sample.cols],node.type="p",node.cex=sample.cex[sample.cols],plot.title=plot.title)
legend("topleft",bty="n",legend=c("HF","EtOH"),text.col=c(hf.color,etoh.color))
legend("topright",bty="n",legend=c(0,1,6,24),text.col="black",pt.cex=c(llm.cex,phx.1h.cex,phx.6h.cex,phx.24h.cex), pch=20,col="black",title="Hours Post PHx")

sample.cols <- c(cho.phx.cols,hf.phx.cols,etoh.phx.cols)
data4tree <- plot.data[,sample.cols]
plot.minspantree(data.matrix=data4tree,node.colors=sample.colors[sample.cols],node.type="p",node.cex=sample.cex[sample.cols],plot.title=plot.title)
legend("topleft",bty="n",legend=c("CHO","HF","EtOH"),text.col=c(cho.color,hf.color,etoh.color))
legend("topright",bty="n",legend=c(1,6,24),text.col="black",pt.cex=c(phx.1h.cex,phx.6h.cex,phx.24h.cex), pch=20,col="black",title="Hours Post PHx")

sample.cols <- c(cho.llm.cols,cho.phx.cols,hf.phx.cols,etoh.phx.cols)
data4tree <- plot.data[,sample.cols]
plot.minspantree(data.matrix=data4tree,node.colors=sample.colors[sample.cols],node.type="p",node.cex=sample.cex[sample.cols],plot.title=plot.title)
legend("topleft",bty="n",legend=c("CHO","HF","EtOH"),text.col=c(cho.color,hf.color,etoh.color))
legend("topright",bty="n",legend=c(0,1,6,24),text.col="black",pt.cex=c(llm.cex,phx.1h.cex,phx.6h.cex,phx.24h.cex), pch=20,col="black",title="Hours Post PHx")

sample.cols <- c(cho.llm.cols,hf.llm.cols,cho.phx.cols,hf.phx.cols,etoh.phx.1h.cols,etoh.phx.6h.cols)
data4tree <- plot.data[,sample.cols]
plot.minspantree(data.matrix=data4tree,node.colors=sample.colors[sample.cols],node.type="p",node.cex=sample.cex[sample.cols],plot.title=plot.title)
legend("topleft",bty="n",legend=c("CHO","HF","EtOH"),text.col=c(cho.color,hf.color,etoh.color))
legend("bottomleft",bty="n",legend=c(0,1,6,24),text.col="black",pt.cex=c(llm.cex,phx.1h.cex,phx.6h.cex,phx.24h.cex), pch=20,col="black",title="Hours Post PHx")

sample.cols <- c(cho.llm.cols,cho.phx.1h.cols,cho.phx.6h.cols,hf.llm.cols,hf.phx.1h.cols,hf.phx.6h.cols,etoh.llm.cols,etoh.phx.1h.cols,etoh.phx.6h.cols)
data4tree <- plot.data[,sample.cols]
plot.minspantree(data.matrix=data4tree,node.colors=sample.colors[sample.cols],node.type="p",node.cex=sample.cex[sample.cols],plot.title=plot.title)
legend("topleft",bty="n",legend=c("CHO","HF","EtOH"),text.col=c(cho.color,hf.color,etoh.color))
legend("topright",bty="n",legend=c(0,1,6),text.col="black",pt.cex=c(llm.cex,phx.1h.cex,phx.6h.cex), pch=20,col="black",title="Hours Post PHx")

dev.off()

