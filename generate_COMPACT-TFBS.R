setwd('/Users/Albert/Dropbox/Research UD 2011/R Workspace')

peak.genes <- read.table("sample.txt",stringsAsFactors=F,header=T,fill=T,sep="\t")

samples <- colnames(peak.genes)

num.genes.uniq <- NULL
for (i in 1:dim(peak.genes)[2]) {
	num.genes.uniq <- c(num.genes.uniq,length(which(unique(peak.genes[,i])!="")))
}

num.genes <- NULL
for (i in 1:dim(peak.genes)[2]) {
	num.genes <- c(num.genes,length(which(peak.genes[,i]!="")))
}

edge.data <- NULL

for (i in 2:dim(peak.genes)[2]) {
	for (j in 1:(i-1)) {
		num.genes.common <- sum(unique(peak.genes[,j]) %in% unique(peak.genes[,i]))
		if (num.genes.common >0) {
				edge.data <- rbind(edge.data, c(samples[i],"ss",samples[j],num.genes.common))
		}
	}
}

all.genes <- NULL

for (i in 1:dim(peak.genes)[2]) {
	all.genes <- c(all.genes, unique(peak.genes[,i]))
}

uniq.genes <- unique(all.genes)
uniq.genes <- uniq.genes[uniq.genes!=""]

peak.genes.uniq <- NULL

for (i in 1:dim(peak.genes)[2]) {
	peak.genes.uniq <- cbind(peak.genes.uniq, uniq.genes %in% peak.genes[,i])
}

peak.genes.uniq.numeric <- cbind(uniq.genes,apply(peak.genes.uniq,2,as.numeric))
colnames(peak.genes.uniq.numeric) <- c("Gene",samples)

write.table(peak.genes.uniq.numeric,file="nfkb-genes-uniq.txt",row.names=F,quote=F,sep="\t")

cho.0h <- apply(peak.genes.uniq[,1:3],1,sum)
cho.1h <- apply(peak.genes.uniq[,4:6],1,sum)
cho.6h <- apply(peak.genes.uniq[,7:10],1,sum)
chow.0h <- apply(peak.genes.uniq[,11:13],1,sum)
chow.1h <- apply(peak.genes.uniq[,14:16],1,sum)
chow.6h <- apply(peak.genes.uniq[,17:19],1,sum)

pattern <- t(rbind(cho.0h, cho.1h, cho.6h, chow.0h, chow.1h, chow.6h))
cpattern <- t(rbind(cho.0h, cho.1h, cho.6h))
chpattern <- t(rbind(chow.0h, chow.1h, chow.6h))
pattern.uniq <- unique(pattern)
cpattern.uniq <- unique(cpattern)
chpattern.uniq <- unique(chpattern)

pattern.write <- cbind(uniq.genes,pattern)
cpattern.write <- cbind(uniq.genes,cpattern)
chpattern.write <- cbind(uniq.genes,chpattern)
colnames(pattern.write) <- c("Gene",colnames(pattern))
colnames(cpattern.write) <- c("Gene",colnames(cpattern))
colnames(chpattern.write) <- c("Gene",colnames(chpattern))
write.table(pattern.write, file="nfkb-pattern.txt",row.names=F,quote=F,sep="\t")
write.table(cpattern.write, file="nfkb-cpattern.txt",row.names=F,quote=F,sep="\t")
write.table(chpattern.write, file="nfkb-chpattern.txt",row.names=F,quote=F,sep="\t")

pattern.binary <- t(apply(pattern,1,function(a){round(a/c(3,3,4,3,3,3))}))
pattern.binary.uniq <- unique(pattern.binary)
pattern.binary.write <- cbind(uniq.genes, pattern.binary)
colnames(pattern.binary.write) <- c("Gene", colnames(pattern.binary))
write.table(pattern.binary.write, file="nfkb-pattern-binary.txt",row.names=F,quote=F,sep="\t")

cpattern.binary <- t(apply(cpattern,1,function(a){round(a/c(3,3,4))}))
cpattern.binary.uniq <- unique(cpattern.binary)
cpattern.binary.write <- cbind(uniq.genes, cpattern.binary)
colnames(cpattern.binary.write) <- c("Gene", colnames(cpattern.binary))
write.table(cpattern.binary.write, file="nfkb-cpattern-binary.txt",row.names=F,quote=F,sep="\t")

chpattern.binary <- t(apply(chpattern,1,function(a){round(a/c(3,3,3))}))
chpattern.binary.uniq <- unique(chpattern.binary)
chpattern.binary.write <- cbind(uniq.genes, chpattern.binary)
colnames(chpattern.binary.write) <- c("Gene", colnames(chpattern.binary))
write.table(chpattern.binary.write, file="nfkb-chpattern-binary.txt",row.names=F,quote=F,sep="\t")


pattern.decimal <- apply(pattern.binary,1,function(x){sum(x*2^(rev(seq_along(x))-1))})
pattern.decimal.write <- cbind(uniq.genes, pattern.decimal)
write.table(pattern.decimal.write, file="nfkb-pattern-decimal.txt",row.names=F,quote=F,sep="\t")

cpattern.decimal <- apply(cpattern.binary,1,function(x){sum(x*2^(rev(seq_along(x))-1))})
cpattern.decimal.write <- cbind(uniq.genes, cpattern.decimal)
write.table(cpattern.decimal.write, file="nfkb-cpattern-decimal.txt",row.names=F,quote=F,sep="\t")

chpattern.decimal <- apply(chpattern.binary,1,function(x){sum(x*2^(rev(seq_along(x))-1))})
chpattern.decimal.write <- cbind(uniq.genes, chpattern.decimal)
write.table(chpattern.decimal.write, file="nfkb-chpattern-decimal.txt",row.names=F,quote=F,sep="\t")

both_pattern.decimal.write <- cbind(uniq.genes, cpattern.decimal, chpattern.decimal)
write.table(both_pattern.decimal.write, file="nfkb-both-pattern-decimal.txt",row.names=F,quote=F,sep="\t")


pattern.binary.ordered <- pattern.binary[order(pattern.decimal,decreasing=T),]
pattern.binary.string <- apply(pattern.binary,1,paste,collapse="")


c.levels <- as.numeric(levels(as.factor(cpattern.decimal)))
e.levels <- as.numeric(levels(as.factor(chpattern.decimal)))
total.levels <- unique(c(c.levels,e.levels))

cor.ce <- matrix(NA,nrow=length(total.levels),ncol=length(total.levels))
rownames(cor.ce) <- paste("c",total.levels,sep='')
colnames(cor.ce) <- paste("ch",total.levels,sep='')
for (i in 1:length(total.levels)) {
	c.i <- (cpattern.decimal == total.levels[i])
	for (j in 1:length(total.levels)) {
		e.j <- (chpattern.decimal == total.levels[j])
		cor.ce[i,j] <- sum(c.i * e.j)
	}
}

write.table(cor.ce,file="cor.ce.txt",sep="\t",col.names=NA)

cor.ce.4circos <- cor.ce

cor.ce.4circos["c0","ch0"] <- 0

diag1 <- 0 * cor.ce
rownames(diag1) <- paste("ch",total.levels,sep='')
colnames(diag1) <- paste("ch",total.levels,sep='')
temp1 <- rbind(diag1, cor.ce.4circos)
temp2 <- rbind(t(cor.ce.4circos),diag1)
cor.ce.4circos <- cbind(temp1,temp2)

write.table(cor.ce.4circos,file="cor.ce.4circos.txt",sep="\t",col.names=NA)

# The following is a horrible way to get the genes in each interaction pattern

# Initialize pattern.interact
pattern.interact <- matrix(NA,(max(cor.ce)),length(cor.ce))
# Put ordered genes into pattern.interact
pattern.decimal.ordered <- pattern.decimal.write[order(pattern.decimal)]
pattern.interact[,1] <- (pattern.decimal.ordered[1:cor.ce[1,1]])
for(i in 2:length(pattern.interact[1,])) {
  pattern.interact[1:cor.ce[i],i] <- (pattern.decimal.ordered[(sum(cor.ce[1:(i-1)])+1):((sum(cor.ce[1:(i-1)]))+cor.ce[i])])
}

# Write interaction patterns to a file
colnames(cor.ce) <- 
colnames(pattern.interact) <- paste(1:64)
write.table(pattern.interact, file="interaction-pattern-gene-list.txt",row.names=F,quote=F,sep="\t")









