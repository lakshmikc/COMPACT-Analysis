#################
#
# Libraries to load at the start
#
#################
library(RCurl)
library(gplots)
library(gtools)
library(rgl)
library(amap)

#################
#
# Useful functions
#
#################


#
# Function to read in the data from the input file
# e.g., gene expression, ChIP-chip, etc
# One row per measure (gene, genomic location, etc)
# One column per sample
# Can have multiple annotation columns that are not quantitative data from samples
#	Specify the column and row numbers for start of the numeric data
#
read.input.data.old <- function(input.data.file,input.data.dir="",
								data.column.start=2,data.row.start=2,header.row=1) {

	input.data.file.full.path <- paste(input.data.dir,input.data.file,sep="/")

	# Sanity check for parameters
	if (header.row > data.row.start) {
		stop("Error: header.row value cannot be larger than data.row.start\n")	
	}

	# Figure out how many columns are in the file.
	skip.rows <- header.row - 1
	header.from.file <- read.table(input.data.file.full.path,skip=skip.rows,nrows=1,
									header=F,stringsAsFactors=F,quote="\"")

	num.columns <- dim(temp.from.file)[2]	
	data.columns <- data.column.start:num.columns

	data.column.classes <- rep("character",times=num.columns)
	data.column.classes[data.columns] <- "numeric"

	temp.from.file <- read.table(input.data.file.full.path,skip=data.row.start-1, 
								row.names=1,colClasses=data.column.classes,
								header=F,stringsAsFactors=F,quote="\"")

	input.data <- NULL
	input.data$samples <- header.from.file[data.columns]
	input.data$table   <- temp.from.file[,data.columns]
	input.data$row.annot <- temp.from.file[,1:data.column.start-1]
	input.data$col.annot <- temp.from.file
	
	return(input.data)
}



read.input.data <- function(input.data.file,input.data.dir=getwd(),
								data.column.start=2,data.row.start=2,header.row=1,skip.rows=0,data.columns=NULL) {

	input.data.file.full.path <- paste(input.data.dir,input.data.file,sep="/")

	# Sanity check for parameters
	if (header.row > data.row.start) {
		stop("Error: header.row value cannot be larger than data.row.start\n")	
	}

	temp.from.file <- read.table(input.data.file.full.path,sep="\t",fill=T,
									header=F,stringsAsFactors=F,skip=skip.rows,,quote="\"")

	num.rows    <- dim(temp.from.file)[1]
	num.columns <- dim(temp.from.file)[2]
	data.rows    <- data.row.start:num.rows
	if (is.null(data.columns)) {
		data.columns <- data.column.start:num.columns
	}

	input.data <- NULL
	input.data$sample.ids <- temp.from.file[header.row,data.columns]
	input.data$table   <- apply((temp.from.file[data.rows,data.columns]),2,as.numeric)
	input.data$col.annot <- temp.from.file[header.row:(data.row.start-1),c(data.column.start-1,data.columns)]
	input.data$col.annot[1,1] <- "Sample.ID"
	input.data$row.annot <- temp.from.file[c(header.row,data.rows),1:(data.column.start-1)]
	
	return(input.data)
}




#
# Function to create an integrated sample info structure from a tab delimited text file
#  or from the col.annot information in the input.data
read.sample.info <- function (sample.info.data,sample.id.column,sample.type.column,
									comparative.condition.column,pattern.condition.column) {
	sample.info <- NULL

	if (is.character(sample.info.data)) {
		temp.annot <- read.table(sample.info.data,header=F,sep="\t",stringsAsFactors=F,quote="\"")
	} else {
		temp.annot <- t(sample.info.data$col.annot)
	}
	
	# Reading the input with the header as part of the data
	# Separate the sample information and the factor names
	sample.info$table       <- temp.annot[-1,]
	sample.info$factor.names <- sample.info.data$col.annot[,1]
	
	sample.info$id.column   <- sample.id.column
	sample.info$type.column <- sample.type.column
	sample.info$comparative.condition.column <- comparative.condition.column
	sample.info$pattern.condition.column <- pattern.condition.column
	
	
	if (length(sample.info$comparative.condition.column) ==1 ) {
		sample.info$comparative.groups <- as.factor(sample.info$table[,sample.info$comparative.condition.column])
	} else {
		sample.info$comparative.groups <- as.factor(apply(sample.info$table[,sample.info$comparative.condition.column],1,paste,collapse='-'))
	}

	if (length(sample.info$pattern.condition.column) ==1 ) {
		sample.info$pattern.groups <- as.factor(sample.info$table[,sample.info$pattern.condition.column])
	} else {
		sample.info$pattern.groups <- as.factor(apply(sample.info$table[,sample.info$pattern.condition.column],1,paste,collapse='-'))
	}

	return(sample.info)
}



# Function to calculate the average data for each comparative group and pattern group
calc.average.per.group.per.gene <- function (input.data,sample.info) {
	
	comparative.groups <- levels(sample.info$comparative.groups)
	pattern.groups     <- levels(sample.info$pattern.groups)
	
	average.data <- NULL
	
	for (cg in comparative.groups) {
		for (pg in pattern.groups) {
			cg.pg.indices <- which((sample.info$comparative.groups == cg)
									& (sample.info$pattern.groups == pg))
			cg.pg.samples.from.info <- sample.info$table[cg.pg.indices,sample.info$id.column]
			
			cg.pg.samples.in.data <- input.data$sample.ids %in% cg.pg.samples.from.info
			
			cg.pg.data <- apply(input.data$table[,cg.pg.samples.in.data],1,mean)
			
			average.data$table <- cbind(average.data$table,cg.pg.data)
			average.data$annot <- cbind(average.data$annot,rbind(cg,pg))
		}
	}
	
	return(average.data)
}


# Function to calculate discretized pattern data for each comparative group, based on pattern group
calc.pattern.discrete.per.group.per.gene <- function (average.data,sample.info,threshold.value=log2(1.5)) {

	comparative.groups.in.info	<- levels(sample.info$comparative.groups)
	pattern.groups.in.info 		<- levels(sample.info$pattern.groups)
	
	pattern.discrete <- NULL
	
	for (cg in comparative.groups.in.info) {
		cg.indices.in.data <- (average.data$annot["cg",] == cg)
		pg.indices.in.data <- average.data$annot["pg",] %in% pattern.groups.in.info
		cg.pg.indices.in.data <- cg.indices.in.data & pg.indices.in.data
		
		cg.pattern <- apply(average.data$table[,cg.pg.indices.in.data],2,
								function(c){(abs(c)>threshold.value)*sign(c)})
		
		pattern.discrete$table <- cbind(pattern.discrete$table, cg.pattern)
		pattern.discrete$annot <- cbind(pattern.discrete$annot, average.data$annot[,cg.pg.indices.in.data])
	}
	
	return(pattern.discrete)
}


calc.patterns.counts.per.group <- function(pattern.data,sample.info,discrete.levels) {
	comparative.groups.in.info	<- levels(sample.info$comparative.groups)
	pattern.groups.in.info 		<- levels(sample.info$pattern.groups)
	
	pattern.count.data <- NULL
	
	for (cg in comparative.groups.in.info) {
		cg.indices.in.data <- (pattern.data$annot["cg",] == cg)
		pg.indices.in.data <- pattern.data$annot["pg",] %in% pattern.groups.in.info
		cg.pg.indices.in.data <- cg.indices.in.data & pg.indices.in.data
		
		cg.patterns <- apply(pattern.data$table[,cg.pg.indices.in.data],1,paste,sep="",collapse="_")
		
		cg.patterns.decimal <- apply(pattern.data$table[,cg.pg.indices.in.data]+1,1,
										function(x,d){sum(x*d^(rev(seq_along(x))-1))}, discrete.levels[cg])
										
		cg.pattern.counts <- table(cg.patterns)
		
		pattern.count.data$patterns[[cg]] <- cg.patterns
		pattern.count.data$patterns.decimal[[cg]] <- cg.patterns.decimal
		pattern.count.data$counts[[cg]] <- cg.pattern.counts
		
	}
	
	return(pattern.count.data)
}


calc.compact.matrix <- function(pattern.data,sample.info,discrete.levels,symmetric=FALSE) {

	patterns.counts.per.cg.pg <- calc.patterns.counts.per.group(pattern.data,sample.info,discrete.levels)
	
	compact.matrix <- NULL
	
	cgs.in.data <- names(patterns.counts.per.cg.pg$patterns)
	
	cg.combinations <- combinations(n=length(cgs.in.data),r=2)
	
	for (cg.c.i in 1:nrow(cg.combinations)) {
		
		cg1 <- cgs.in.data[cg.combinations[cg.c.i,1]]
		cg2 <- cgs.in.data[cg.combinations[cg.c.i,2]]
	
		cg1.patterns <- rownames(patterns.counts.per.cg.pg$counts[[cg1]])
		cg2.patterns <- rownames(patterns.counts.per.cg.pg$counts[[cg2]])
	
		if (symmetric) {
			cg1.patterns <- sort(unique(c(cg1.patterns,cg2.patterns)))
			cg2.patterns <- cg1.patterns	
		}
	
		compact.matrix.per.cg1.cg2 <- matrix(NA,nrow=length(cg1.patterns),ncol=length(cg2.patterns))
		rownames(compact.matrix.per.cg1.cg2) <- cg1.patterns
		colnames(compact.matrix.per.cg1.cg2) <- cg2.patterns

		for (i in 1:length(cg1.patterns)) {
			cg1.matches <- (patterns.counts.per.cg.pg$patterns[[cg1]] == cg1.patterns[i])
			for (j in 1:length(cg2.patterns)) {
				cg2.matches <- (patterns.counts.per.cg.pg$patterns[[cg2]] == cg2.patterns[j])
				compact.matrix.per.cg1.cg2[i,j] <- sum(cg1.matches & cg2.matches)
			}
		}
	
		compact.matrix[[paste(cg1,cg2,sep=":")]] <- compact.matrix.per.cg1.cg2
	}
	
	return(compact.matrix)
}


# Convert COMPACT matrix to Circos input to generate CIMPACT view
compact2circos <- function (compact.matrix, color.palette,
							circos.input.file.prefix="circos-input", pattern.counts=NULL) {
	
	circos.matrix <- NULL
	
	cg.pairs <- names(compact.matrix)
	
	for (cgp in cg.pairs) {
		cg1.name <- cgp[1]
		cg1.patterns <- rownames(compact.matrix[[1]])
		cg1.pattern.names <- paste(cg1.name,":",cg1.patterns,sep="")
		
		cg2.patterns <- colnames(compact.matrix[[1]])

		unique.patterns <- unique(c(cg1.patterns,cg2.patterns))
	
		color.panel <- color.palette[1:length(unique.patterns)]
		names(color.panel) <- unique.patterns

		if (is.null(pattern.counts)) {
			
			
	
		
		}

	
		circos.save.file.per.cg <- paste(paste(circos.input.file.prefix, cgp, sep="-"),".txt",sep="")
		write.table(circos.matrix[[cgp]], file=circos.save.file.per.cg, sep="\t",col.names=F,row.names=F)
	}
	
	
	return (circos.matrix)	
}





