#USAGE: Rscript partitionCSV.R 'CSVfromAlignedfas' 'folderName'
library(seqinr)
library(caret)

createInput <- function(seq, newDir) {
	# @param seq:  FIXME: What is this expected input?
	# @param newDir:  
	# @return:
	
	# TODO: write a check that dataframe <seq> has field <year>
	
	seqInput <- split(seq, f=seq$year)
	dir.create(newDir)

	for (yearSeq in seqInput) {
		yearFile <- paste(newDir, "/", yearSeq[1, 'year'], ".fas", sep="")
		write.fasta(sequences=as.list(yearSeq$sequence), names=yearSeq$id, file.out=yearFile, nbchar=60, open="w")
	}		
}

args = commandArgs(trailingOnly=TRUE)

# set stringsAsFactor to FALSE so they are handled as character objects
seq <- read.csv(args[1], stringsAsFactor=F)

set.seed(8266)

#FIXME: implement the next line in native R
index <- createDataPartition(seq$year, p=0.5, list=FALSE, times=1)

# FIXME: comment needed
train <- seq[index, ]
test <- seq[-index, ]


createInput(train, 'train')
createInput(test, 'test')

