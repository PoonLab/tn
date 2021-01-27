require("ape")

#A few Basic checks to ensure the quality of an alignment
##-TO-DO: This should be optional in final package -##
alignPrePro <- function(iFile, ambiThresh=0.05, gapThresh=0.15) {
  #@param iFile: Input fasta file.
  #@param gapThresh: The threshold for gapped regions. If greater than this percentage of a given sequence are gaps, reject
  #@param ambiThresh: The threshold for sequence ambiguity. If greater than this percentage of a given sequence (without gaps) are ambiguous, reject
  
  #Take fasta file input,
  seqs <- read.dna(iFile, format = "fasta", as.character = T)
  
  #Collect Frequencies of ambiguity
  ambi <- sapply(1:nrow(seqs), function(i) {
    iSeq <- seqs[i,]
    lGaps <- length(iSeq[iSeq%in%"-"])
    lAmbi <- length(iSeq[iSeq%in%c("y","r","w","s","k","m","d","v","h","b","n",
                                   "Y","R","W","S","K","M","D","V","H","B","N")])
    lAmbi/(length(iSeq)-lGaps)
  })
  
  #Collect Lengths of sequences
  gaps <- sapply(1:nrow(seqs), function(i) {
    iSeq <- seqs[i,]
    lGaps <- length(iSeq[iSeq%in%"-"])
    lGaps/length(iSeq)
  })
  
  #Sequences that are labelled "Bad", do not meet one or more of the threshold requirements
  badSeq <- union(which(ambi>=ambiThresh),
                  which(gaps>=gapThresh))
  
  if(length(badSeq)==0) {
    print("No unnacceptable Sequences")
  } else {
    #Write only those sequences with ambiguity below 1.5% and sequence length above 85%
    write.dna(seqs[-badSeq,], colsep = "", 
              gsub(".fasta$", "_PRO.fasta", iFile), "fasta")
  }
}

#Print a graph of positions and their ambiguity as well as their gapped rate
positScan <- function(iFile, showPlot=T) {
  #@param iFile: Input fasta file.
  #@param showPlot: True false option for printing the plot of data
  #@return: A data frame of results that the plot is based
  
  #Take fasta file inpjut,
  seqs <- read.dna(iFile, format = "fasta", as.character = T)
  
  #Check highly gapped regions for trimming
  positGaps <- sapply(1:ncol(seqs), function(i) {
    iSeq <- seqs[,i]
    lGaps <- length(iSeq[iSeq%in%"-"])
    lGaps/length(iSeq)
  })
  
  #Collect Frequencies of ambiguity
  positAmbi <- sapply(1:ncol(seqs), function(i) {
    iSeq <- seqs[,i]
    lAmbi <- length(iSeq[iSeq%in%c("y","r","w","s","k","m","d","v","h","b","n",
                                   "Y","R","W","S","K","M","D","V","H","B","N")])
    lAmbi/(length(iSeq))
  })
  
  df <- data.frame(Gaps=positGaps, Amiguity=positAmbi)
  
  if(showPlot) {
    plot(df, type="n", xlim=c(1, nrow(df)), ylim = c(0,1), xlab="Position", ylab="Proportion")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "antiquewhite")
    abline(h=seq(0,1,0.2), col="white", lwd=2.5)
    abline(h=seq(0.1,0.9,0.2), col="white", lwd=1)
    lines(df$Gaps, lty=2)
    lines(df$Amiguity, lty=1)
  }
  
  return(df)
}

#Annotates headers with metadata
headerPrePro <- function(iFile, metData){
  #@param iFile: Input fasta file.
  #@param metData: A csv containing both sequence ID's and MetaData
  
  #Take fasta file input,
  seqs <- read.dna(iFile, format = "fasta", as.character = T)
  
}