require("ape")

iFile <- "~/Data/Seattle/SeattleB.fas"
seqs <- read.dna(iFile, format = "fasta", as.character = T)

#Collect Frequencies of ambiguity
fs <- sapply(1:nrow(seqs), function(i) {
  iSeq <- seqs[i,]
  lGaps <- length(iSeq[iSeq%in%"-"])
  lAmbi <- length(iSeq[iSeq%in%c("y","r","w","s","k","m","d","v","h","b","n")])
  lAmbi/(length(iSeq)-lGaps)
})

#Collect Lengths of sequences
lngs <- sapply(1:nrow(seqs), function(i) {
  iSeq <- seqs[i,]
  lGaps <- length(iSeq[iSeq%in%"-"])
  (length(iSeq) - lGaps)/length(iSeq)
})

#Write only those sequences with ambiguity below 1.5% and sequence length above 85%
write.dna(seqs[-union(which(fs>0.05),which(lngs<0.85) ),],  gsub(".fas$", "_PRO.fas", iFile), "fasta")
