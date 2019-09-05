require("ape")

iFile <- "Data/Seattle/SeattleB.fas"
seqs <- read.FASTA(iFile, type = "DNA")

#Collect Frequencies of ambiguity
fs <- sapply(1:length(seqs), function(i) {
  iSeq <- seqs[i]
  f <- base.freq(iSeq, all=TRUE)
  sum(f[c("r", "m", "w", "s", "k", "v", "h", "d", "b", "n", "?")]>0.015)
})

#Write only those sequences with ambiguity below 1.5%
write.FASTA(seqs[-which(fs>0.015)],  gsub(".fas$", "_PRO.fas", iFile)  )
