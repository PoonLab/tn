require("ape")

#Take input      
iFile <- "~/Data/Tennessee/TennesseeB_PRO.fasta"
seqs <- read.FASTA(iFile, type="DNA")

metD <- read.delim("~/Data/Tennessee/sourceData/TnMetD/tnMD.csv")


idsTOT <- unname(sapply(labels(seqs), function(x) strsplit(x, '_')[[1]])[1,])
times <- sapply(idsTOT, function(id) {metD$YEAR_OF_HIV_DX[which(metD$CFAR_PID%in%id)]})
times <- times[!is.na(times)]
ids <- names(times)
times <- unname(times)

seqs <- seqs[which(idsTOT%in%ids)]

names(seqs) <- sapply(1:length(ids), function(i){paste0(ids[i],"_",times[i])})

seqs <- seqs[which(times<2012)]
times <- times[which(times<2012)]

write.FASTA(seqs, "~/Data/Tennessee/TennesseeB_PRO_Diag.fasta")

seqs <- seqs[which(times<2011)]
times <- times[which(times<2011)]

write.FASTA(seqs, "~/Data/Tennessee/TennesseeB_PRO_Diag_Filt.fasta")
