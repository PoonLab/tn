args = commandArgs(trailingOnly = T)

fileNm <- paste0(args[1], "Search.txt")
start <- as.numeric(args[2])
end <- as.numeric(args[3])

fileConn <- file(fileNm)
s <- vector("character")

for (i in start:end) {
  ln <- paste0(args[1], i, ".1")
  s <- c(s, ln)
}

writeLines(s, fileConn)
close(fileConn)

