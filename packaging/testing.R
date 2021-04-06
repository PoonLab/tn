#Check Out: https://r-pkgs.org/ 

setwd("~/git/tn/packaging/")
source("analysis.R")
source("tree.setup.R")
source("pplacer.utils.R")
source("hub.functions.R")
source("tree.clustering.R")
source("tn93.utils.R")

#### TN93 TESTING
seqs.full <- ape::read.FASTA("~/Data/NAlberta/naFullTree/refpkg/naFull.fasta")
seq.info <- pull.headers(seqs.full,var.names = c("ID", "Time", "Subtype"),
                         var.transformations =list(as.Date, as.Date, as.factor))

edge.info <- run.tn93(seqs.full, seq.info)

#### PPLACER TESTING

stats.json <- translate.log(log.file = "~/Data/NAlberta/IqTree_Complete/NorthAlbertaB_PRO.fasta.log", program = "IQ-TREE")
seqs.full <- ape::read.FASTA("~/Data/NAlberta/naFullTree/refpkg/naFull.fasta")
t <- ape::read.tree("~/Data/NAlberta/naFullTree/old.treefile")

refpkg <- taxit.create(t, seqs.full, stats.json)
ts <- run.pplacer_guppy(refpkg)

#### ANALYSIS TESTING

seqs <- ape::read.FASTA("~/Data/NAlberta/naFullTree/refpkg/naFull.fasta")
t <- ape::read.tree("~/Data/NAlberta/naFullTree/old.treefile")
t.grown <- ts
mc.cores <- 4

seq.info <- pull.headers(seqs=seqs, var.names = c("ID", "CollectionDate", "Subtype"), 
                   var.transformations =list(as.character, as.Date, as.factor), sep="_")
seq.info$ID <- NULL

t <- extend.tree(t, seq.info, mc.cores=mc.cores)
t$growth.info <- annotate.growth(t, t.grown, mc.cores=mc.cores)

clusters1 <- step.cluster(t, 0.007, 0.3)
clusters2 <- mono.pat.cluster(t, 0.07, 0.3)

param.list <- lapply(seq(0,0.04,0.001), function(x){list("t"=t, "branch.thresh"=x)})
cluster.range <- multi.cluster(step.cluster, param.list, mc.cores = 4)

cluster.data <- cluster.range
data.table::setnames(cluster.data, "CollectionDate", "Time")

res <- fit.analysis(cluster.data)
nf <- sapply(res$NullFit, function(x){x$aic})
ff <- sapply(res$FullFit, function(x){x$aic})

ff - nf
plot(ff-nf)




