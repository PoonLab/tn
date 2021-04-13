#Check Out: https://r-pkgs.org/ 

setwd("~/git/tn/packaging/")
source("analysis.R")
source("tree.setup.R")
source("pplacer.utils.R")
source("tree.clustering.R")
source("sequence.setup.R")
source("graph.setup.R")
source("graph.clustering.R")

###SEQ/TREE SETUP TESTING
seqs.full <- ape::read.FASTA("test_seq.fasta")
seq.info <- pull.headers(seqs.full,var.names = c("ID", "CollectionDate", "Subtype"),
                         var.transformations =list(as.character, as.Date, as.factor))
data.table::setnames(seq.info, "CollectionDate", "Time")
seq.info$ID <- NULL
seq.info <- annotate.new(seq.info)
t.old <- ape::read.tree("test_old_tree.nwk")
t.full <- ape::read.tree("test_full_tree.nwk")

#### GRAPH TESTING
edge.info.tn93 <- ape::dist.dna(seqs.full, pairwise.deletion = T, as.matrix = T, model = "TN93", )
edge.info.patristic <- ape::cophenetic.phylo(t.full)

g.tn93 <- create.graph(seq.info, edge.info.tn93)
g.patristic <- create.graph(seq.info, edge.info.patristic)

clusters.tn93 <- component.cluster(g.tn93, 0.007)
clusters.patristic <- component.cluster(g.patristic, 0.007)

param.list.tn93 <- lapply(seq(0,0.01,0.0001), function(x){list("g"=g.tn93, "dist.thresh"=x)})
param.list.patristic <- lapply(seq(0,0.08,0.001), function(x){list("g"=g.patristic, "dist.thresh"=x)})
cluster.range.tn93 <- multi.cluster(component.cluster, param.list.tn93, mc.cores = 4)
cluster.range.patristic <- multi.cluster(component.cluster, param.list.patristic, mc.cores = 4)

res.tn93 <- fit.analysis(cluster.range.tn93)
res.patristic <- fit.analysis(cluster.range.patristic)
  
plot.aic.diff(res.tn93)
param.list.tn93[[which.min()]]
plot.aic.diff(res.patristic)
param.list.patristic[which.min]

#### PPLACER/TREE TESTING
stats.json.ft.test <- translate.log(log.file = "test_full_log_FastTree.txt", program = "FastTree")
stats.json.rml.test <- translate.log(log.file = "test_different_log_RAxML.txt", program = "RAxML")
stats.json <- translate.log(log.file = "test_old_log_IQTREE.txt", program = "IQ-TREE")
refpkg <- taxit.create(t.old, seqs.full, stats.json)
t.grown  <- run.pplacer_guppy(refpkg)
mc.cores <- 4

t <- extend.tree(t, seq.info, mc.cores=mc.cores)
t$growth.info <- annotate.growth(t, t.grown, mc.cores=mc.cores)

clusters.step <- step.cluster(t, 0.007, 0.3)
clusters.mono <- mono.pat.cluster(t, 0.07, 0.3)

param.list <- lapply(seq(0,0.04,0.001), function(x){list("t"=t, "branch.thresh"=x)})
cluster.range <- multi.cluster(step.cluster, param.list, mc.cores = 4)

cluster.data <- cluster.range
data.table::setnames(cluster.data, "CollectionDate", "Time")

res <- fit.analysis(cluster.data)
plot.aic.diff(res)





