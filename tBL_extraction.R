## extract terminal branch lengths of a tree

setwd("/local/workdir/sps257/dspr_piRNA/TEfam_CN_trees/terminalBL_extraction/")

trees <- list.files(path = ".", pattern = "\\.nwk$", full.names = FALSE)
for (i in 1:length(trees)) assign(trees[i], read.tree(file = trees[i]))
library(dplyr)

tBL_dfs <- lapply(trees, function(a) {
  z <- read.tree(a)
  z$edge.length <- round(z$edge.length,4)
  n <- length(z$tip.label)
  z1 <- setNames(z$edge.length[sapply(1:n,function(x,y) which(y==x),y=z$edge[,2])],z$tip.label)
  z2 <- as.data.frame(z1)
  return(z2)
})

names(tBL_dfs) <- trees

for (i in 1:length(tBL_dfs)) {
  write.table(tBL_dfs[i], file=paste0(names(tBL_dfs)[i],".txt"), row.names = T, col.names = F, sep = "\t", quote = F)
}
