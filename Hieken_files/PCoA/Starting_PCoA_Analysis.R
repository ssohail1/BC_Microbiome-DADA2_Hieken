# Starting PCoA analysis
# updated Oct 26 2020
BiocManager::install("hopach")
#distmat<-distancematrix(seqtab.nochim, d = "euclid", na.rm=TRUE)
distmat<-dist(seqtab.nochim, method = "euclidean", diag = FALSE, upper = FALSE, p = 2) #distance matrix with ASV table as input
pc_plot<-pcoa(distmat)
biplot(pc_plot)
