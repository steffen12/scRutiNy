cluster_Sc3 <- function(sce, numClusters) {
  clusterString = paste("sc3_", numClusters, "_clusters", sep="")
  sce = sc3(sce, ks=numClusters, biology=TRUE)
  clusterIndices = vector("list", numClusters)
  p_data = pData(sce)
  for (i in 1:numClusters)
  {
    clusterIndices[[i]] = which(p_data[ , clusterString] == i)
  }
  return(clusterIndices);
}

hierarchical_cluster <- function(sce) {
  librarySizes = colSums(counts(sce))
  keepCells <- librarySizes > 0
  sce <- sce[, keepCells]
  sum(keepCells)
  ave.counts <- rowMeans(counts(sce))
  keep <- ave.counts >= 0.1 #Use 1 for read data and 0.1 for UMI data
  sum(keep)
  sce <- sce[keep,] 
  
  
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, cluster=clusters, positive=TRUE)
  
  sce <- normalize(sce)
  
  var.fit <- trendVar(sce, trend="loess", use.spikes=FALSE, span=0.2)
  var.out <- decomposeVar(sce, var.fit)
  plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
       ylab="Variance of log-expression")
  
  hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
  hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),]
  nrow(hvg.out)
  
  var.cor <- correlatePairs(sce, subset.row=rownames(hvg.out), per.gene=TRUE)
  
  sig.cor <- var.cor$FDR <= 0.05
  chosen <- var.cor$gene[sig.cor]
  
  chosen = rownames(hvg.out)
  chosen.exprs <- norm_exprs(sce)[chosen,]
  my.dist <- dist(t(chosen.exprs))
  my.tree <- hclust(my.dist, method="ward.D2")
  
  library(dynamicTreeCut)
  my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))
  
  heat.vals <- chosen.exprs - rowMeans(chosen.exprs)
  clust.col <- rainbow(max(my.clusters))
  heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.3,
            ColSideColors=clust.col[my.clusters], Colv=as.dendrogram(my.tree),
            breaks=seq(-5, 5, length.out=21))
  
  markers <- findMarkers(sce, my.clusters)
  
  top.markers <- marker.set$Gene[marker.set$Top <= 10]
  top.exprs <- norm_exprs(sce)[top.markers,,drop=FALSE]
  heat.vals <- top.exprs - rowMeans(top.exprs)
  heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.6,
            ColSideColors=clust.col[my.clusters], Colv=as.dendrogram(my.tree), dendrogram='none')
  legend("bottomleft", col=clust.col, legend=sort(unique(my.clusters)), pch=16)
}