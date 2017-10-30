cluster_CIDR <- function(sce) {
  library("cidr")
  sce_CIDR <- scDataConstructor(counts(sce))
  sce_CIDR <- determineDropoutCandidates(sce_CIDR)
  sce_CIDR <- wThreshold(sce_CIDR)
  sce_CIDR <- scDissim(sce_CIDR)
  sce_CIDR <- scPCA(sce_CIDR)
  sce_CIDR <- nPC(sce_CIDR)
  nCluster(sce_CIDR)
  sce_CIDR <- scCluster(sce_CIDR)
  plot(sce_CIDR@PC[,c(1,2)], pch=sce_CIDR@clusters, main="CIDR", xlab="PC1", ylab="PC2")
  pData(sce)$finalCluster = sce_CIDR@clusters
  plotTSNE(sce, perplexity=30, colour_by = "finalCluster", rand_seed=100) + fontsize
  return(sce)
}

cluster_seurat <- function(sce) {
  library("Seurat")
  sce_seurat <- new("seurat", raw.data = Matrix(counts(sce), sparse=T))
  sce_seurat <- Setup(sce_seurat, min.cells = 0, min.genes = 0, do.logNormalize = T, total.expr = 1e4, project = "Matri-seq")
  #pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.high = 2500)
  sce_seurat <- MeanVarPlot(sce_seurat ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
  sce_seurat <- PCA(sce_seurat, pc.genes = sce_seurat@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
  sce_seurat <- ProjectPCA(sce_seurat)
  sce_seurat <- FindClusters(sce_seurat, pc.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = T)
  #sce_seurat <- RunTSNE(sce_seurat, dims.use = 1:10, do.fast = T)
  #TSNEPlot(sce_seurat)
  sce_seurat.markers <- FindAllMarkers(sce_seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  sce_seurat.markers %>% group_by(cluster)
  
  return(sce_seurat)
}

seurat_find_proj_const <- function(sce, sce_seurat) { #Put in original SCE
  clusterAssignments = as.numeric(sce_seurat@ident)
  names(clusterAssignments) = colnames(sce)
  numClusters = max(clusterAssignments)
  pData(sce)$finalCluster = clusterAssignments
  resultList = hierarchically_cluster_cell_types(sce, numClusters)
  clusterParents = resultList$clusterParents
  clusterHeights = resultList$clusterHeights
  totalNumClusters = resultList$totalClusters
  
  projMatrix = matrix(1, nrow = totalNumClusters, ncol = nrow(counts(sce)))
  colnames(projMatrix) <- rownames(sce)
  constMatrix = matrix(0, nrow = totalNumClusters, ncol = nrow(counts(sce)))
  colnames(constMatrix) <- rownames(sce)
  
  clusterHeightsCopy = clusterHeights
  currClustHeight = min(clusterHeightsCopy)
  clusterHeightsCopy[which(clusterHeights == currClustHeight)] = Inf
  
  #Iterate hierarchically, and start from the bottom to fill in the genes unique at each level. Then iterate back
  #down to make sure that proj and const in the child cell types inherit from the parent cell types
  for(clusterID in 1:numClusters) {
    markerGeneIndexes = which(sce_seurat.markers$cluster == (clusterID - 1))
    markerGenes = sce_seurat.markers$gene[markerGeneIndexes]
    avgDiff = sce_seurat.markers$avg_diff[markerGeneIndexes]
    projMatrix[clusterID, markerGenes] = 0
    constMatrix[clusterID, markerGenes] = avgDiff
  }
  
  while(min(clusterHeightsCopy) != Inf) {
    currClustHeight = min(clusterHeightsCopy)
    currClusters = which(clusterHeights == currClustHeight)
    clusterParent = clusterParents[currClusters[1]] #Should be the same for all clusters in currClusters
    
    clusterAssignments[which(clusterAssignments == currClusters)] = clusterParent
    sce_seurat@ident = factor(clusterAssignments)
    clusterMarkers <- FindMarkers(sce_seurat, ident.1 = (clusterParent-1), min.pct = 0.25)
    
    markerGenes = rownames(clusterMarkers)
    avgDiff = clusterMarkers$avg_diff
    projMatrix[clusterParent, markerGenes] = 0
    constMatrix[clusterParent, markerGenes] = avgDiff
    
    clusterHeightsCopy[currClusters] = Inf
  }
  
  clusterHeightsCopy = clusterHeights
  clusterHeightsCopy[which(clusterHeights == 0)] = -Inf
  
  while(max(clusterHeightsCopy) != -Inf) {
    currClustHeight = max(clusterHeightsCopy)
    currClusterParents = which(clusterHeights == currClustHeight)
    for(j in 1:length(currClusterParents)) {
      clusterParent = currClusterParents[j]
      
      clusterChildren = which(clusterParents == clusterParent)
      
      for(k in 1:length(clusterChildren)) {
        clusterChild = clusterChildren[k]
        parentGeneDifferences = which(projMatrix[clusterParent] == 0 && projMatrix[clusterChild] == 0)
        projMatrix[clusterChild, parentGeneDifferences] = 0
        constMatrix[clusterChild, parentGeneDifferences] = constMatrix[clusterParent, parentGeneDifferences]
      }
      
      clusterHeightsCopy[clusterParent] = -Inf
    }
  }
  return(list(projMatrix = projMatrix, constMatrix = constMatrix))
}

reverse_engineer_clusters_seurat <- function(sce) {
  numCells = dim(counts(sce))[1]#1000
  numGenes = dim(counts(sce))[2]#2718
  
  na.rows <- which(rowSums(is.na(counts(sce))) > 0)
  if (length(na.rows) > 0) {
    sce <- sce[-na.rows, ]
  }
  sce <- sce[, colSums(counts(sce)) > 0]
  if(ncol(sce) > numCells) {
    sce <- sce[, sample(1:ncol(sce), numCells)]
  }
  sce <- sce[rowSums(counts(sce)) > 0, ]
  if(nrow(sce) > numGenes) {
    sce <- sce[sample(1:nrow(sce), numGenes), ]
  }
  
  #HANDLE EMPTY CELLS EARLIER
  n = dim(counts(sce))[1]
  cells = dim(counts(sce))[2]
  
  #I don't think it is returning this properly
  sce_seurat = cluster_seurat(sce)
  resultList = seurat_find_proj_const(sce, sce_seurat)
  projMatrix = resultList$projMatrix
  constMatrix = resultList$constMatrix
}

cluster_Sc3 <- function(sce) {
  library("SC3")
  #librarySizes = colSums(counts(sce))
  #keepCells <- librarySizes > 0
  #sce <- sce[, keepCells]
  #sum(keepCells)
  
  sce <- calculateQCMetrics(sce)
  sce <- sc3_estimate_k(sce)
  k_estimate = sce@sc3$k_estimation
  print(paste("K estimate:", k_estimate)) #What if one cluster is empty - is this possible? Would crease error with reading string
  sce <- sc3(sce, ks = k_estimate, biology=TRUE, k_estimator=FALSE)
  sc3_cluster_string = paste("sc3", k_estimate, "clusters", sep="_")
  tsne = plotTSNE(sce, perplexity=30, colour_by=sc3_cluster_string, rand_seed=100) + fontsize #Make sure to rename to cluster number!
  plot(tsne)
  pData(sce)$finalCluster = pData(sce)[[sc3_cluster_string]]
  return(sce)
  #return(pData(sce)[[sc3_cluster_string]])
}

find_sc3_proj_const <- function(sce, numClusters, keepGenes, origNumGenes) {
  #markerAUROCThresh = 0.50
  sc3_marker_auroc_string = paste("sc3", numClusters, "markers_auroc", sep="_")
  #markerAUROCThresh = sort(fData(sce)[[sc3_marker_auroc_string]], decreasing = TRUE)[numConstGenes]
  
  markerAUROCThresh = 0.70
  propMarkerKept = 1 #By Cluster?
  keepGeneIndexes = which(keepGenes)
  
  sc3_marker_auroc_string = paste("sc3", numClusters, "markers_auroc", sep="_")
  numMarkerGenes = sum(fData(sce)[[sc3_marker_auroc_string]] >= markerAUROCThresh, na.rm = TRUE)
  numFinalMarkerGenes = as.integer(propMarkerKept*numMarkerGenes)
  
  finalMarkerAUROCThresh = sort(fData(sce)[[sc3_marker_auroc_string]], decreasing = TRUE)[numFinalMarkerGenes]
  markerGeneIndexes = which(fData(sce)[[sc3_marker_auroc_string]] >= finalMarkerAUROCThresh)
  print(paste("Number of marker genes: ", length(markerGeneIndexes)))
  
  sc3_marker_cluster_string = paste("sc3", numClusters, "markers_clusts", sep="_")
  markerGeneClusters = fData(sce)[[sc3_marker_cluster_string]][markerGeneIndexes]
  
  projMatrix = matrix(1, nrow = numClusters, ncol = origNumGenes)
  constMatrix = matrix(0, nrow = numClusters, ncol = origNumGenes)
  
  for(i in 1:numClusters) {
    clusterMarkers = which(markerGeneClusters == i)
    if(length(clusterMarkers) > 0) {
      clusterMarkerIndexes = markerGeneIndexes[clusterMarkers]
      projMatrix[i, keepGeneIndexes[clusterMarkerIndexes]] = 0
      for(k in 1:length(clusterMarkerIndexes)) {
        clusterMarkerIndex = clusterMarkerIndexes[k]
        markerClusterMeanExpression = mean(exprs(sce)[clusterMarkerIndex, which(pData(sce)$finalCluster == i)])
        markerTotalMeanExpression = fData(sce)$mean_exprs[clusterMarkerIndex]
        if(markerClusterMeanExpression > markerTotalMeanExpression) {
          constMatrix[i, keepGeneIndexes[clusterMarkerIndex]] = markerClusterMeanExpression - markerTotalMeanExpression #1
        } else {
          constMatrix[i, keepGeneIndexes[clusterMarkerIndex]] = markerClusterMeanExpression - markerTotalMeanExpression #-1
        }
        print(paste(i, ":", markerClusterMeanExpression - markerTotalMeanExpression, fData(sce)[[sc3_marker_auroc_string]][clusterMarkerIndex]))
      }
    }
  }
  return(list(projMatrix = projMatrix, constMatrix = constMatrix))
}

find_diff_expr_sc3 <- function(sce) {
  #Preprocess Data]
  origNumGenes = dim(counts(sce))[1]
  ave.counts <- rowMeans(counts(sce))
  keepGenes <- ave.counts >= 0.1 #Use 1 for read data and 0.1 for UMI data
  sum(keepGenes)
  sce <- sce[keepGenes,] 
  
  sce = sc3_cluster(sce)
  clusterAssignments = pData(sce)$finalCluster
  numClusters = max(clusterAssignments) #Find a better way to get this
  
  resultList1 = find_sc3_proj_const(sce, numClusters, keepGenes, origNumGenes)
  resultList1$clusterAssignments = clusterAssignments 
  
  diffExprPValThresh = 0.01
  sc3_diff_expr_pval_string = paste("sc3", numClusters, "de_padj", sep="_")
  diffExprGenes = which(fData(sce)[[sc3_diff_expr_pval_string]] < diffExprPValThresh)
  resultList2 = hierarchically_cluster_cell_types(sce[diffExprGenes, ], numClusters)
  
  resultList = c(resultList1, resultList2)
  return(resultList)
}

reverse_engineer_clusters_sc3 <- function(sce) {
  numCells = 1000
  numGenes = 2718
  
  na.rows <- which(rowSums(is.na(counts(sce))) > 0)
  if (length(na.rows) > 0) {
    sce <- sce[-na.rows, ]
  }
  sce <- sce[, colSums(counts(sce)) > 0]
  if(ncol(sce) > numCells) {
    sce <- sce[, sample(1:ncol(sce), numCells)]
  }
  sce <- sce[rowSums(counts(sce)) > 0, ]
  if(nrow(sce) > numGenes) {
    sce <- sce[sample(1:nrow(sce), numGenes), ]
  }
  
  #HANDLE EMPTY CELLS EARLIER
  n = dim(counts(sce))[1]
  cells = dim(counts(sce))[2]
  
  #Cluster Data and Analyze Hierarchy
  resultList = find_diff_expr_sc3(sce) #Insert favorite clustering method here - cannot remove cells
  projMatrix = resultList$projMatrix
  constMatrix = resultList$constMatrix
  clusterAssignments = resultList$clusterAssignments
  clusterParents = resultList$clusterParents
  clusterHeights = resultList$clusterHeights
  totalClusters = resultList$totalClusters
  
  numClusters = max(clusterAssignments) #Find a better way to get this
  pData(sce)$finalCluster = clusterAssignments
  
  #Compare Clusters
  resultList <- compareClusters(sce, totalClusters)
  clusterMeans = resultList$clusterMeans
  
  #Sample to decrease number of genes and number of cells
  #new_n = 1000
  #new_cells = 1000
  
  #cellIndices = sample(1:cells, new_cells)
  #geneIndices = sample(1:n, new_n)
  
  #projMatrix = projMatrix[geneIndices]
  #constMatrix = constMatrix[geneIndices]
  #clusterAssignments = clusterAssignments[cellIndices]
  
  #Save Meta information
  metaInfo <- vector(mode="integer", length=3)
  metaInfo[1] = n
  metaInfo[2] = cells
  metaInfo[3] = numClusters
  
  clusterParents = clusterParents - 1 #Translate to Python indices
  save(clusterParents, file="clusterParents.RData")
  save(clusterHeights, file="clusterHeights.RData")
  save(clusterMeans, file="clusterMeans.RData")
  save(projMatrix, file="projMatrix.RData")
  save(constMatrix, file="constMatrix.RData")
  clusterAssignments = clusterAssignments - 1 #Subtract 1 for python indexing
  save(clusterAssignments, file="clusterAssignments.RData") 
  save(metaInfo, file="metaInfo.RData")
}

hierarchically_cluster_cell_types <- function(sce, numClusters) {
  clusterCentroids <- NULL
  for(i in 1:numClusters) {
    clusterCentroids <- rbind(clusterCentroids, rowMeans(exprs(sce)[, pData(sce)$finalCluster == i]))
  }
  clusterCentroids = t(clusterCentroids)
  distanceMatrix = dist(t(clusterCentroids))
  hierarchTree <- hclust(distanceMatrix, method="ward.D2")
  #plot(hierarchTree)
  
  totalClusters = dim(hierarchTree$merge)[1] + numClusters
  clusterParents <- vector(mode = "integer", length = totalClusters)
  clusterHeights <- vector(mode = "double", length = totalClusters)
  for(i in 1:dim(hierarchTree$merge)[1]) {
    for(j in 1:dim(hierarchTree$merge)[2]) {
      currClust = hierarchTree$merge[i, j]
      if(currClust > 0) {
        clusterParents[currClust + numClusters] = i + numClusters
        clusterHeights[currClust + numClusters] = hierarchTree$height[i]
      } else {
        clusterParents[abs(currClust)] = i + numClusters
        clusterHeights[abs(currClust)] = hierarchTree$height[i]
      }
    }
  }
  return(list(clusterParents = clusterParents, clusterHeights = clusterHeights, totalClusters = totalClusters))
}