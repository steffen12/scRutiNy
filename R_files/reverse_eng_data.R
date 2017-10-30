#source("/home/scornwell/Internship/R_files/data_test_methods.R")

find_cell_size_scalars <- function(sce) {

  numCells = dim(counts(sce))[2]
  ave_counts <- rowSums(counts(sce))/numCells
  cellSizeScales <- vector(mode="double", length=numCells)
  for(cell in 1:numCells) {
    cellScales = counts(sce[, cell]) / ave_counts
    #print(cellScales)
    cellScales = cellScales[which(!is.nan(cellScales) & cellScales != 0)]
    cellSizeScale = median(cellScales)
    print(cellSizeScale)
    cellSizeScales[cell] = cellSizeScale
  }
  hist(cellSizeScales)
}

rowVar <- function(x) {
  numNotNA = vector(mode="double", length=dim(x)[1])
  for(i in 1:dim(x)[1]) {
    numNotNA[i] = length(which(!is.na(x[i,])))
  }
  return(rowSums((x - rowMeans(x, na.rm = TRUE))^2, na.rm = TRUE)/numNotNA)
}

reverse_engineer_gene_means <- function(sce) {
  library("emg")
  
  sce = computeSumFactors(sce)
  counts(sce) = t(t(counts(sce)) / pData(sce)$size_factor)
  
  geneMeans = rowMeans(counts(sce))
  geneMeans = geneMeans[which(geneMeans != 0)]
  geneMeans = log(geneMeans)
  
  posModel = tryCatch({
    emg.mle(geneMeans)
  }, error = function(err) {
    NULL
  })
  
  negModel = tryCatch({
    emg.mle(-1*geneMeans)
  }, error = function(err) {
    NULL
  })
  
  if(is.null(posModel)) {
    paramEstimates = coef(negModel)
  } else if(is.null(negModel)) {
    paramEstimates = coef(posModel)
  }
  
  mu = paramEstimates['mu']
  sigma = paramEstimates['sigma']
  lambda = paramEstimates['lambda']
  lambdaSign = -1

  y <- remg(length(geneMeans), mu = mu, sigma = sigma, lambda = lambda)
  #y <- log(remg(100, mu = 2, sigma = 1, lambda = 0.2))
  
  hist(geneMeans, main="Gene Means", breaks=100)
  hist(-y, main="EMG Gene Means", breaks=100)
}

calculate_conditional_var <- function(x, x_mean) {
  var = sum((x - x_mean)^2) / length(x)
  return(var)
}

reverse_engineer_input_parameters <- function(sce) {
  sce = calculateQCMetrics(sce)
  sce = computeSumFactors(sce)
  
  #Determine Cell Radius Disp
  sizeFactors = pData(sce)$size_factor
  sizeFactors = sizeFactors^(1/3)
  
  hist(sizeFactors)
  cellRadiusDisp = sd(sizeFactors)
  counts(sce) = t(t(counts(sce)) / pData(sce)$size_factor)
  
  #Estimate Gene Scale
  ave_counts <- rowMeans(counts(sce))
  log_ave_counts <- log(ave_counts)
  geneScale = estimateDensityPeak(log_ave_counts)
  geneScale = 2/3 * geneScale
  
  #Estimate Exponential Sign - maybe base it on median or double check to make sure it isn't wrong?
  geneMeanTotal = mean(log_ave_counts)
  geneMeanTotal = 2*geneMeanTotal
  if(geneMeanTotal <= geneScale) { #Mean below density peak
    exponentialSign = -1
  } else {
    exponentialSign = 1
  }
  
  #Estimate Exponent Scalar
  if(exponentialSign == 1) {
    conditionalSection = log_ave_counts[which(log_ave_counts < geneScale)]
  } else {
    conditionalSection = log_ave_counts[which(log_ave_counts > geneScale)]
  }
  conditionalVar = calculate_conditional_var(conditionalSection, geneScale)
  
  alphaVar = 0.05
  exponentScalar = sqrt(conditionalVar / alphaVar)
  
  #Estimate Exponential Lambda
  exponentialLambda = abs(geneMeanTotal - geneScale)/exponentScalar
  
  
  print(paste("Cell Radius Disp: ", cellRadiusDisp))
  print(paste("Gene Scale: ", geneScale))
  print(paste("Exponential Lambda: ", exponentialLambda))
  print(paste("Exponential Sign: ", exponentialSign))
  print(paste("Exponent Scalar: ", exponentScalar))
}

reverse_engineer_data <- function(sce) {
  sce <- sce[which(rowSums(counts(sce)) > 0)]
  sce = computeSumFactors(sce)
  counts(sce) = t(t(counts(sce)) / pData(sce)$size_factor)
  originalGeneMeans = rowMeans(counts(sce))
  #for (gene in 1:dim(counts(sce))[1]) {
  #  counts(sce)[gene, which(counts(sce)[gene, ] == 0)] = originalGeneMeans[gene]
  #}
  originalDropoutByGene = vector(mode="double", length=dim(counts(sce))[1])
  for(i in 1:dim(counts(sce))[1]) {
    originalDropoutByGene[i] = length(which(counts(sce)[i,] == 0)) / dim(counts(sce))[2]
  }
  plot(originalGeneMeans, originalDropoutByGene, pch = 20, main="Original Gene Dropout vs Mean", xlab="Mean Expression", ylab="Dropout by Gene", xlim=c(0,25))
  curve(exp(-7/10*x), add = TRUE, col = "red")
  curve(1/(2*x + 1), add = TRUE, col = "red")
  hist(log(originalGeneMeans), breaks=50, main="ln Gene Means", col="grey80", xlab=expression("ln(average count)"))
  numCells = dim(counts(sce))[2]
  #for(cell in 1:numCells) {
  #  counts(sce)[, cell] = counts(sce)[, cell] / pData(sce)$size_factor[cell]
  #}
  #counts(sce)[which(counts(sce) == 0)] = NA
  
  zeroMeans = rowMeans(counts(sce), na.rm = TRUE)
  nonzeroMeans = vector(mode="double", length=dim(counts(sce))[1])
  for (row in 1:dim(counts(sce))[1]) {
    nonzeroMeans[row] = sum(counts(sce)[row,]) / nnzero(counts(sce)[row,])
  }
  outlierGenes = isOutlier(zeroMeans, nmads=3, type="higher", log=TRUE)
  zeroMeans = zeroMeans[!outlierGenes]
  nonzeroMeans = nonzeroMeans[!outlierGenes]
  plot(zeroMeans, nonzeroMeans, xlim=c(0,3), ylim=c(0,3))
  
  counts(sce) = log(counts(sce))
  counts(sce)[which(is.infinite(counts(sce)))] = NA
  sce = calculateQCMetrics(sce)
  sceRowMeans = rowMeans(counts(sce), na.rm = TRUE)
  sceRowVar = rowVar(counts(sce))

  hist(sceRowMeans, xlab="Mean Expression", main="RevEng Mean Gene Expression", breaks=20, col="grey80", ylab="Number of genes")
  hist(sceRowVar, xlab="Variance of Expression", main="RevEng Variance of Gene Expression", breaks=20, col="grey80", ylab="Number of genes")
  print(paste("Average Row Mean:", mean(sceRowMeans)))
  print(paste("Average Row Var:", mean(sceRowVar)))
  dropoutByGene = vector(mode="double", length=dim(counts(sce))[1])
  dropoutByCell = vector(mode="double", length=dim(counts(sce))[2])
  for(i in 1:dim(counts(sce))[1]) {
    dropoutByGene[i] = length(which(is.na(counts(sce)[i,]))) / dim(counts(sce))[2]
  }
  for(j in 1:dim(counts(sce))[2]) {
    dropoutByCell[j] = length(which(is.na(counts(sce)[,j]))) / dim(counts(sce))[1]
  }
  hist(dropoutByGene, xlab="Dropout %", main="RevEng Dropout by Gene", breaks=20, col="grey80", ylab="Number of genes")
  hist(dropoutByCell, xlab="Dropout %", main="RevEng Dropout by Cell", breaks=20, col="grey80", ylab="Number of cells")
  print(paste("Average Dropout by Gene:", mean(dropoutByGene)))
  print(paste("Average Dropout by Cell:", mean(dropoutByCell)))
  geneRanges <- vector(mode="double", dim(counts(sce))[1])
  for(gene in 1:dim(counts(sce))[1]) {
    maxMin = range(counts(sce)[gene, ], na.rm = TRUE)
    geneRanges[gene] = maxMin[2] - maxMin[1]
  }
  plot(sceRowMeans, sceRowVar, pch = 20, main="RevEng Var vs Mean", xlab="Mean Expression", ylab="Variance of Expression")
  plot(sceRowMeans, dropoutByGene, pch = 20, main="RevEng Gene Dropout vs Mean", xlab="Mean Expression", ylab="Dropout by Gene")
  #curve(exp(-x), add = TRUE, col = "red")
  #hist(geneRanges)
  #for(gene in 1:dim(counts(sce))[1]) {
  #   print(gene)
  #   hist(counts(sce)[gene, ], breaks = 10)
  #   Sys.sleep(2)
  #}
  #plot(pData(sce_sim)$simCellSizeFactors, pData(sce_sim)$size_factor)
  #lm(pData(sce_sim)$simCellSizeFactors ~ pData(sce_sim)$size_factor)
}

test_row_transformations <- function(sce) {
  expScalar = 3
  #rowScalars = 2* runif(dim(counts(sce))[1], 0, 1)^3
  rowScalars = exp(runif(dim(counts(sce))[1], 0, 1)) - 1
  counts(sce) = counts(sce) * rowScalars
  hist(log(rowMeans(counts(sce)))/expScalar, breaks=50)
}

convertDataToSForm <- function(sce) {
  num_genes = dim(counts(sce))[1]
  num_cells = dim(counts(sce))[2]
  
  S = counts(sce) #Just to initialize matrix with right size
  
  sce = computeSumFactors(sce)
  sce_counts = t(t(counts(sce)) / pData(sce)$size_factor)
  
  for(gene in 1:num_genes) {
    rowRanks = rank(sce_counts[gene, ])
    S[gene, ] = 2 * (rowRanks / (num_cells + 1)) - 1
  }
  
  sce_counts = log(sce_counts)
  
  sce_row_means = rowMeans(sce_counts, na.rm = TRUE)
  sce_row_var = rowVar(sce_counts)
  
  sce_counts = sce_counts - rowMeans(sce_counts)
  sce_counts = sce_counts / sce_row_var
  
  save(S, file="S.RData")
}

compareClusters <- function(sce) {
  sce <- calculateQCMetrics(sce)
  numClusters = max(sce$finalCluster)
  clusterIndices = lapply(1:numClusters, function(cluster_num) which(sce$finalCluster == cluster_num))
  #clusterIndices = lapply(1:numClusters, function(cluster_num) sample(1:dim(counts(sce))[2], as.integer(dim(counts(sce))[2] / numClusters)))
  clusterSizes = sapply(1:numClusters, function(cluster_num) length(clusterIndices[[cluster_num]]))
  
  display_cols = as.integer(sqrt(numClusters))
  display_rows = as.integer(numClusters / display_cols)
  par(mfrow=c(display_rows, display_cols))
  
  clusterMeans <- vector(mode = "double", length = numClusters)
  
  for (i in 1:numClusters)
  {
    sce_cluster = sce[, clusterIndices[[i]]]
    #logCounts = log(counts(sce_cluster))
    #logCounts[which(is.infinite(logCounts))] = NA
    #ave_counts_cluster <- rowMeans(logCounts, na.rm = TRUE)
    ave_counts_cluster <- rowMeans(counts(sce_cluster))
    print(paste(i, "ln Mean Counts: ", log(mean(ave_counts_cluster))))
    clusterMeans[i] = mean(log(ave_counts_cluster))
    print(paste("Num Zero Genes for Cluster: ", length(which(ave_counts_cluster != 0))))
    ave_counts_cluster = ave_counts_cluster[which(ave_counts_cluster != 0)]
    hist(log(ave_counts_cluster), breaks=100, main="", col="grey80",
         xlab=expression("Ln average count"))
    abline(v=log10(1), col="blue", lwd=2, lty=2)
  }
  
  for (i in 1:numClusters)
  {
    sce_cluster = sce[, clusterIndices[[i]]]
    #totalCountsTransformed = 7*(sce_cluster$total_counts) - 5*mean(sce_cluster$total_counts)
    #hist(totalCountsTransformed/1e6, xlab="Library sizes (millions)", main="", 
    #     breaks=20, col="grey80", ylab="Number of cells")
    hist(sce_cluster$total_counts/1e6, xlab="Library sizes (millions)", main="", breaks=20, col="grey80", ylab="Number of cells")
    print(paste(i, "Mean Library Size (millions): ", mean(sce_cluster$total_counts/1e6)))
  }
  
  for (i in 1:numClusters)
  {
    sce_cluster = sce[, clusterIndices[[i]]]
    hist(sce_cluster$total_features, xlab="Number of expressed genes", main="", 
         breaks=20, col="grey80", ylab="Number of cells")
    print(paste(i, "Mean Number of Expressed Genes: ", mean(sce_cluster$total_features)))
  }
  
  for (i in 1:numClusters)
  {
    sce_cluster = sce[, clusterIndices[[i]]]
    var_sim = apply(counts(sce_cluster), 1, var)
    var_sim = var_sim[which(var_sim != 0)]
    hist(log(var_sim), xlab=expression("Ln Variance of Gene Expression"), main="", breaks=20, col="grey80", ylab="Number of genes")
    print(paste(i, "Mean ln Variance by Gene: ", mean(log(var_sim))))
  }
    
  # for (i in 1:numClusters)
  # {
  #   sce_cluster = sce[, clusterIndices[[i]]]
  #   plot(sce_cluster$total_features, sce_cluster$total_counts, xlab="Number of expressed genes", main="", ylab="Library sizes")
  # }
  
  # for (i in 1:numClusters)
  # {
  #   sce_cluster = sce[, clusterIndices[[i]]]
  #   pca = plotPCA(sce_cluster, pca_data_input="pdata") + fontsize
  #   plot(pca)
  # }
  
  return(list(clusterMeans = clusterMeans))
}

# compare_clusters <- function(sce) {
#   experimentName = sce$expName[1]
#   num_clusters = max(sce$finalCluster)
#   display_cols = as.integer(sqrt(num_clusters))
#   display_rows = as.integer(num_clusters / display_cols)
#   par(mfrow=c(display_rows,display_cols))
#   cluster_indices = lapply(1:num_clusters, function(cluster_num) which(sce$finalCluster == cluster_num))
#   for(i in 1:num_clusters) {
#     clusterMembers = cluster_indices[[i]]
#     print(paste("Cluster", i))
#     sce_cluster = sce[clusterMembers]
#     sce_cluster$expName = paste(experimentName, "Cluster", i);
#     plotDistribution(sce_cluster)
#     print("")
#   }
# }