qualityControl <- function(sce) {
  ave.counts <- rowSums(counts(sce))/dim(counts(sce))[2]
  keep <- ave.counts >= 1
  sum(keep)
  
  sce <- sce[keep,]
  sce <- calculateQCMetrics(sce)
  
  libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="higher", log=TRUE)
  feature.drop <- isOutlier(sce$total_features, nmads=3, type="higher", log=TRUE)
  sce <- sce[,!(libsize.drop | feature.drop)]
  data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), Remaining=ncol(sce))
  
  return(sce);
}

estimateDensityPeak <- function(x) {
  den <- density(x)
  return(den$x[which.max(den$y)])
}

compareDistributions <- function(sce1, sce2) {
  par(mfrow=c(1,2))
  
  sce1Name = sce1$expName[1] #Take first index because the there is an array of the same name
  sce2Name = sce2$expName[1]
  sce1Title = paste(sce1Name, "Data")
  sce2Title = paste(sce2Name, "Data")
  
  sce1 = calculateQCMetrics(sce1)
  sce2 = calculateQCMetrics(sce2)
  
  #Can cause differences
  sce1 <- sce1[which(rowSums(counts(sce1)) > 0)]
  sce1 = computeSumFactors(sce1)
  counts(sce1) = t(t(counts(sce1)) / pData(sce1)$size_factor)
  
  sce2 <- sce2[which(rowSums(counts(sce2)) > 0)]
  sce2 = computeSumFactors(sce2)
  counts(sce2) = t(t(counts(sce2)) / pData(sce2)$size_factor)
  
  #Compare Ln Mean Gene Expression Levels
  ave_counts_sim <- rowSums(counts(sce1))/dim(counts(sce1))[2]
  hist(log(ave_counts_sim), breaks=100, main=sce1Title, col="grey80", xlab=expression(Ln~"average count"))
  abline(v=log(1), col="blue", lwd=2, lty=2)
  print(paste("Simulated ln(Gene Mean)", mean(log(ave_counts_sim))))
  print(paste("Simulated Density Peak", estimateDensityPeak(log(ave_counts_sim))))
  ave_counts_real <- rowSums(counts(sce2))/dim(counts(sce2))[2]
  hist(log(ave_counts_real), breaks=100, main=sce2Title, col="grey80", xlab=expression(Ln~"average count"))
  abline(v=log(1), col="blue", lwd=2, lty=2)
  print(paste("Real ln(Gene Mean)", mean(log(ave_counts_real))))
  print(paste("Real Density Peak", estimateDensityPeak(log(ave_counts_real))))
  
  #Compare Library Sizes
  hist(sce1$total_counts/1e6, xlab="Library sizes (millions)", main=sce1Title, breaks=20, col="grey80", ylab="Number of cells")
  hist(sce2$total_counts/1e6, xlab="Library sizes (millions)", main=sce2Title, breaks=20, col="grey80", ylab="Number of cells")

  #Compare Number of Expressed Genes
  hist(sce1$total_features, xlab="Number of expressed genes", main=sce1Title, breaks=20, col="grey80", ylab="Number of cells")
  hist(sce2$total_features, xlab="Number of expressed genes", main=sce2Title, breaks=20, col="grey80", ylab="Number of cells")
  
  #Compare Number of Nonzero Cells
  hist(rowSums(counts(sce1) != 0), xlab="Number of active cells", main=sce1Title, breaks=20, col="grey80", ylab="Number of genes")
  hist(rowSums(counts(sce2) != 0), xlab="Number of active cells", main=sce2Title, breaks=20, col="grey80", ylab="Number of genes")

  #Compare Gene Variances
  var_sim = apply(counts(sce1), 1, var)
  hist(log(var_sim), xlab=expression(Ln~"Variance of Gene Expression"), main=sce1Title, breaks=20, col="grey80", ylab="Number of genes")
  var_real = apply(counts(sce2), 1, var)
  hist(log(var_real), xlab=expression(Ln~"Variance of Gene Expression"), main=sce2Title, breaks=20, col="grey80", ylab="Number of genes")
  
  #Compare Mean-Variance Trend
  var.fit_sim <- trendVar(sce1, trend="loess", use.spikes=FALSE, span=0.2)
  var.out_sim <- decomposeVar(sce1, var.fit_sim)
  plot(var.out_sim$mean, var.out_sim$total, pch=16, cex=0.6, main=sce1Title, xlab="Mean log-expression", 
       ylab="Variance of log-expression")
  o <- order(var.out_sim$mean)
  lines(var.out_sim$mean[o], var.out_sim$tech[o], col="dodgerblue", lwd=2)
  
  var.fit_real <- trendVar(sce2, trend="loess", use.spikes=FALSE, span=0.2)
  var.out_real <- decomposeVar(sce2, var.fit_real)
  plot(var.out_real$mean, var.out_real$total, pch=16, cex=0.6, main=sce2Title, xlab="Mean log-expression", 
       ylab="Variance of log-expression")
  o <- order(var.out_real$mean)
  lines(var.out_real$mean[o], var.out_real$tech[o], col="dodgerblue", lwd=2)
  
  #Compare Gene Mean-Dropout Trend
  gene_dropout_sim = vector(mode="double", length=dim(counts(sce1))[1])
  for(i in 1:dim(counts(sce1))[1]) {
    gene_dropout_sim[i] = length(which(counts(sce1)[i,] == 0)) / dim(counts(sce1))[2]
  }
  plot(log(ave_counts_sim), gene_dropout_sim, pch = 20, main=sce1Title, xlab="ln(Mean Expression)", ylab="Dropout by Gene", xlim=c(-5,5))
  
  gene_dropout_real = vector(mode="double", length=dim(counts(sce2))[1])
  for(i in 1:dim(counts(sce2))[1]) {
    gene_dropout_real[i] = length(which(counts(sce2)[i,] == 0)) / dim(counts(sce2))[2]
  }
  plot(log(ave_counts_real), gene_dropout_real, pch = 20, main=sce2Title, xlab="ln(Mean Expression)", ylab="Dropout by Gene", xlim=c(-5,5))

  #Compare Size Factor Distributions
  hist(pData(sce1)$size_factor, breaks=20, main=sce1Title, xlab="Size Factors", col="grey80")
  hist(pData(sce2)$size_factor, breaks=20, main=sce2Title, xlab="Size Factors", col="grey80")
}

plotDistribution <- function(sce1) {
  sce1Name = sce1$expName[1] #Take first index because the there is an array of the same name
  sce1Title = paste(sce1Name, "Data")

  sce1 = calculateQCMetrics(sce1)
  
  #Can cause differences
  sce1 <- sce1[which(rowSums(counts(sce1)) > 0)]
  sce1 = computeSumFactors(sce1)
  counts(sce1) = t(t(counts(sce1)) / pData(sce1)$size_factor)

  #Log[10] Mean Gene Expression Levels
  ave_counts_sim <- rowSums(counts(sce1))/dim(counts(sce1))[2]
  hist(log10(ave_counts_sim), breaks=50, main=sce1Title, col="grey80", xlab=expression(Log[10]~"average count"))
  abline(v=log10(1), col="blue", lwd=2, lty=2)
  print(paste("Mean log10 Counts: ", mean(log10(ave_counts_sim))))
  
  #Library Sizes
  hist(sce1$total_counts/1e6, xlab="Library sizes (millions)", main=sce1Title, breaks=20, col="grey80", ylab="Number of cells")
  print(paste("Mean Library Size (millions): ", mean(sce1$total_counts/1e6)))
  
  #Number of Expressed Genes
  hist(sce1$total_features, xlab="Number of expressed genes", main=sce1Title, breaks=20, col="grey80", ylab="Number of cells")
  print(paste("Mean Number of Expressed Genes: ", mean(sce1$total_features)))
  
  #Number of Nonzero Cells
  hist(rowSums(counts(sce1) != 0), xlab="Number of active cells", main=sce1Title, breaks=20, col="grey80", ylab="Number of genes")
  print(paste("Mean Number of Active Cells: ", mean(rowSums(counts(sce1) != 0))))
  
  #Compare Gene Variances
  var_sim = apply(counts(sce1), 1, var)
  hist(log10(var_sim), xlab=expression(Log[10]~"Variance of Gene Expression"), main=sce1Title, breaks=20, col="grey80", ylab="Number of genes")
  print(paste("Mean log10 Variance by Gene: ", mean(log10(var_sim))))
  
  #Compare Mean-Variance Trend
  var.fit_sim <- trendVar(sce1, trend="loess", use.spikes=FALSE, span=0.2)
  var.out_sim <- decomposeVar(sce1, var.fit_sim)
  plot(var.out_sim$mean, var.out_sim$total, pch=16, cex=0.6, main=sce1Title, xlab="Mean log-expression", 
       ylab="Variance of log-expression")
  o <- order(var.out_sim$mean)
  lines(var.out_sim$mean[o], var.out_sim$tech[o], col="dodgerblue", lwd=2)
  
  #Compare Gene Mean-Dropout Trend
  gene_dropout_sim = vector(mode="double", length=dim(counts(sce1))[1])
  for(i in 1:dim(counts(sce1))[1]) {
    gene_dropout_sim[i] = length(which(counts(sce1)[i,] == 0)) / dim(counts(sce1))[2]
  }
  plot(log(ave_counts_sim), gene_dropout_sim, pch = 20, main=sce1Title, xlab="ln(Mean Expression)", ylab="Dropout by Gene", xlim=c(-5,5))
  
  #Compare Size Factor Distributions
  hist(pData(sce1)$size_factor, breaks=20, main=sce1Title, xlab="Size Factors", col="grey80")
}

# compareClusterDistributions <- function(sce) {
#   sce1Name = sce1$expName[1] #Take first index because the there is an array of the same name
#   sce1Title = paste(sce1Name, "Data")
#   
#   sce1 = calculateQCMetrics(sce1)
#   
#   #Log[10] Mean Gene Expression Levels
#   ave_counts_sim <- rowSums(counts(sce1))/dim(counts(sce1))[2]
#   hist(log10(ave_counts_sim), breaks=50, main=sce1Title, col="grey80", xlab=expression(Log[10]~"average count"))
#   abline(v=log10(1), col="blue", lwd=2, lty=2)
#   print(paste("Mean log10 Counts: ", mean(log10(ave_counts_sim))))
#   
#   #Library Sizes
#   hist(sce1$total_counts/1e6, xlab="Library sizes (millions)", main=sce1Title, breaks=20, col="grey80", ylab="Number of cells")
#   print(paste("Mean Library Size (millions): ", mean(sce1$total_counts/1e6)))
#   
#   #Number of Expressed Genes
#   hist(sce1$total_features, xlab="Number of expressed genes", main=sce1Title, breaks=20, col="grey80", ylab="Number of cells")
#   print(paste("Mean Number of Expressed Genes: ", mean(sce1$total_features)))
#   
#   #Number of Nonzero Cells
#   hist(rowSums(counts(sce1) != 0), xlab="Number of active cells", main=sce1Title, breaks=20, col="grey80", ylab="Number of genes")
#   print(paste("Mean Number of Active Cells: ", mean(rowSums(counts(sce1) != 0))))
#   
#   #Compare Gene Variances
#   var_sim = apply(counts(sce1), 1, var)
#   hist(log10(var_sim), xlab=expression(Log[10]~"Variance of Gene Expression"), main=sce1Title, breaks=20, col="grey80", ylab="Number of genes")
#   print(paste("Mean log10 Variance by Gene: ", mean(log10(var_sim))))
#   
#   #Compare Mean-Variance Trend
#   var.fit_sim <- trendVar(sce1, trend="loess", use.spikes=FALSE, span=0.2)
#   var.out_sim <- decomposeVar(sce1, var.fit_sim)
#   plot(var.out_sim$mean, var.out_sim$total, pch=16, cex=0.6, main=sce1Title, xlab="Mean log-expression", 
#        ylab="Variance of log-expression")
#   o <- order(var.out_sim$mean)
#   lines(var.out_sim$mean[o], var.out_sim$tech[o], col="dodgerblue", lwd=2)
# }

inspectData <- function(sce) {
  sce = computeSumFactors(sce)
  counts(sce) = t(t(counts(sce)) / pData(sce)$size_factor)
  
  par(mfrow=c(2,2))
  for (i in 1:100)
  {
    hist(log(counts(sce)[i,]), breaks=20, main="", col="grey80",
         xlab="Gene count value", ylab="Number of cells")
    Sys.sleep(2)
  }
  
  for (j in 100:200)
  {
    #hist(counts(sce)[,j], breaks=100, main="", col="grey80",
    #     xlab="Cell count value", ylab="Number of Genes", xlim=c(0,50))
    plot((counts(sce)[1:100,j+1] - counts(sce)[1:100,j]))
    Sys.sleep(2)
  }
}

compareAnalyses <- function(sce_sim, sce_real) {
  sce = sce_sim
  sce_sim = calculateQCMetrics(sce_sim)
  
  par(mfrow=c(1,2))
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
  
  pca = plotPCA(sce_sim, pca_data_input="pdata") + fontsize
  plot(pca)
  
  pca = plotPCA(sce_real, pca_data_input="pdata") + fontsize
  plot(pca)
  
  #####
  
  tsne = plotTSNE(sce_sim, perplexity=30, rand_seed=100) + fontsize
  plot(tsne)
  
  tsne = plotTSNE(sce_real, perplexity=30, rand_seed=100) + fontsize
  plot(tsne)
  
  #####
  
  plotQC(sce, type = "highest-expression", n=50) + fontsize
  
  numcells <- nexprs(sce, byrow=TRUE)
  smoothScatter(log10(ave.counts), numcells, xlab=expression(Log[10]~"average count"), 
                ylab="Number of expressing cells")
  sce <- sce[keep,]
  
  #Analyze mean-variance trend
  var.fit <- trendVar(sce, trend="loess", use.spikes=FALSE, span=0.2)
  var.out <- decomposeVar(sce, var.fit)
  plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
       ylab="Variance of log-expression")
  o <- order(var.out$mean)
  lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
  
  #Identify Highly Variable Genes
  hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
  hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
  nrow(hvg.out)
  
  write.table(file="hsc_hvg.tsv", hvg.out, sep="\t", quote=FALSE, col.names=NA)
  head(hvg.out)
  
  plotExpression(sce, rownames(hvg.out)[1:10]) + fontsize
  
  #Find correlations between highly variable genes
  set.seed(100)
  var.cor <- correlatePairs(sce, subset.row=rownames(hvg.out))
  write.table(file="hsc_cor.tsv", var.cor, sep="\t", quote=FALSE, row.names=FALSE)
  head(var.cor)
  
  sig.cor <- var.cor$FDR <= 0.05
  summary(sig.cor)
  
  #Make a graph between highly correlated genes and display clusters
  g <- ftM2graphNEL(cbind(var.cor$gene1, var.cor$gene2)[sig.cor,], 
                    W=NULL, V=NULL, edgemode="undirected")
  cl <- highlyConnSG(g)$clusters
  cl <- cl[order(lengths(cl), decreasing=TRUE)]
  head(cl)
  
  #Compute p-values for each gene to see if it is correlated to at least one other gene
  #Doesn't work right now
  var.cor <- correlatePairs(sce, subset.row=rownames(hvg.out), per.gene=TRUE) #use per.gene = TRUE
  head(var.cor)
  
  sig.cor <- var.cor$FDR <= 0.05
  summary(sig.cor)
  
  #Make a correlation matrix
  chosen <- unique(var.cor$gene1[sig.cor])
  norm.exprs <- exprs(sce)[chosen,,drop=FALSE]
  
  heat.vals <- norm.exprs - rowMeans(norm.exprs)
  heat.out <- heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.6)
  
  #Plot a PCA plot based on the correlation matrix from the HVG's
  plotPCA(sce, exprs_values="exprs", colour_by="total_features",
          feature_set=chosen) + fontsize
}

countZeroRows <- function(sce) {
  zeroRows = which(rowSums(counts(sce)) == 0)
  return(length(zeroRows))
}

showMeanVarTrend <- function(sce) {
  sce <- sce[which(rowSums(counts(sce)) > 0)] #Not original indices now
  ave_counts <- rowMeans(counts(sce))
  row_var <- apply(counts(sce), 1, var)
  var_outliers <- (isOutlier(row_var, nmads=3, type="higher", log=TRUE))
  mean_outliers <- (isOutlier(ave_counts, nmads=3, type="higher", log=TRUE))
  row_var = row_var[!(var_outliers | mean_outliers)]
  ave_counts = ave_counts[!(var_outliers | mean_outliers)]
  log_ave_counts = log(ave_counts)
  log_row_var = log(row_var)
  
  gene_dataframe = data.frame(log_ave_counts, log_row_var)
  fun.1 <- function(x) x
  plot(ggplot(aes(x=log_ave_counts, y=log_row_var), data=gene_dataframe) + geom_point() + stat_smooth(color="blue") + ggtitle("Log Gene Variance vs Log Gene Mean") + stat_function(fun = fun.1))
  
  log_var_resid = log_row_var - log_ave_counts
  var_resid_dataframe = data.frame(log_ave_counts, log_var_resid)
  plot(ggplot(aes(x=log_ave_counts, y=log_var_resid), data=var_resid_dataframe) + geom_point() + stat_smooth(color="blue") + ggtitle("Log Gene Variance Residual vs Log Gene Mean"))
  
  paste(lm(row_var ~ ave_counts))
  meanVarCorr = cor(x=log_ave_counts, y=log_row_var)
  paste("Log Mean-Var Correlation: ", meanVarCorr)
}