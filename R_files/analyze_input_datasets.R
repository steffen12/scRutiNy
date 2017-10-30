# RNA-seq
library("scater")
library("scran")
library("splatter")
library("ccRemover")
library("RBGL")

# Parallel
library("BiocParallel")

# Plotting
library("cowplot")
library("gplots")

# Tables
library("knitr")

# Tidyverse
library("tidyverse")

library("GEOquery")

source("/home/scornwell/Internship/R_files/load_datasets.R")
source("/home/scornwell/Internship/R_files/simulate_datasets.R")
source("/home/scornwell/Internship/R_files/utils.R")

root <- "/home/scornwell/Internship/R_files/data"
datasets <- read_tsv(file.path(root, "datasets.txt"),
                     col_types = cols(.default = col_character(),
                                      NumCells = col_integer()
                     )
)

sce1 = newSCESet(countData=loadDataset(datasets[1,], root))
sce2 = newSCESet(countData=loadDataset(datasets[2,], root))
sce3 = newSCESet(countData=loadDataset(datasets[3,], root))
sce4 = newSCESet(countData=loadDataset(datasets[4,], root))
sce5 = newSCESet(countData=loadDataset(datasets[5,], root))
sce_list = list(sce1, sce2, sce3, sce4, sce5)

for (x in 1:5) 
{
  sce = sce_list[[x]]
  print(sce)
  
  sce <- calculateQCMetrics(sce)
  libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="higher", log=TRUE)
  feature.drop <- isOutlier(sce$total_features, nmads=3, type="higher", log=TRUE)
  sce <- sce[,!(libsize.drop | feature.drop)]
  data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), Remaining=ncol(sce))
  
  sce_list[[x]] <- sce
}

for (x in 1:5)
{
  sce = sce_list[[x]]
  
  ave.counts <- rowSums(counts(sce))/dim(counts(sce))[2]
  keep <- ave.counts >= 1
  sum(keep)
  
  par(mfrow=c(2,2))
  for (i in 1000:1100)
  {
    hist(log(counts(sce))[i,], breaks=100, main="", col="grey80",
         xlab="Gene count value", ylab="Number of cells")
    Sys.sleep(2)
  }
  hist(counts(sce)[200,], breaks=100, main="", col="grey80",
       xlab="Gene count value", ylab="Number of cells")
  hist(counts(sce)[1000,], breaks=100, main="", col="grey80",
       xlab="Gene count value", ylab="Number of cells")
  hist(counts(sce)[2000,], breaks=100, main="", col="grey80",
       xlab="Gene count value", ylab="Number of cells")
  Sys.sleep(10)
  
  par(mfrow=c(1,1))
  hist(log10(ave.counts), breaks=100, main="", col="grey80",
       xlab=expression(Log[10]~"average count"))
  abline(v=log10(1), col="blue", lwd=2, lty=2)
  Sys.sleep(10)
  
  par(mfrow=c(1,2))
  hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", 
       breaks=20, col="grey80", ylab="Number of cells")
  hist(sce$total_features, xlab="Number of expressed genes", main="", 
       breaks=20, col="grey80", ylab="Number of cells")
  Sys.sleep(10)
}

fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))

sce = sce_list[[1]]
plotPCA(sce, pca_data_input="pdata") + fontsize

sce = sce_list[[2]]
plotPCA(sce, pca_data_input="pdata") + fontsize

sce = sce_list[[3]]
plotPCA(sce, pca_data_input="pdata") + fontsize

sce = sce_list[[4]]
plotPCA(sce, pca_data_input="pdata") + fontsize

sce = sce_list[[5]]
plotPCA(sce, pca_data_input="pdata") + fontsize

#####

sce = sce_list[[1]]
plotTSNE(sce, perplexity=10, rand_seed=100) + fontsize

sce = sce_list[[2]]
plotTSNE(sce, perplexity=10, rand_seed=100) + fontsize

sce = sce_list[[3]]
plotTSNE(sce, perplexity=10, rand_seed=100) + fontsize

sce = sce_list[[4]]
plotTSNE(sce, perplexity=10, rand_seed=100) + fontsize

sce = sce_list[[5]]
plotTSNE(sce, perplexity=10, rand_seed=100) + fontsize

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
