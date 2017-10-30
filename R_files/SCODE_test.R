run_SCODE <- function(fdata, ftime, dir, tfnum, pnum, cnum, maxite) { #Transform from command-line arguments
  library(MASS)
  
  maxB <- 2.0
  minB <- -10.0
  
  system(paste("mkdir", dir, sep=" "))
  
  X <- as.matrix(read.table(fdata, sep="\t"))[1:tfnum,1:cnum]
  W <- matrix(rep(0,tfnum*pnum), nrow=tfnum, ncol=pnum)
  Z <- matrix(rep(0,pnum*cnum), nrow=pnum, ncol=cnum)
  WZ <- matrix(nrow=tfnum, ncol=cnum)
  
  #read pseudo-time and normalize pseudo-time
  pseudotime <- read.table(ftime, sep="\t")[1:cnum,2]
  pseudotime <- pseudotime/max(pseudotime)
  
  new_B <- rep(0, pnum)
  old_B <- rep(0, pnum)
  
  #initialization
  RSS <- Inf
  for(i in 1:pnum){
    new_B[i] <- runif(1, min=minB, max=maxB)
    old_B[i] <- new_B[i]
  }
  
  #function to sample Z
  sample_Z <- function(){
    for(i in 1:pnum){
      for(j in 1:cnum){
        Z[i,j] <<- exp(new_B[i]*pseudotime[j]) + runif(1, min=-0.001, max=0.001)
      }
    }
  }
  
  #optimize W and B iteratively
  for(ite in 1:maxite){
    #sampling B
    target <- floor(runif(1, min=1, max=pnum+1))
    new_B[target] <- runif(1, min=minB, max=maxB)
    
    #for last calculation
    if(ite == maxite){
      for(i in 1:pnum){
        new_B[i] <- old_B[i]
      }
    }
    
    #sample Z from new B
    sample_Z()
    
    #regression
    for(i in 1:tfnum){
      X.lm <- lm(X[i,] ~ t(Z)-1)
      for(j in 1:pnum){
        W[i,j] <- X.lm$coefficients[j]
      }
      WZ[i,] <- W[i,] %*% Z
    }
    
    #RSS
    tmp_RSS <- sum((X-WZ)**2)
    if(tmp_RSS < RSS){
      RSS <- tmp_RSS
    }
    else{
      new_B[target] <- old_B[target]
    }
  }
  
  #output RSS
  write.table(RSS, paste(dir,"/RSS.txt",sep=""), row.names=F, col.names=F, sep="\t")
  
  #output W
  write.table(W, paste(dir,"/W.txt",sep=""), row.names=F, col.names=F, sep="\t")
  
  #infer A
  B <- matrix(rep(0,pnum*pnum), nrow=pnum, ncol=pnum)
  for(i in 1:pnum){
    B[i,i] <- new_B[i]
  }
  invW <- ginv(W)
  A <- W %*% B %*% invW
  
  #output A and B
  write.table(A, paste(dir,"/A.txt",sep=""), row.names=F, col.names=F, sep="\t")
  write.table(B, paste(dir,"/B.txt",sep=""), row.names=F, col.names=F, sep="\t")
}

graphGenePseudotime <- function(counts, pseudotimes, gene_num) {
  library(ggplot2)
  
  gene_counts <- counts[gene_num, ]
  title = paste("Gene", gene_num, "Counts over Pseudotime")
  gene_dataframe = data.frame(pseudotimes, gene_counts)
  p = ggplot(aes(x=pseudotimes, y=gene_counts), data=gene_dataframe, main=title) + geom_point() + geom_line() + stat_smooth(color="blue") + ggtitle(title)
  return(p)
}

convertTimeSeriesDataToSForm <- function(sce) {
  num_genes = dim(counts(sce))[1]
  num_cells = dim(counts(sce))[2]
  
  S = exprs(sce) #Just to initialize matrix with right size
  
  sce = calculateQCMetrics(sce)
  sce = computeSumFactors(sce, positive=TRUE)
  normalize(sce)
  
  for(gene in 1:num_genes) {
    rowRanks = rank(exprs(sce)[gene, ])
    S[gene, ] = 2 * (rowRanks / (num_cells + 1)) - 1
  }
  
  #for(cell in 1:num_cells) {
  #  cellRanks = rank(exprs(sce)[, cell])
  #  S[, cell] = 2 * (cellRanks / (num_genes + 1)) - 1
  #}
  
  sce_counts = log(exprs(sce))
  
  sce_row_means = rowMeans(sce_counts, na.rm = TRUE)
  sce_row_var = rowVar(sce_counts)
  
  sce_counts = sce_counts - rowMeans(sce_counts)
  sce_counts = sce_counts / sce_row_var
  
  #for(gene_num in 1:num_genes) {
  #  p = graphGenePseudotime(S, pseudotimes, gene_num)
  #  plot(p)
  #  Sys.sleep(7)
  #}
  
  save(S, file="S.RData")
  return(S)
}

write_SCODE_test_data_to_python <- function() {
  root <- "/home/scornwell/Internship/R_files/SCODE-master"
  home_dir = "/home/scornwell"
  setwd(paste(home_dir, "Internship", sep="/"))
  data_dir = "data2"
  data_file <- paste(root, data_dir, "data.txt", sep="/")
  pseudotime_file <- paste(root, data_dir, "time.txt", sep="/")
  tf_file <- paste(root, data_dir, "tf.txt", sep="/")
  output_directory <- paste(root, data_dir, "predicted_SCODE", sep="/")
  num_transcription_factors = 100
  num_z = 4
  num_cells = 373
  max_iterations = 1000
  pseudotimeScalar = 4
  
  #run_SCODE(data_file, pseudotime_file, output_directory, num_transcription_factors, num_z, num_cells, max_iterations)
  
  counts <- as.matrix(read.table(data_file, sep="\t"))[1:num_transcription_factors,1:num_cells]
  sce <- newSCESet(countData=counts)
  pseudotimes <- read.table(pseudotime_file, sep="\t")[1:num_cells,2]
  pseudotimes <- pseudotimeScalar*(pseudotimes/max(pseudotimes))
  transcriptionFactors <- read.table(tf_file)
  S = convertTimeSeriesDataToSForm(sce)
  
  regulator_gene = "Myc"
  regulated_gene = "Plagl1"
  regulator_gene_index = which(transcriptionFactors == regulator_gene)
  regulated_gene_index = which(transcriptionFactors == regulated_gene)
  p1 <- graphGenePseudotime(S, pseudotimes, regulator_gene_index)
  p2 <- graphGenePseudotime(S, pseudotimes, regulated_gene_index)
  plot_grid(p1, p2, align="h")
  
  write.table(S, file='S_cell_states.tsv', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
  write.table(pseudotimes, file='realPseudotimes.tsv', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
  write.table(transcriptionFactors, file='transcriptionFactors.txt', quote=FALSE, col.names = FALSE, row.names = FALSE)
}

graphSimulatedGenePseudotimes <- function() {
  sce_sim = load_sce_sim()
  S = convertTimeSeriesDataToSForm(sce_sim)
  for(gene_num in 1:5) {
    graphGenePseudotime(S, pData(sce_sim)$synthPseudotimes, gene_num)
    ggsave(paste("Pseudotime_Gene", gene_num, ".tiff", sep=""), plot = last_plot(), device="tiff")
  }
}