find_pseudotimes_slingshot <- function(sce) {
  library("slingshot")
  pca <- prcomp(t(counts(sce)), center=TRUE, scale=TRUE)
  plot(pca$x, pch=19)
  pcNum = 2
  for(cellType in 0:max(pData(sce)$cellType)) {
    cellTypeCells = which(pData(sce)$cellType == cellType)
    plot(main = cellType, pca$x[cellTypeCells, pcNum], pData(sce)$totalCellIterations[cellTypeCells], pch=19)
    cor(pca$x[cellTypeCells, pcNum], pData(sce)$totalCellIterations[cellTypeCells])
  }
}

getInitialState <- function(sce) {
  cellTypeStateTable = table(pData(sce)$State, pData(sce)$cellType)
  initialState = which(cellTypeStateTable == max(cellTypeStateTable[, 1]))
  print(cellTypeStateTable)
  return(initialState)
}

order_by_pseudotime <- function(sce, simulated) {
  mon_sce <- newCellDataSet(counts(sce), phenoData = phenoData(sce), featureData = featureData(sce), expressionFamily = negbinomial.size())
  mon_sce <- estimateSizeFactors(mon_sce)
  mon_sce <- estimateDispersions(mon_sce)
  disp_table <- dispersionTable(mon_sce)
  ordering_genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
  mon_sce <- setOrderingFilter(mon_sce, ordering_genes)
  plot_ordering_genes(mon_sce)
  
  #plot_pc_variance_explained(mon_sce, return_all = F)
  
  mon_sce <- reduceDimension(mon_sce, max_components = 2)
  
  mon_sce <- orderCells(mon_sce)
  #mon_sce <- clusterCells(mon_sce)
  #plot_cell_clusters(mon_sce)
  if(simulated) {
    initialState = getInitialState(mon_sce)
    mon_sce <- orderCells(mon_sce, root_state=initialState)
  }
  
  #plot(plot_cell_trajectory(mon_sce))
  return(mon_sce)
  
  #pData(mon_sce)$State = factor(pData(sce)$cellType)
  #plot_cell_trajectory(mon_sce)
}

analyze_cell_branches <- function(sce) {
  #nonZeroRows = which(rowSums(counts(sce)) != 0)
  #names(nonZeroRows) = NULL
  #sce = sce[which(rowSums(counts(sce)) != 0), ]
  
  BEAM_res <- BEAM(sce, branch_point = 1)
  BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
  
  plot_genes_branched_heatmap(sce[row.names(subset(BEAM_res, qval < 0.05)),], branch_point = 1, num_clusters = 2, show_rownames = T)
  
  #branch_to_analyze = 2
  #plot_genes_branched_pseudotime(sce[branch_genes[[branch_to_analyze]],], branch_point = 1, color_by = "Time", ncol = 1)
  
  analysis_branch_genes <- as.integer(rownames(BEAM_res[which(BEAM_res$qval < 0.05), ]))
  return(analysis_branch_genes);
}

analyze_pseudotimes_monocle <- function(sce, simulated, directory) {
  library("monocle")
  library("ggplot2")
  library("gplots")
  outputDir <- paste(directory, "Monocle_Analysis", sep="/")
  dir.create(outputDir)

  #sce = computeSumFactors(sce, positive=TRUE)
  #sce = normalize(sce)
  #counts(sce) = t(t(counts(sce)) / pData(sce)$size_factor)
  #cellStates = counts(sce)
  
  sce = order_by_pseudotime(sce, simulated)
  cellTraj <- plot_cell_trajectory(sce)
  pdf(paste(outputDir, 'Monocle_Cell_Trajectory.pdf', sep="/"))
  plot(cellTraj)
  dev.off()
  pseudotimes = pData(sce)$Pseudotime
  
  BEAM_res <- BEAM(sce, branch_point = 1)
  BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
  
  pdf(paste(outputDir, 'BEAM_Analysis.pdf', sep="/"))
  plot_genes_branched_heatmap(sce[row.names(subset(BEAM_res, qval < 0.05)),], branch_point = 1, num_clusters = 3, show_rownames = F)
  dev.off()
  
  analysis_branch_genes <- as.integer(rownames(BEAM_res[which(BEAM_res$qval < 0.05), ]))
  
  if(length(pseudotimes) != 0 && simulated) {
    branch_genes <- get_branch_1_genes() #Load genes determining the branch
    trueDiffExprGenes <- union(branch_genes[[2]], branch_genes[[3]])
    
    cellTypeStateTable = table(pData(sce)$State, pData(sce)$cellType)
    write.table(cellTypeStateTable, paste(outputDir, "cellStateTable.tsv", sep="/"), sep="\t")
    
    pdf(paste(outputDir, 'Branch_Genes_Monocle_Venn.pdf', sep="/"))
    plot(venn(list(True_Branch_Genes = trueDiffExprGenes, Monocle_Branch_Genes = analysis_branch_genes)))
    dev.off()
    
    onlySimGenes <- setdiff(trueDiffExprGenes, analysis_branch_genes)
    bothGenes <- intersect(trueDiffExprGenes, analysis_branch_genes)
    onlyMonocleGenes <- setdiff(analysis_branch_genes, trueDiffExprGenes)
    
    print(paste("Branch genes only found in simulation: ", onlySimGenes))
    print(paste("Branch genes found in both simulation and analysis: ", bothGenes))
    print(paste("Branch genes only found in Monocle analysis: ", onlyMonocleGenes))
    
    print(paste("Num Branch genes only found in simulation: ", length(onlySimGenes)))
    print(paste("Num Branch genes found in both simulation and analysis: ", length(bothGenes)))
    print(paste("Num Branch genes only found in Monocle analysis: ", length(onlyMonocleGenes)))
    
    save(pseudotimes, file="pseudotimes.RData")
    #save(cellStates, file="cellStates.RData")
    
    states = unique(pData(sce)$State)
    correlations = vector(mode="double", length=length(states))
    for(stateIndex in 1:length(states)) {
      state = states[stateIndex]
      currStateCells = which(pData(sce)$State == state)
      xlabel = paste("Actual Pseudotime - State", state)
      ylabel = paste("Monocle Pseudotime - State", state)
      
      simulatedPseudotimes = (pData(sce)$synthPseudotimes)[currStateCells]
      monoclePseudotimes = pseudotimes[currStateCells]
      pseudotime_dataframe = data.frame(simulatedPseudotimes, monoclePseudotimes)
      ggplot(aes(x=simulatedPseudotimes, y=monoclePseudotimes), data=pseudotime_dataframe) + geom_point() + stat_smooth(color="blue")
      ggsave(paste("Pseudotime_Corelation_State_", state, ".pdf", sep=""), plot = last_plot())
      
      correlation = cor((pData(sce)$synthPseudotimes)[currStateCells], pseudotimes[currStateCells], method="spearman")
      print(paste("State", state, "Pseudotime Correlation:", correlation))
      correlations[stateIndex] = correlation
    }
    write.table(correlations, paste(outputDir, "pseudotimeCorrelations.tsv", sep="/"), sep="\t")
  }
  plot(plot_genes_branched_pseudotime(sce[1:4,], branch_point = 1, ncol=1))
  return(pseudotimes)
}

find_pseudotimes_cellTree <- function(sce) {
  library("cellTree")
  
  lda.results = compute.lda(counts(sce), k.topics=3:8, method="maptpx")
  sce_lda_model = compute.lda(counts(sce), k.topics=lda.results$K)
  dists = get.cell.dists(sce_lda_model)
  mst.tree = compute.backbone.tree(sce_lda_model, only.mst = TRUE)
  mst.tree.with.layout = ct.plot.topics(mst.tree)
  b.tree = compute.backbone.tree(sce_lda_model)
  b.tree.with.layout = ct.plot.topics(b.tree)
}