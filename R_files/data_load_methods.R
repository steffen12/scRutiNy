#source("/home/scornwell/Internship/R_files/data_test_methods.R")

load_sce_sim <- function(directory) {
  load(paste(directory, "synthscRNAseq.gzip", sep="/"))
  #synthscRNAseq <- synthscRNAseq[rowSums(synthscRNAseq) > 0, ] #This is currently done to the real data too
  matriSeq = newSCESet(countData=synthscRNAseq)
  sce_sim = calculateQCMetrics(matriSeq)
  
  load(paste(directory, "cellTypesRecord.gzip", sep="/"))
  pData(sce_sim)$cellType = cellTypesRecord
  
  load(paste(directory, "cellSizeFactors.gzip", sep="/"))
  pData(sce_sim)$simCellSizeFactors = cellSizeFactors
  sce_sim$expName = "Simulated"
  
  load(paste(directory, "synthPseudotimes.gzip", sep="/"))
  pData(sce_sim)$synthPseudotimes = synthPseudotimes
  return(sce_sim);
}

get_branch_1_genes <- function() {
  load("/home/scornwell/Internship/R_files/projMatrix.gzip")
  
  branch_1_genes = which(projMatrix[1,] == 0)
  branch_2_genes = which(projMatrix[1,] != 0 & projMatrix[2,] == 0)
  branch_3_genes = which(projMatrix[1,] != 0 & projMatrix[3,] == 0)
  branch_genes = list(branch_1_genes, branch_2_genes, branch_3_genes)
  return(branch_genes);
}

load_sims <- function(datasetNum, matriSeq) {
  root <- "/home/scornwell/Internship/R_files/data"
  datasets <- read_tsv(file.path(root, "datasets.txt"),
                       col_types = cols(.default = col_character(),
                                        NumCells = col_integer()
                       )
  )
  
  dataset = datasets[datasetNum,]
  numGenes = dim(counts(matriSeq))[1]
  numCells = dim(counts(matriSeq))[2]
  sims = simCompDataset(dataset, root, numGenes, numCells)
  sce_real = sims$Real
  sce_real$expName = "Real"
  sims$Real = sce_real
  return(sims);
}

load_sce_real_from_file <- function(countFile) {
  dataFolder = "/home/scornwell/Internship/R_files/data"
  dataLocation = paste(dataFolder, countFile, sep="/")
  counts = read.table(dataLocation, header=TRUE)
  sce_real = newSCESet(countData=counts)
  sce_real = calculateQCMetrics(sce_real)
  return(sce_real);
}

compare_sims <- function(sims) {
  dataList = compareAllDatasets(sims)
  savePlots(dataList)
  dataList
  
  comp = dataList$Comp
  diff = dataList$Diff
  diffSum = summariseDiff(diff)
  totalRanks = cbind(diffSum$Ranks[1:2])
  for(i in 1:nrow(diffSum$Ranks)) {
    totalRanks[i, 2] = sum(diffSum$Ranks[i, 2:ncol(diffSum$Ranks)]) / (ncol(diffSum$Ranks) - 1)
  }
  orderedTotalRanks = totalRanks[order(subset(totalRanks, select=c("MeanRank"))), ]
  orderedTotalRanks
  return(list(Comp = comp, Diff = diff, AllRanks = diffSum$Ranks, OrderedTotalRanks = orderedTotalRanks))
}

load_sce_list <- function(matriSeq) {
  numGenes = dim(counts(matriSeq))[1]
  numCells = dim(counts(matriSeq))[2]
  
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
  sce6 = newSCESet(countData=loadDataset(datasets[6,], root))
  sce7 = newSCESet(countData=loadDataset(datasets[7,], root))
  
  #Special transformations
  sce1 = sce1[,2:dim(sce1)[2]] #First cell has way too many counts for some reason
  
  sce_list = list(sce1, sce2, sce3, sce4, sce5, sce6, sce7)
  
  for (x in 1:7) 
  {
    sce = sce_list[[x]]
    
    #COMMENT this out when I want the dimensions of the real data to match that of the simulated data
    #numCells = ncol(counts(sce))
    #numGenes = nrow(counts(sce))
    
    #if(dim(counts(sce))[2] > numCells) {
    #  sce <- sce[, sample(1:dim(counts(sce))[2], numCells)]
    #}
    #if(dim(counts(sce))[1] > numGenes) {
    #  sce <- sce[sample(1:dim(counts(sce))[1], numGenes), ]
    #}
    sce_list[[x]] <- calculateQCMetrics(sce)
  }
  
  return(sce_list);
}

#Adapted from Seurat Package
read10x <- function(data.dir = NULL){
  full_data <- list()
  for(i in seq_along(data.dir)){
    run <- data.dir[i]
    if (!dir.exists(run)){
      stop("Directory provided does not exist")
    }
    
    if(!grepl("\\/$", run)){
      run <- paste(run, "/", sep = "")
    }
    
    barcode.loc <- paste(run, "barcodes.tsv", sep ="")
    gene.loc <- paste(run, "genes.tsv", sep ="")
    matrix.loc <- paste(run, "matrix.mtx", sep ="")
    
    if (!file.exists(barcode.loc)){
      stop("Barcode file missing")
    }
    if (!file.exists(gene.loc)){
      stop("Gene name file missing")
    }
    if (!file.exists(matrix.loc)){
      stop("Expression matrix file missing")
    }
    
    data <- readMM(matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if(all(grepl("\\-1$", cell.names)) == TRUE) {
      cell.names <- as.vector(as.character(sapply(cell.names, function(cellName) strsplit(cellName, "-")[[1]][1])))
    }
    rownames(data) <- make.unique(as.character(sapply(gene.names, function(geneName) strsplit(geneName, "\\t")[[1]][2]))) 
    
    if(is.null(names(data.dir))){
      if(i < 2){
        colnames(data) <- cell.names
      }
      else {
        colnames(data) <- paste0(i, "_", cell.names, sep = "") 
      }
    } else {
      colnames(data) <- paste0(names(data.dir)[i],"_",cell.names) 
    }
    full_data <- append(full_data, data)
  }
  full_data <- do.call(cbind, full_data)
  full_data <- full_data[, which(colSums(full_data) > 0)]
  #sce_real <- newSCESet(countData=full_data)
  return(full_data)
}