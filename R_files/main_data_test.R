# RNA-seq
library("scater") #sce library
library("scran")
## library("splatter")
#library("ccRemover")
#library("RBGL")

# Parallel
#library("BiocParallel")

# Plotting
library("cowplot")
library("gplots")

# Tables
library("knitr")

# Tidyverse
library("tidyverse")

#Clustering
#library("dbscan")

library("Matrix")
#library("GEOquery")

#Directory where R files and objects are stored
directory = "/home/steffen12/NIH_Internship/R_files"

source(paste(directory, "load_datasets.R", sep="/"))
source(paste(directory, "simulate_datasets.R", sep="/"))
source(paste(directory, "utils.R", sep="/"))
source(paste(directory, "data_test_methods.R", sep="/"))
source(paste(directory, "reverse_eng_data.R", sep="/"))
source(paste(directory, "data_load_methods.R", sep="/"))
source(paste(directory, "clustering.R", sep="/"))
source(paste(directory, "pseudotimeInference.R", sep="/"))
source(paste(directory, "SCODE_test.R", sep="/"))

###############Parameters################
datasetNum = 3
numClusters = 8

################Load Data################
#Load sce_sim and sce_real
sce_sim = load_sce_sim(directory);
#sce_sim = sce_sim[which(rowSums(counts(sce_sim)) != 0), ]
sims = load_sims(datasetNum, sce_sim)
sce_real = sims$Real
sims$matriSeq = sce_sim;

#Load all sce_real
sce_list = load_sce_list(sce_sim);

#Load sce_sim
sce_sim = load_sce_sim();
sims$matriSeq = sce_sim;

#Run Splatter Comparison
sce_sim = load_sce_sim(directory);
sims = load_sims(datasetNum, sce_sim)
sce_real = sims$Real
sims$matriSeq = sce_sim; #NEED TO DO THIS FIRST!
compare_sims(sims);
###

countFile = "Myoblast_fpkm_matrix.txt"
sce_real = load_sce_real_from_file(countFile);

#sce_sim = calculateQCMetrics(sims$Splat)

############### Analyze Data ###################
#sce = qualityControl(sce);

#clusterIndices = generateClusters(sce, numClusters)

params <- splatter::splatEstimate(counts(sce_real))

inspectData(sce_sim)

#Plot distributions for simulated data
sce_sim = load_sce_sim(directory);
plotDistribution(sce_sim)
###

compareDistributions(sce_sim, sce_real)
+#compareClusters(sce, numClusters, clusterIndices)

compareAnalyses(sce_sim, sce_real)


branch_genes = get_branch_genes();
diff_branch_genes_list = unlist(branch_genes[2:length(branch_genes)])
analysis_branch_genes = analyze_cell_branches(sce_sim, branch_genes);
shared_genes = intersect(diff_branch_genes_list, analysis_branch_genes)
original_only = setdiff(diff_branch_genes_list, analysis_branch_genes)
analysis_only = setdiff(analysis_branch_genes, diff_branch_genes_list)

sce_clustered = cluster_CIDR(sce)
compareClusters(sce_clustered)

write_SCODE_test_data_to_python()
for(gene_num in 1:100) {
  graphGenePseudotime(exprs(sce), pseudotimes, gene_num)
  Sys.sleep(7)
}

#Monocle Analysis
sce_sim = load_sce_sim(directory);
simulated = TRUE
pseudotimes = analyze_pseudotimes_monocle(sce_sim, simulated, directory)
###


simulated = FALSE
pseudotimes = analyze_pseudotimes_monocle(sce_real, simulated)

reverse_engineer_data(sce_real)
reverse_engineer_data(sce_sim)

reverse_engineer_input_parameters(sce_real)

test_row_transformations(sce_sim)
sce_subset = sce_sim[,which(pData(sce_sim)$cellType != 0)]
