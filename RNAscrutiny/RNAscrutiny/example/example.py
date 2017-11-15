import numpy as np

# import the generative model module and initialize
# the output directory is scRutiNy/MatriSeq by default
# but can be specified by passing parameter outputDir to init
import RNAscrutiny.MatriSeq as ms
ms.init(verbose=True)

# generate gene correlation (gc) matrix
# and transcription factors (tf) i.e. genes with positive outdegree
# networkStructure should be one of
# - 'SimulatedPowerLaw'
# - 'TRRUST' (default)
# - 'Random'
# optional parameters to generateGeneCorrleationMatrix:
# - maxAlphaBeta (default: 10)
# - normalizeWRows (default: False)
# - powerLawExponent (default: 2.05)
# - networkSparsity (default: 0.999)
networkStructure = 'SimulatedPowerLaw'
n = 500  # number of genes
gc, tf = ms.generateGeneCorrelationMatrix(
    n=n, networkStructure=networkStructure)

# form a tree to represent cell types
cells = 200
cellTypeTree, projMatrix, constMatrix = ms.formCellTypeTree(
    n=n,
    cells=cells,
    cellTypeParents=[-1, -1],
    cellTypeNumCells=[100, 100],
    cellTypeConstProps=[0.0, 0.0],
    transcriptionFactors=tf,
    cellTypeMeans=[0, 0])

# rescale gc to ensure convergence
# optional paramters to finOptimalAlpha
# - timeStep (default: 0.01)
# - convergenceThreshold (default: 0.1)
# - noiseDisp (default: 0.05)
# - maxNumIterations (default: 500)
# - convDistAvgNum (default: 20)
# - initialAlpha (default: 0.01)
# - alphaStep (default: 0.02)
# - showAlphaGeneMeans (default: False)
# - geneMeanLevels (default: 0)
initialStates = ms.generateInitialStates(n, cells)
alpha = ms.findOptimalAlpha(n, gc, cellTypeTree, initialStates, numToSample=50)
gc /= alpha

# recurse over cell type tree to develop cells
# optional parameters to findAllNextStates:
# - timeStep (default: 0.01)
# - noiseDisp (default: 0.05)
# - convergenceThreshold (default: 0.1)
# - maxNumIterations (default: 500)
# - convDistAvgNum (default: 20)
# - cellsToSample (default: 50)
# - showPlots (default: False)
# - cellDevelopmentMode (default: True)
# - geneMeanLevels (default: 0)
currentStates = initialStates
cellTypesRecord = np.zeros(shape=cells, dtype="int64")
totalCellIterations = np.zeros(shape=cells, dtype="int64")
ms.findAllNextStates(gc,
                     cellTypeTree.head(), cellTypeTree, currentStates,
                     cellTypesRecord, totalCellIterations)
finalCellStates = currentStates

# get cells' final states, plot and save data
# optional paramter to analyzeSCRNAseqData:
# - perplexityParam (default: 30)
cellSizeFactors = np.zeros(shape=cells)
timeStep = 0.01
pseudotimes = timeStep * totalCellIterations
ms.analyzeDataBeforeTransformation(
    n,
    cells,
    finalCellStates,
    alpha,
    pseudotimes,
    networkStructure=networkStructure)
finalCellStates, cellSizeFactors = ms.transformData(n, cells, finalCellStates)
ms.analyzeSCRNAseqData(n, cells, cellTypesRecord, finalCellStates)
ms.saveData(finalCellStates, cellTypesRecord, gc, projMatrix, constMatrix,
            cellSizeFactors, pseudotimes)

# import inference module and initialize
# by default, the output directory is RNAscrutiny/RegNetInference
# but can be specified by passing parameter outputDir to init
# init expects the input directory to be RNAscrutiny/MatriSeq
# but can be specified by passing parameter inputDir to init
import RNAscrutiny.RegNetInference as rni
rni.init(verbose=True)

# load cell states from input directory
# rescale pseudotime, smooth cell states,
# plot gene expression over pseudotime
n, cells, cellStates, pseudotimes, cellTypesRecord, W_real = rni.readCellStates(
)
pseudotimes /= np.max(pseudotimes)
cellStates = rni.smoothCellData(cells, cellStates, pseudotimes)
rni.plotGenesOverPseudotime(n, pseudotimes, cellStates)

# cluster genes
# optional parameters to clusterGenes
# - perplexity (default: 30)
# - epsilon (default: 0.05)
# - min_samples (default: 2)
geneClusters = rni.clusterGenes(n, cells, cellStates)

# find optimal gene correlation matrix W, save, and plot
# optional parameters to crossValidateData:
# - numPartitions (default: 5)
# - alphaMin (default: 1e-4)
# - alphaMax (default: 1e-3)
# - alphaMultStep (default: 0.01)
# - pseudotimeThresh (default: 0.1)
optW = rni.crossValidateData(n, cells, cellStates, pseudotimes,
                             cellTypesRecord, W_real, geneClusters)
rni.saveW(optW)
rni.plotWGraphs(n, optW, W_real)
