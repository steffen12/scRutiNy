import numpy as np
import RNAscrutiny.MatriSeq as ms

import datetime
ts = datetime.datetime.now().strftime('%d_%m_%Y-%H:%M:%S')

# set output directory and seed random number generator
ms.init(outputDir=ts, verbose=True)

# Generate correlation matrix and transcription factors
# parameter networkStructure should be one of
# - 'SimulatedPowerLaw'
# - 'TRRUST' (default)
# - 'Random'
# n: number of genes
networkStructure = 'SimulatedPowerLaw'
n = 700
gc, tf = ms.generateGeneCorrelationMatrix(
    n=n, networkStructure=networkStructure)

# Form cell type tree
cells = 200
cellTypeTree, projMatrix, constMatrix = ms.formCellTypeTree(
    n=n,
    cells=cells,
    cellTypeParents=[-1],
    cellTypeNumCells=[200],
    cellTypeConstProps=[0.0],
    transcriptionFactors=tf,
    cellTypeMeans=[0])

# Determine alpha based on optimal development for progenitor cell type (parent = -1)
initialStates = ms.generateInitialStates(n, cells)
alpha = ms.findOptimalAlpha(n, gc, cellTypeTree, initialStates, numToSample=50)
gc /= alpha  #Update gc with new alpha

# Recurse over cell type tree to develop cells
currentStates = initialStates
cellTypesRecord = np.zeros(shape=cells, dtype="int64")
totalCellIterations = np.zeros(shape=cells, dtype="int64")
ms.findAllNextStates(gc,
                     cellTypeTree.head(), cellTypeTree, currentStates,
                     cellTypesRecord, totalCellIterations)
finalCellStates = currentStates

# output plots to outputDir
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

# get final cell states
finalCellStates, cellSizeFactors = ms.transformData(n, cells, finalCellStates)

# output plots to outputDir
ms.analyzeSCRNAseqData(n, cells, cellTypesRecord, finalCellStates)
