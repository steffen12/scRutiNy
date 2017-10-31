import numpy as np
import RNAscrutiny as scrutiny

import datetime
ts = datetime.datetime.now().strftime('%d_%m_%Y-%H:%M:%S')

# set output directory and seed random number generator
scrutiny.init(outputDir=ts, verbose=True)

# Generate correlation matrix and transcription factors
# parameter networkStructure should be one of
# - 'SimulatedPowerLaw'
# - 'TRRUST' (default)
# - 'Random'
# n: number of genes
networkStructure = 'SimulatedPowerLaw'
n = 700
gc, tf = scrutiny.generateGeneCorrelationMatrix(
    n=n, networkStructure=networkStructure)

# Form cell type tree
cells = 200
cellTypeTree, projMatrix, constMatrix = scrutiny.formCellTypeTree(
    n=n,
    cells=cells,
    cellTypeParents=[-1],
    cellTypeNumCells=[200],
    cellTypeConstProps=[0.0],
    transcriptionFactors=tf,
    cellTypeMeans=[0])

# Determine alpha based on optimal development for progenitor cell type (parent = -1)
initialStates = scrutiny.generateInitialStates(n, cells)
alpha = scrutiny.findOptimalAlpha(
    n, gc, cellTypeTree, initialStates, numToSample=50)
gc /= alpha  #Update gc with new alpha

# Recurse over cell type tree to develop cells
currentStates = initialStates
cellTypesRecord = np.zeros(shape=cells, dtype="int64")
totalCellIterations = np.zeros(shape=cells, dtype="int64")
scrutiny.findAllNextStates(gc,
                           cellTypeTree.head(), cellTypeTree, currentStates,
                           cellTypesRecord, totalCellIterations)
finalCellStates = currentStates

# output plots to outputDir
cellSizeFactors = np.zeros(shape=cells)
timeStep = 0.01
pseudotimes = timeStep * totalCellIterations
scrutiny.analyzeDataBeforeTransformation(
    n,
    cells,
    finalCellStates,
    alpha,
    pseudotimes,
    networkStructure=networkStructure)

# get final cell states
finalCellStates, cellSizeFactors = scrutiny.transformData(
    n, cells, finalCellStates)

# output plots to outputDir
scrutiny.analyzeSCRNAseqData(n, cells, cellTypesRecord, finalCellStates)
