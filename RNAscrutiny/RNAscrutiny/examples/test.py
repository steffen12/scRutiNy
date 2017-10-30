import numpy as np
import random
import RNAscrutiny as scrutiny

saveFolder = './test'
try:
    os.makedirs(saveFolder)
except:
    pass

# Generate correlation matrix
n = 500
networkStructure = 'SimulatedPowerLaw'
n, gc, tf = scrutiny.generateGeneCorrelationMatrix(
    n=500, networkStructure=networkStructure, saveFolder=saveFolder)

# Form cell type tree
cells = 200
cellTypeTree, projMatrix, constMatrix = scrutiny.formCellTypeTree(
    n,
    cells=cells,
    cellTypeParents=[-1],
    cellTypeNumCells=[cells],
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

cellSizeFactors = np.zeros(shape=cells)
timeStep = 0.01
pseudotimes = timeStep * totalCellIterations
scrutiny.analyzeDataBeforeTransformation(
    n,
    cells,
    finalCellStates,
    alpha,
    pseudotimes,
    networkStructure=networkStructure,
    saveFolder=saveFolder)

finalCellStates, cellSizeFactors = scrutiny.transformData(
    n, cells, finalCellStates)

# This calls r commands
# saveData(finalCellStates, cellTypesRecord, gc, projMatrix, constMatrix,
#          cellSizeFactors, pseudotimes, saveFolder, targetDirectory,
#          R_Directory)

scrutiny.analyzeSCRNAseqData(
    n, cells, cellTypesRecord, finalCellStates, saveFolder=saveFolder)
