from scipy.special import expit
from scipy.misc import derivative
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from rpy2.robjects import r
from rpy2.robjects.numpy2ri import numpy2ri
from collections import Counter
from numpy import arange
from tree import Tree

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

import random
import scipy
import time
import math
import os
import datetime

def main():
	###Define Input Parameters###

	###Step 0 - choose a directory to work in###
	targetDirectory = "/home/steffen12/NIH_Internship/"
	R_Directory = "/home/steffen12/NIH_Internship/R_files/"
	os.chdir(targetDirectory)

	###Step 1 - Select the parameters to determine W###

	#n is the number of genes, W will be a nxn matrix
	n = 700 #19027
	#cells = 100 #864

	#SimulatedPowerLaw - Simulated power law network with powerLawExponent as parameter
	#TRRUST - Real network from TRRUST database of transcription factor interactions
	#Random - Random network with networkSparsity as a parameter

	networkStructure = "SimulatedPowerLaw"

	powerLawExponent = 2.05 #Increase for more connections, but must stay between 2-3
	networkSparsity = 0.999

	###Step 2 - Choose the minimum value of alpha (alpha must be greater than or equal to 0)#
	#and the maximum value of the Beta distribution alpha and beta parameters###

	initialAlpha = 0.01 #Set to 0 otherwise
	maxAlphaBeta = 10

	###Choose whether to divide each row of W by the number of connections:###

	normalizeWRows = False

	###Step 3 - Choose the cell development structure:###

	sameInitialCell = True
	cellDevelopmentMode = True

	cellTypeNumCells = [200] #[100, 100, 100, 300, 300, 300, 300]
	cellTypeParents = [-1]#[-1, 0, 0, 1, 1, 2, 2]
	cellTypeConstProps = [0.00] #[0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05] #Inherited
	cellTypeMeans = [0] #Not inherited, feature not available in current version

	###Step 4 - Choose the timeStep, noise, convergence, and iteration sampling parameters:###

	#Controlling gene trajectory smoothness
	timeStep = 0.01 #0.01
	noiseDisp = 0.05

	#Parameters to calculate alpha
	convergenceThreshold = 0.1
	alphaStep = 0.02 #How much to increase alpha by in each step - DO NOT decrease below 0.01
	maxNumIterations = 500
	cellsToSample = 50 #How many cells to sample - must be less than num cells
	convDistAvgNum = 20

	###Step 5 - Choose the mean and cell size mapping parameters:###

	#Select exponential base
	expScale = 1#4#6.81

	#Select whether to use a Gamma distributed gene mean or an Exponentially-Modified Gaussian gene mean
	#Option must be set to True in current version
	gammaDistMean = True

	#Gamma mean parameters
	shape = 0.35#0.429743565557383
	rate = 0.19

	#Exp modified mean parameters - not available in current version
	expLambda = 0#0.50
	geneScale = 0#0.47

	#Select cell size distribution
	cellRadiusDisp = 0.14

	###(Optional) Step 6 - Optional loading and graphing parameters###
	loadData = False
	showCellPlots = False
	showAlphaGeneMeans = False
	saveFolder = None #Set to None if you want to make it dynamic
	randomSeed = 393#123

	###Unused parameters###
	dropoutProportion = 0.0
	perplexityParam = 30


	###Start running commands###

	#Initialize variables
	cells = sum(cellTypeNumCells)
	cellTypesRecord = np.zeros(shape=cells, dtype="int64")
	totalCellIterations = np.zeros(shape=cells, dtype="int64")
	maxNumIterations = maxNumIterations + convDistAvgNum

	#Make folder for results
	if saveFolder == None:
		saveFolder = makeFolder(n, cells)

	#Set environment
	np.random.seed(randomSeed)
	random.seed(randomSeed)
	np.set_printoptions(threshold=20)

	#Determine important underlying objects
	n, gc, transcriptionFactors = generateGeneCorrelationMatrix(n, maxAlphaBeta, normalizeWRows, networkStructure, powerLawExponent, networkSparsity, saveFolder)
	geneMeanLevels = np.zeros(shape=n) #Unused in current version
	cellTypeTree, projMatrix, constMatrix = formCellTypeTree(n, cells, cellTypeParents, cellTypeNumCells, cellTypeConstProps, transcriptionFactors, cellTypeMeans)
	#geneMeanLevels = generateGeneMeanLevels(n, expScale)

	initialStates = generateInitialStates(n, cells, geneMeanLevels, sameInitialCell)
	print("Initial States: ", initialStates)
	#alpha determined based on optimal development for progenitor cell type (parent = -1)
	alpha = findOptimalAlpha(n, gc, timeStep, cellTypeTree, geneMeanLevels, convergenceThreshold, noiseDisp, initialStates, maxNumIterations, convDistAvgNum, cellsToSample, initialAlpha, alphaStep, showAlphaGeneMeans)
	gc = gc / alpha #Update gc with new alpha
	print("Alpha: ", alpha)

	if(loadData):
		gc = loadSavedData()
	print("Loaded Data: ", loadData)
	print("GC :", gc)

	#currentStates = developInitialStates(n, cells, initialStates, gc, timeStep, iterations, cellTypesRecord, totalCellIterations, projMatrix, constMatrix, noiseDisp, geneMeanLevels, showCellPlots, keepProportion, currentCellSplit, numCellSplits, cellDevelopmentMode)
	#print("Branch 0 Iterations: ", iterations)
	#print("Current State at Branch 0: \n", currentStates)
	#np.save("initialDevelopedStates_saved", currentStates)

	currentStates = initialStates
	#findAllNextStates recurses over cell type tree and develops cells
	findAllNextStates(gc, n, timeStep, cellTypeTree.head(), cellTypeTree, currentStates, cellTypesRecord, noiseDisp, convergenceThreshold, maxNumIterations, convDistAvgNum, cellsToSample, geneMeanLevels, showCellPlots, cellDevelopmentMode, totalCellIterations)
	finalCellStates = currentStates
	print("Final Cell States Before Transformations: \n", finalCellStates)

	cellSizeFactors = np.zeros(shape=cells)
	pseudotimes = timeStep * totalCellIterations
	randomGenes = np.random.randint(0, n, size=5)
	analyzeDataBeforeTransformation(n, cells, finalCellStates, saveFolder, alpha, randomGenes, pseudotimes, networkStructure)
	finalCellStates, cellSizeFactors = transformData(n, cells, finalCellStates, dropoutProportion, cellRadiusDisp, geneScale, expScale, expLambda, saveFolder, shape, rate, gammaDistMean)
	print("Final States After Transformations: \n", finalCellStates)

	saveData(finalCellStates, cellTypesRecord, gc, projMatrix, constMatrix, cellSizeFactors, pseudotimes, saveFolder, targetDirectory, R_Directory)
	analyzeSCRNAseqData(n, cells, perplexityParam, cellTypesRecord, finalCellStates, saveFolder)


#Method for later on when the R scripts can interface well the the Python Script

# def loadRInputs(R_Directory):

# 	r_object = 'metaInfo'
# 	robjects.r['load'](R_Directory + r_object + '.RData')
# 	metaInfo = np.array(robjects.r[r_object])

# 	n = int(metaInfo[0])
# 	cells = int(metaInfo[1])
# 	numRealClusters = int(metaInfo[2])

# 	r_object = 'clusterParents'
# 	robjects.r['load'](R_Directory + r_object + '.RData')
# 	clusterParents = np.array(robjects.r[r_object])

# 	r_object = 'projMatrix'
# 	robjects.r['load'](R_Directory + r_object + '.RData')
# 	projMatrix = np.array(robjects.r[r_object])

# 	r_object = 'constMatrix'
# 	robjects.r['load'](R_Directory + r_object + '.RData')
# 	constMatrix = np.array(robjects.r[r_object])

# 	r_object = 'clusterHeights'
# 	robjects.r['load'](R_Directory + r_object + '.RData')
# 	clusterHeights = np.array(robjects.r[r_object])

# 	r_object = 'clusterAssignments'
# 	robjects.r['load'](R_Directory + r_object + '.RData')
# 	clusterAssignments = np.array(robjects.r[r_object])

# 	r_object = 'clusterMeans'
# 	robjects.r['load'](R_Directory + r_object + '.RData')
# 	clusterMeans = np.array(robjects.r[r_object])

# 	cellTypeTree = formCellTypeTree(clusterParents, numRealClusters, clusterAssignments, projMatrix, constMatrix, clusterHeights, clusterMeans)
# 	maxCluster = int(max(clusterParents)) #Highest cluster number corresponds to the highest cluster
# 	for node in cellTypeTree.traverse(maxCluster):
# 		for field in node:
# 			print(field)

# 	return(n, cells, cellTypeTree)

#Generate the W matrix
def generateGeneCorrelationMatrix(n, maxAlphaBeta, normalizeWRows, networkStructure, powerLawExponent, networkSparsity, saveFolder):
	networkMatrix = np.zeros(shape=(n,n))

	#Generate the network structure
	if(networkStructure == "SimulatedPowerLaw"):
		#Generate number of in and out connections per gene
		inConnections = [int((1/2.0) * math.pow(1 - np.random.uniform(), -1/(powerLawExponent - 1)) + 1/(2.0)) for i in range(n)]
		outConnections = [int((1/2.0) * math.pow(1 - np.random.uniform(), -1/(powerLawExponent - 1)) + 1/(2.0)) for i in range(n)]
		
		#Make sure number of in and out connections is not greater than the number of genes
		for i in range(n):
			if(inConnections[i] > n):
				inConnections[i] = n
			if(outConnections[i] > n):
				outConnections[i] = n

		#Set up a multiset of the number of outconnections each gene has, then for each gene sample the inconnections from
		#the multiset of outconnections and remove duplicates 
		outNodesSamplePool = []
		for col in range(n):
			numOutConnection = outConnections[col]
			outNodesSamplePool += [col] * numOutConnection

		for row in range(n):
			numInConnections = inConnections[row]
			connectionIndexes = random.sample(outNodesSamplePool, numInConnections)
			connectionIndexes = list(set(connectionIndexes))
			networkMatrix[row, connectionIndexes] = 1

	elif(networkStructure == "TRRUST"):
		#This file must be downloaded, or it can be generated by using the final AnalyzeTRRUST.py and the original file
		networkMatrix = np.load("TRRUST_Network.npy")
		n = networkMatrix.shape[0]

	elif(networkStructure == "Random"):
		#Randomly generate network using sparsity and binomial distribution for the number of in connections per gene
		for row in range(n):
			numInConnections = np.random.binomial(n, (1-networkSparsity))
			connectionIndexes = random.sample(range(n), numInConnections)
			networkMatrix[row, connectionIndexes] = 1
		
	else:
		raise Exception('Invalid Network Structure')

	analyzeNetworkStructure(n, networkMatrix, networkStructure, saveFolder)

	valueMatrix = np.zeros(shape=(n,n))

	#Generate values for matrix from a Beta distribution using randomly generated alpha and beta in the range [0, 1, 2, ... , maxAlphaBeta]
	#Then transform the actual values from the (0, 1) range to the (-1, 1) range
	for col in range(n):
		alpha = int(maxAlphaBeta * np.random.uniform() + 1) 
		beta = int(maxAlphaBeta * np.random.uniform() + 1) 
		valueMatrix[col, :] = (2 * np.random.randint(0, 2) - 1) * np.random.beta(a=alpha, b=beta, size=n)

	#Element-wise multiply network structure matrix and value matrix
	gc = np.multiply(networkMatrix, valueMatrix)

	#If set to true, will divide each value of each row (in connections) by the number of connections
	if(normalizeWRows):
		for row in range(n):
			if(np.count_nonzero(gc[row, :]) != 0):
				gc[row, :] = gc[row, :] / np.count_nonzero(gc[row, :]) #Normalize by number of connections per row

	#Transcription Factors are defined as genes with a nonzero outdegree
	outDegrees = np.sum(networkMatrix, axis = 0)
	transcriptionFactors = np.nonzero(outDegrees)[0]

	return(n, gc, transcriptionFactors)

#Analysis of the indegree and outdegree distributions of the chosen network
def analyzeNetworkStructure(n, networkMatrix, networkStructure, saveFolder):
	os.chdir(saveFolder)

	inDegrees = np.sum(networkMatrix, axis = 1)
	outDegrees = np.sum(networkMatrix, axis = 0)

	#inDegreesXVals, inDegreesFreq = np.unique(inDegrees, return_counts=True)
	#outDegreesXVals, outDegreesFreq = np.unique(outDegrees, return_counts=True)

	plt.hist(inDegrees, bins = 100, range=[0,100])
	plt.xlabel("Indegree")
	plt.ylabel("Frequency")
	inDegreeString = networkStructure + "_Network_Indegree.tiff"
	plt.savefig(inDegreeString)
	plt.clf()
	plt.cla()

	plt.hist(outDegrees, bins = 100, range=[0,100])
	plt.xlabel("Outdegree")
	plt.ylabel("Frequency")
	outDegreeString = networkStructure + "_Network_Outdegree.tiff"
	plt.savefig(outDegreeString)
	plt.clf()
	plt.cla()

	print("Number of Connections: ", np.count_nonzero(networkMatrix.flatten()))
	sparsityProportion = float(np.count_nonzero(networkMatrix.flatten())) / (n**2)
	print("Network Density: ", sparsityProportion)

# def generateCellTypeDynamics(cellTypeParents, cellTypeConstProps, cellTypeNumCells, transcriptionFactors, cellTypeMeans):
# 	cellTypeTree = formCellTypeTree(cellTypeParents, cellTypeNumCells, cellTypeConstProps, transcriptionFactors, cellTypeMeans)
# 	return(n, cells, cellTypeTree)

def formCellTypeTree(n, cells, cellTypeParents, cellTypeNumCells, cellTypeConstProps, transcriptionFactors, cellTypeMeans):
	cellTypeTree = Tree()

	parentCellType = -1
	projInit = np.ones(shape=n)
	constInit = np.zeros(shape=n)
	cellTypeMember = []
	cellTypeChildMember = list(range(cells))
	cellTypeMean = 0
	numCellTypes = len(cellTypeParents)

	projMatrix = np.zeros(shape = (numCellTypes, n))
	constMatrix = np.zeros(shape = (numCellTypes, n))

	headNode = cellTypeTree.add_head(parentCellType, cellTypeMember, constInit, projInit, cellTypeMean)
	headNode.setChildCellTypeMembers(cellTypeChildMember)

	numCellTypes = len(cellTypeParents)
	cellTypeMembers = [[] for i in range(numCellTypes)]

	cellNumPool = range(cells)
	for cellType in range(numCellTypes):
		numCells = cellTypeNumCells[cellType]
		cellTypeIndices = random.sample(cellNumPool, numCells)
		cellTypeMembers[cellType] = cellTypeIndices
		cellNumPool = list(set(cellNumPool) - set(cellTypeIndices))

	addTreeChildren(n, cellTypeTree, cellTypeParents, cellTypeMembers, cellTypeConstProps, projInit, constInit, projMatrix, constMatrix, cellTypeMeans, parentCellType, transcriptionFactors)
	
	#maxCellType = int(max(cellTypeParents)) #Highest cell type number corresponds to the highest cell type
	for node in cellTypeTree.traverse(-1):
		for field in node:
			print(field)

	return cellTypeTree, projMatrix, constMatrix

def addTreeChildren(n, cellTypeTree, cellTypeParents, cellTypeMembers, cellTypeConstProps, projInit, constInit, projMatrix, constMatrix, cellTypeMeans, parentCellType, transcriptionFactors):
	childCellTypes = np.where(np.array(cellTypeParents) == parentCellType)[0]
	#print("Child Clusters: ", childClusters)
	parentCellTypeMember = []
	for i in range(len(childCellTypes)):
		childCellType = childCellTypes[i]

		constProp = cellTypeConstProps[childCellType]
		proj, const = getNextProjConst(n, projInit, constInit, constProp, transcriptionFactors)
		cellTypeMean = cellTypeMeans[childCellType]
		cellTypeMember = cellTypeMembers[childCellType]

		projMatrix[childCellType, :] = proj
		constMatrix[childCellType, :] = const

		cellTypeNode = cellTypeTree.add_node(childCellType, cellTypeMember, const, proj, cellTypeMean, parentCellType)

		childCellTypeMember = addTreeChildren(n, cellTypeTree, cellTypeParents, cellTypeMembers, cellTypeConstProps, proj, const, projMatrix, constMatrix, cellTypeMeans, childCellType, transcriptionFactors)
		parentCellTypeMember += childCellTypeMember
	cellTypeTree[parentCellType].setChildCellTypeMembers(parentCellTypeMember)

	return parentCellTypeMember + cellTypeMembers[parentCellType]

#Generate projection and constant vectors for next cell type
def getNextProjConst(n, projInit, constInit, constProp, transcriptionFactors):
	numZeros = int(len(transcriptionFactors)*constProp) #np.random.binomial(n, constProp, 1)[0]
	zeroGenes = random.sample(np.ndarray.tolist(transcriptionFactors), numZeros)

	proj = np.copy(projInit)
	proj[zeroGenes] = 0

	const = np.copy(constInit)
	for i in range(0,n):
		if (proj[i] == 0 and projInit[i] != 0): #Only if cell is zero in only new proj
			const[i] = rangeTransformation(random.randint(0, 1))
	return proj, const

#Dummy function
def generateGeneMeanLevels(n, expScale):
	geneMeanLevels = np.zeros(shape = (n))
	return(geneMeanLevels)

#Randomly generate an initial gene expression state, and use it for all the cells
def generateInitialStates(n, cells, geneMeanLevels, sameInitialCell):
	initialStates = np.zeros(shape = (n, cells))
	baseInitialState = rangeTransformation(np.random.rand(n)) + geneMeanLevels
	for j in range(cells):
		if(sameInitialCell):
			initialStates[:, j] = baseInitialState
		else:
			initialStates[:, j] = rangeTransformation(np.random.rand(n)) + geneMeanLevels
	return(initialStates)

#Transform cell states from (0, 1) range to (-1, 1) range
def rangeTransformation(cellStates):
	cellStates = 2*cellStates - 1 
	return (cellStates)

#Find the optimal alpha for W and the number of optimal iterations given that alpha
def findOptimalAlpha(n, gc, timeStep, cellTypeTree, geneMeanLevels, convergenceThreshold, noiseDisp, currentStates, maxNumIterations, convDistAvgNum, numToSample, initialAlpha, alpha_step, showAlphaGeneMeans):

	#Initialize variables (not to be changed)
	cells = currentStates.shape[1]
	successRateNeeded = 1
	numSuccessesNeeded = int(numToSample * successRateNeeded)
	successRate = 0
	if(initialAlpha >= alpha_step):
		alpha = initialAlpha - alpha_step
	else:
		alpha = 0
	avgIterations = 0
	normDiffVals = [0] * maxNumIterations

	treeHead = cellTypeTree.head()
	if(len(treeHead.children) == 1): #If only one cell lineage, use head
		origCellType = cellTypeTree[(treeHead.children)[0]]
		proj = origCellType.proj
		const = origCellType.const
	else: #Multiple cell lineages
		proj = np.zeros(shape=n)
		const = np.zeros(shape=n)

	#While proportion of cells that have converged is less than the required proportion
	while(successRate < successRateNeeded):
		geneMeans = np.zeros(shape=n)

		#Sample cells and save them in a new array so they don't overwrite old values
		cellsToSample = random.sample(range(cells), numToSample)
		simulatedStates = np.zeros(shape = (n, numToSample))
		for i in range(numToSample):
			simulatedStates[:, i] = currentStates[:, cellsToSample[i]]

		#Initialize these parameters at the start of every run of each alpha
		numSuccesses = 0
		if(successRate > 0): #Slow down alpha_step once at least one cell converges
			alpha_step = 0.01
		alpha = alpha + alpha_step
		avgIterations = 0

		#For each cell being sampled
		for i in range(numToSample):
			initialState = simulatedStates[:, i]
			cellIteration = 0
			normDiffVals = [0] * maxNumIterations
			convDist = float("inf")
			gc_alpha = gc / alpha

			#Develop each cell until the maximum number of iterations is reached or until it converges
			while(cellIteration < maxNumIterations and convDist > convergenceThreshold):
				currentState = developCell(n, gc_alpha, timeStep, initialState, proj, const, noiseDisp, geneMeanLevels)
				normDiffVals[cellIteration] = np.mean(np.absolute(currentState - initialState)/timeStep) #np.linalg.norm((currentState - initialState)) / n #vector norm
				
				#Average average value of derivative over previous "convDistAvgNum" iterations
				if(cellIteration >= (convDistAvgNum - 1)):
					convDist = np.mean(normDiffVals[(cellIteration - convDistAvgNum + 1):(cellIteration + 1)])
				initialState = currentState
				cellIteration += 1
			if(cellIteration != maxNumIterations): #If the cell has converged
				numSuccesses += 1
			avgIterations += cellIteration
			geneMeans += currentState
		geneMeans /= numToSample
		if(showAlphaGeneMeans):
			plt.hist(geneMeans)
			plt.title("Gene Means for alpha: " + str(alpha))
			plt.show()
		successRate = numSuccesses / numToSample
		avgIterations = int(avgIterations / numToSample)
		avgIterations = avgIterations - convDistAvgNum #Subtract iterations where it is below convergence distance
		print("Alpha: " , alpha)
		print("Success Rate: ", successRate)
		print("Average Iterations: ", avgIterations)
	#xVals = [x for x in range(0, cellIteration)]
	#plt.scatter(xVals, normDiffVals[0:cellIteration])
	#plt.plot(xVals, normDiffVals[0:cellIteration])
	#plt.show()
	return(alpha)

#Given a W matrix and the projection/constant vectors, find the minimum number of iterations until convergence
def findOptimalIterations(n, gc, timeStep, cellTypeMembers, proj, const, geneMeanLevels, convergenceThreshold, noiseDisp, currentStates, maxNumIterations, convDistAvgNum, numToSample):

	#Initialize variables (not to be changed)
	avgIterations = 0
	normDiffVals = [0] * maxNumIterations

	#Sample cells and save them in a new array so they don't overwrite old values
	colsToSample = random.sample(cellTypeMembers, numToSample)
	simulatedStates = np.zeros(shape = (n, numToSample))
	for i in range(numToSample):
		simulatedStates[:, i] = currentStates[:, colsToSample[i]]

	#For each cell being sampled
	for i in range(numToSample):
		initialState = simulatedStates[:, i]
		cellIteration = 0
		normDiffVals = [0] * maxNumIterations
		convDist = float("inf")

		#Develop each cell until the maximum number of iterations is reached or until it converges
		while(cellIteration < maxNumIterations and convDist > convergenceThreshold):
			currentState = developCell(n, gc, timeStep, initialState, proj, const, noiseDisp, geneMeanLevels)
			normDiffVals[cellIteration] = np.mean(np.absolute(currentState - initialState)/timeStep) #np.linalg.norm((currentState - initialState)) / n #vector norm
			
			#Average average value of derivative over previous "convDistAvgNum" iterations
			if(cellIteration >= (convDistAvgNum - 1)):
				convDist = np.mean(normDiffVals[(cellIteration - convDistAvgNum + 1):(cellIteration + 1)])
			initialState = currentState
			cellIteration += 1
		avgIterations += cellIteration
	avgIterations = int(avgIterations / numToSample)
	avgIterations = avgIterations - convDistAvgNum #Subtract iterations where it is below convergence distance
	print("Average Iterations: ", avgIterations)
	#xVals = [x for x in range(0, cellIteration)]
	#plt.scatter(xVals, normDiffVals[0:cellIteration])
	#plt.plot(xVals, normDiffVals[0:cellIteration])
	#plt.show()
	return(avgIterations)

#Equation to develop each cell for a time step
def developCell(n, gc, timeStep, initialState, proj, const, noiseDisp, geneMeanLevels):
	tau = np.ones(shape=n)
	nextState = np.tanh(gc.dot(initialState - geneMeanLevels) + np.random.normal(loc=0, scale=noiseDisp, size=(n)))
	derivative = (nextState - initialState) / tau 
	currentState = initialState + timeStep*derivative
	currentState = np.multiply(proj, currentState) + const + geneMeanLevels

	return currentState

#Develop the first "cell type" - the progenitor of all other "cell types"
# def developInitialStates(n, cells, initialStates, gc, timeStep, iterations, cellTypesRecord, totalCellIterations, projMatrix, constMatrix, noiseDisp, geneMeanLevels, showPlots, keepProportion, currentCellSplit, numCellSplits, cellDevelopmentMode):
# 	global currentCellType
# 	print("Current Cell Type: ", currentCellType)
# 	numKeepCellType = int(keepProportion*cells) #How many cells to keep at the first cell type 
# 	numDiffCellType = cells - numKeepCellType #How many cells to differentiate later
	
# 	cellTypesRecord[0:cells] = int(currentCellType) #Fill in cell types
# 	proj1 = projMatrix[currentCellType, :]
# 	const1 = constMatrix[currentCellType, :]
# 	currentCellType = currentCellType + 1

# 	for i in range(0, cells):
# 		initialState = initialStates[:, i]
# 		if(not cellDevelopmentMode): #If cells are supposed to all be fully developed
# 			cellIterations = np.random.poisson(iterations) #The number of iterations should be centered around the equilibrium value
# 		else: #Else if the cells are captured at different developmental time points
# 			if(i < numKeepCellType or currentCellSplit == numCellSplits): #If cell's final cell type is the current cell type
# 				cellIterations = np.random.randint(0, iterations + 1)
# 			else: #Else if the cell's cell type is not the current cell type
# 				cellIterations = iterations
# 		totalCellIterations[i] += cellIterations
# 		initialStates[:, i] = findCellNextState(gc, n, timeStep, i, proj1, const1, initialState, cellIterations, noiseDisp, geneMeanLevels, showPlots)

# 	return(initialStates)

def findCellNextState(gc, n, timeStep, i, proj, const, initialState, cellIterations, noiseDisp, geneMeanLevels, showPlots):
	normDiffVals = [0] * cellIterations
	currentState = initialState

	gene = int(np.random.uniform() * n)
	geneExpr = np.zeros(cellIterations)
	for x in range(0, cellIterations):
		yVals = initialState
		geneExpr[x] = initialState[gene]
		currentState = developCell(n, gc, timeStep, initialState, proj, const, noiseDisp, geneMeanLevels)
		normDiffVals[x] = np.mean(np.absolute(currentState - initialState)) #np.linalg.norm((currentState - initialState))
		initialState = currentState

	if(showPlots):
		xVals = [x for x in range(0, cellIterations)]

		plt.title("Convergence Graph")
		plt.plot(xVals, normDiffVals)
		plt.show()

		title = "Gene " + str(gene) + " Over Time"
		plt.scatter(xVals, geneExpr)
		plt.title(title)
		plt.plot(xVals, geneExpr)
		plt.show()

		plt.title("Current Distribution of Gene Expression")
		plt.hist(yVals, bins = int(n/25))
		plt.show()

	return(currentState)

def findAllNextStates(gc, n, timeStep, cellTypeNode, cellTypeTree, currentStates, cellTypesRecord, noiseDisp, convergenceThreshold, maxNumIterations, convDistAvgNum, numToSample, geneMeanLevels, showPlots, cellDevelopmentMode, totalCellIterations):
	currentCellType = cellTypeNode.identifier
	cellTypeMembers = cellTypeNode.cellTypeMembers
	proj = cellTypeNode.proj
	const = cellTypeNode.const
	childCellTypeMembers = cellTypeNode.childCellTypeMembers

	cellTypesRecord[cellTypeMembers] = int(currentCellType)
	print("Current Cell Type: ", currentCellType)

	if len(cellTypeMembers) == 0:
		cellTypeIterations = 0
	else:
		#print("Constant Gene Proportion: ", 1 - float(np.count_nonzero(proj))/np.size(proj))
		cellTypeIterations = findOptimalIterations(n, gc, timeStep, cellTypeMembers, proj, const, geneMeanLevels, convergenceThreshold, noiseDisp, currentStates, maxNumIterations, convDistAvgNum, numToSample)
		for i in range(len(cellTypeMembers)):
			if cellDevelopmentMode:
				cellIterations = np.random.randint(0, cellTypeIterations + 1)
			else:
				cellIterations = np.random.poisson(cellTypeIterations)
			cellTypeMember = cellTypeMembers[i]
			initialState = currentStates[:, cellTypeMember]
			currentStates[:, cellTypeMember] = findCellNextState(gc, n, timeStep, i, proj, const, initialState, cellIterations, noiseDisp, geneMeanLevels, showPlots)
			totalCellIterations[cellTypeMember] += cellIterations
		for i in range(len(childCellTypeMembers)):
			cellIterations = cellTypeIterations
			childCellTypeMember = childCellTypeMembers[i]
			initialState = currentStates[:, childCellTypeMember]
			currentStates[:, childCellTypeMember] = findCellNextState(gc, n, timeStep, i, proj, const, initialState, cellIterations, noiseDisp, geneMeanLevels, showPlots)
			totalCellIterations[childCellTypeMember] += cellIterations

	cellTypeNodeChildrenIDs = cellTypeNode.children
	for cellTypeNodeChildID in cellTypeNodeChildrenIDs:
		cellTypeNodeChild = cellTypeTree[cellTypeNodeChildID]
		findAllNextStates(gc, n, timeStep, cellTypeNodeChild, cellTypeTree, currentStates, cellTypesRecord, noiseDisp, convergenceThreshold, maxNumIterations, convDistAvgNum, numToSample, geneMeanLevels, showPlots, cellDevelopmentMode, totalCellIterations)

def analyzeDataBeforeTransformation(n, cells, finalCellStates, saveFolder, alpha, randomGenes, pseudotimes, networkStructure):
	os.chdir(saveFolder)

	#What do you want to name the graphs?
	testingAlpha = False
	testingNetworkStructure = True

	if(testingAlpha):
		testingString = "_Alpha_" + str(alpha) #What is being tested
	elif(testingNetworkStructure):
		testingString = "_" + networkStructure + "_Network" #What is being tested
	else:
		testingString = ""

	fig = plt.figure()
	SGeneMeans = np.mean(finalCellStates, axis = 1)
	geneMeansString = "S_Gene_Means" + testingString + ".tiff"
	plt.hist(SGeneMeans, bins = int(n/50))
	ax = fig.add_subplot(111)
	ax.set_xlabel("S Gene Means")
	ax.set_ylabel("Frequency")
	plt.savefig(geneMeansString)
	plt.clf()
	plt.cla()

	fig = plt.figure()
	SGeneVars = np.var(finalCellStates, axis = 1)
	geneVarsString = "S_Gene_Vars" + testingString + ".tiff"
	plt.hist(SGeneVars, bins = int(n/50))
	ax = fig.add_subplot(111)
	ax.set_xlabel("S Gene Variances")
	ax.set_ylabel("Frequency")
	plt.savefig(geneVarsString)
	plt.clf()
	plt.cla()

	fig = plt.figure()
	SCellMeans = np.mean(finalCellStates, axis = 0)
	cellMeansString = "S_Cell_Means" + testingString + ".tiff"
	plt.hist(SCellMeans, bins = int(cells/10))
	ax = fig.add_subplot(111)
	ax.set_xlabel("S Cell Means")
	ax.set_ylabel("Frequency")
	plt.savefig(cellMeansString)
	plt.clf()
	plt.cla()

	print("Random Genes Being Graphed: ", randomGenes)

	for randomGene in randomGenes:
		geneExpr = finalCellStates[randomGene, :]
		fig = plt.figure()
		genePseudotimeString = "Gene_" + str(randomGene) + "_Pseudotime" + testingString + ".tiff"
		plt.scatter(pseudotimes, geneExpr)
		ax = fig.add_subplot(111)
		ax.set_xlabel("Pseudotime")
		ax.set_ylabel("S Gene Expression")
		plt.savefig(genePseudotimeString)
		plt.clf()
		plt.cla()


def transformData(n, cells, finalCellStates, dropoutProportion, cellRadiusDisp, geneScale, expScale, expLambda, saveFolder, shape, rate, gammaDistMean):

	originalGeneMeans = np.mean(finalCellStates, axis = 1)
	geneMeanRankings = np.argsort(originalGeneMeans)

	if(gammaDistMean):
		newGeneMeans = np.log(np.random.gamma(shape, scale=1/rate, size=n))/expScale
		newGeneMeansSorted = np.sort(newGeneMeans)

		for i in range(len(geneMeanRankings)):
			geneIndex = geneMeanRankings[i]
			finalCellStates[geneIndex,:] += newGeneMeansSorted[i] - np.mean(finalCellStates[geneIndex,:])
	else:
		newGeneMeanOffsets = (geneScale + expLambda)/expScale + -1*np.random.exponential(expLambda, size=n) 
		newGeneMeansSorted = np.sort(newGeneMeanOffsets)

		for i in range(len(geneMeanRankings)):
			geneIndex = geneMeanRankings[i]
			finalCellStates[geneIndex,:] += newGeneMeansSorted[i]

	finalCellStates = np.exp(expScale*finalCellStates)

	cellSizeFactors = (np.random.normal(loc=1, scale=cellRadiusDisp, size=cells)**3)
	for j in range(cells):
		finalCellStates[:,j] *= cellSizeFactors[j]
	
	for i in range(n):
		for j in range(cells):
			if finalCellStates[i,j] < 0:
				finalCellStates[i,j] = 0
			finalCellStates[i,j] = np.random.poisson(lam=finalCellStates[i,j])

		#numCellKeep = numCellsKeep[i]
		#zeroCells = random.sample(range(cells), cells - numCellKeep)
		#finalCellStates[i, zeroCells] = 0

	return finalCellStates, cellSizeFactors

def plotSCRNAseqData(n, cells, cellStates, cellTypeString, saveFolder):
	numBinsScalar = 1
	newpath = saveFolder + "Cell_Type_" + cellTypeString
	if not os.path.exists(newpath):
		os.makedirs(newpath)

	os.chdir(newpath)

	geneExpr = cellStates[int(np.random.uniform() * n), :]
	plt.hist(geneExpr, bins = int(cells/numBinsScalar))
	plt.title("Expression of Random Gene")
	plt.savefig("Expression_Of_Gene")
	plt.clf()
	plt.cla()

	cellExpr = cellStates[:, int(np.random.uniform() * cells)]
	plt.hist(cellExpr, bins = int(n/numBinsScalar))
	plt.title("Expression of Random Cell")
	plt.savefig("Expression_Of_Cell")
	plt.clf()
	plt.cla()

	geneMeansString = "Gene_Means_" + "Cell_Type_" + cellTypeString
	librarySizesString = "Library_Sizes_" + "Cell_Type_" + cellTypeString
	numExpressedGenesString = "Num_Expressed_Genes_" + "Cell_Type_" + cellTypeString
	geneMeans = np.mean(cellStates, axis = 1)
	plt.hist(geneMeans, bins = int(n/numBinsScalar))
	plt.title("Gene Means") #Change back to Log[10]
	plt.savefig(geneMeansString)
	plt.clf()
	plt.cla()

	librarySizes = np.sum(cellStates, axis = 0)
	plt.hist(librarySizes, bins = int(n/numBinsScalar))
	plt.title("Library Sizes")
	plt.savefig(librarySizesString)
	plt.clf()
	plt.cla()

	numExpressedGenes = np.zeros(shape=int(cells))
	for cellNum in range(cells):
		numExpressedGenes[cellNum] = np.count_nonzero(cellStates[:, cellNum])
	plt.hist(numExpressedGenes, bins = int(n/numBinsScalar))
	plt.title("Num Expressed Genes")
	plt.savefig(numExpressedGenesString)
	plt.clf()
	plt.cla()
	
def analyzeSCRNAseqData(n, cells, perplexityParam, cellTypesRecord, finalCellStates, saveFolder):
	cellTypeString = "All"
	plotSCRNAseqData(n, cells, finalCellStates, cellTypeString, saveFolder)

	cellTypesSet = list(set(cellTypesRecord))
	for cellTypeNum in cellTypesSet:
		cellTypeString = str(cellTypeNum)
		cellTypeIndices = [i for i in range(cells) if cellTypesRecord[i] == cellTypeNum]
		finalCellStatesOfType = finalCellStates[:, cellTypeIndices]
		plotSCRNAseqData(n, len(cellTypeIndices), finalCellStatesOfType, cellTypeString, saveFolder)

	os.chdir(saveFolder)
	pca = PCA(n_components=15)
	pcaFit = pca.fit(finalCellStates)
	pcaResult = pca.components_
	print("PCA Explained Var: ", pca.explained_variance_ratio_) 
	np.save("Variance_Explained", np.asarray(pca.explained_variance_ratio_))
	xVals = pcaResult[0,:]
	yVals = pcaResult[1,:]
	plt.plot(xVals, yVals, 'bo')
	plt.savefig("PCA")

	# model = TSNE(n_components=2, perplexity=perplexityParam)
	# np.set_printoptions(suppress=True)
	# tsneResult = model.fit_transform(np.transpose(finalCellStates)) #Transpose because cells are data we are clustering

	# palette = np.array(sns.color_palette("hls", max(cellTypesRecord)+1))

	# # We create a scatter plot.
	# f = plt.figure(figsize=(8, 8))
	# ax = plt.subplot(aspect='equal')

	# colorHandles = []
	# for cellTypeNum in cellTypesSet:
	# 	colorHandles.append(mpatches.Patch(color=palette[cellTypeNum], label=cellTypeNum))

	# plt.legend(handles=colorHandles, borderaxespad=0)

	# sc = ax.scatter(tsneResult[:,0], tsneResult[:,1], lw=0, s=40, c=palette[cellTypesRecord.astype(np.int)])
	# ax.axis('off')
	# ax.axis('tight')
	# #plt.show()
	# plt.savefig("t-SNE")

def makeFolder(n, cells):
	cwd = os.getcwd()
	dataPath = cwd + "/Saved_Runs"
	experimentName = datetime.datetime.now().strftime("%m-%d-%Y+%H:%M:%S") + "_N:" + str(n) + "_C:" + str(cells)
	copyFolder = dataPath + "/" + experimentName + "/"
	if not os.path.exists(copyFolder):
		os.makedirs(copyFolder)

	try:
		import shutil #Doesn't work on MAC
		shutil.copy(__file__, copyFolder+"MatrixGenerator.py")
	except:
		print("Cannot copy file into folder")

	return (copyFolder)
	
def loadSavedData():
	gc = np.load("gc.npy")
	#projMatrix = np.load("projMatrix.npy")
	#constMatrix = np.load("constMatrix.npy")

	return gc#, projMatrix, constMatrix

def saveData(synthscRNAseq, cellTypesRecord, gc, projMatrix, constMatrix, cellSizeFactors, synthPseudotimes, copyFolder, targetDirectory, R_Directory):
	
	os.chdir(copyFolder)

	itemsToSave = [synthscRNAseq, cellTypesRecord, gc, projMatrix, constMatrix, cellSizeFactors, synthPseudotimes]
	itemNames = ["synthscRNAseq", "cellTypesRecord", "gc", "projMatrix", "constMatrix", "cellSizeFactors", "synthPseudotimes"]

	for i in range(len(itemsToSave)):
		np.save(itemNames[i], itemsToSave[i])
		np.save(targetDirectory+itemNames[i], itemsToSave[i])
		ro = numpy2ri(itemsToSave[i])
		r.assign(itemNames[i], ro)
		rSaveString = "save(" + itemNames[i] + ", file='"+ itemNames[i] +".gzip', compress=TRUE)"
		r(rSaveString)
		rSaveString = "save(" + itemNames[i] + ", file='" + R_Directory + ".gzip', compress=TRUE)"
		r(rSaveString)

main()