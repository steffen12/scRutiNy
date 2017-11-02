import matplotlib.pyplot as plt
import numpy as np
import os
import random
from scipy.stats import rankdata
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from sklearn.linear_model import Lasso
import seaborn as sns
import matplotlib.patches as mpatches
import math
from collections import Counter
import networkx as nx
from sklearn import metrics

# from scipy.stats import spearmanr
#from sklearn.cluster import AgglomerativeClustering
#from sklearn.cluster import SpectralClustering
#from sklearn.linear_model import LinearRegression

__inputDir = None
__outputDir = None
__verbose = None


def init(inputDir='./scRutiNy/MatriSeq',
         outputDir='./scRutiNy/RegNetInference',
         seed=1337,
         verbose=True):

    global __inputDir, __outputDir, __verbose

    # set directories
    __inputDir = inputDir
    __outputDir = outputDir
    try:
        os.makedirs(outputDir)
    except:
        pass

    # seed random number generators
    np.random.seed(seed)
    random.seed(seed)

    # set verbosity
    __verbose = verbose


def readCellStates():

    cellStates = np.load(os.path.join(__inputDir, "synthscRNAseq.npy"))
    pseudotimes = np.load(os.path.join(__inputDir, "synthPseudotimes.npy"))
    cellTypesRecord = np.load(os.path.join(__inputDir, "cellTypesRecord.npy"))
    W_real = np.load(os.path.join(__inputDir, "gc.npy"))

    n, cells = cellStates.shape

    cellStates = convertToS(n, cells, cellStates)

    if __verbose:
        print("Num genes: ", n, ", num cells: ", cells)
        print("Cell States: ", cellStates)
        print("Num TF Connections: ", np.flatnonzero(W_real).size)
        print("Pseudotimes: ", pseudotimes)

    return n, cells, cellStates, pseudotimes, cellTypesRecord, W_real


def convertToS(n, cells, cellStates):
    for cellNum in range(cells):
        cellRankings = rankdata(cellStates[:, cellNum])
        cellStates[:, cellNum] = (2 * ((cellRankings + 1) / (n)) - 1)

    for geneNum in range(n):
        geneRankings = rankdata(cellStates[geneNum, :])
        cellStates[geneNum, :] = (2 * ((geneRankings + 1) / (cells)) - 1)

    return cellStates


def smoothCellData(cells, cellStates, pseudotimes, windowSize=30):
    #Please choose only an odd window size
    cellStatesSmooth = np.copy(cellStates)
    pseudotimeOrderedIndexes = np.argsort(pseudotimes)
    for j in range(len(pseudotimeOrderedIndexes)):
        pseudotimeWindowIndexes = []
        pseudotimeOriginalIndex = pseudotimeOrderedIndexes[j]
        windowStart = j - int(windowSize / 2)
        if (windowStart < 0):
            windowStart = 0
        windowEnd = j + 1 + int(windowSize / 2)  #Add one because not inclusive
        if (windowEnd > cells):
            windowEnd = cells
        for k in range(windowStart, windowEnd):
            pseudotimeWindowIndexes.append(pseudotimeOrderedIndexes[k])
        cellStatesSmooth[:, pseudotimeOriginalIndex] = np.mean(
            cellStates[:, pseudotimeWindowIndexes], axis=1)
    return cellStatesSmooth


def clusterGenes(n,
                 cells,
                 cellStates,
                 perplexity=30,
                 epsilon=0.05,
                 min_samples=2):
    #numClusters = 100

    #numNeighbors = 5

    distanceMatrix = getCorrelationMatrix(n, cells, cellStates)
    if __verbose:
        print("Distance Matrix: ", distanceMatrix)
        #print(distanceMatrix)

    #trace = go.Heatmap(z=distanceMatrix)
    #data = [trace]
    #heatmapFile = "heatmap_" + datetime.datetime.now().strftime("%m-%d-%Y+%H:%M:%S")

    #plotly.offline.plot(data, filename=heatmapFile)
    #py.iplot(data, filename=heatmapFile)

    model = TSNE(n_components=2, perplexity=perplexity)
    np.set_printoptions(suppress=True)
    tsneResult = model.fit_transform(
        cellStates)  #Transpose because cells are data we are clustering

    dbClust = DBSCAN(
        eps=epsilon, min_samples=min_samples,
        metric="precomputed").fit(distanceMatrix)
    labels = dbClust.labels_  #Labels of -1 are noisy

    #dbClust = DBSCAN(eps=epsilon, min_samples=min_samples, metric="precomputed").fit(distanceMatrix)
    #labels = dbClust.labels_ #Labels of -1 are noisy

    #aggClust = AgglomerativeClustering(n_clusters=numClusters, affinity="precomputed", linkage="average").fit(distanceMatrix)
    #labels = aggClust.labels_

    #specClust =  SpectralClustering(n_clusters=numClusters, affinity="precomputed", n_neighbors=numNeighbors).fit(distanceMatrix)
    #labels = specClust.labels_

    if __verbose:
        print("Gene Clusters: ", labels)

    newClustNumber = max(labels) + 1
    if min(labels) == -1:  #There are noise labels
        for i in range(len(labels)):
            if (labels[i] == -1):
                labels[i] = newClustNumber
                newClustNumber += 1

    labelsSet = list(set(labels))
    palette = np.array(sns.color_palette("hls", len(labelsSet)))

    if __outputDir is not None:
        f = plt.figure(figsize=(8, 8))
        ax = plt.subplot(aspect='equal')

        colorHandles = []
        for cellTypeNum in labelsSet:
            colorHandles.append(
                mpatches.Patch(color=palette[cellTypeNum], label=cellTypeNum))

        # plt.legend(handles=colorHandles, borderaxespad=0)

        sc = ax.scatter(
            tsneResult[:, 0],
            tsneResult[:, 1],
            lw=0,
            s=40,
            c=palette[labels.astype(np.int)])
        ax.axis('off')
        ax.axis('tight')
        plt.savefig(os.path.join(__outputDir, 'cluster.tiff'))
        plt.close()

    return labels


def getCorrelationMatrix(n, cells, cellStates):
    cellStatesNorm = np.copy(cellStates)

    geneMeans = np.mean(cellStatesNorm, axis=1)
    for cell in range(cells):
        cellStatesNorm[:, cell] -= geneMeans

    geneStds = np.std(cellStatesNorm, axis=1)
    for gene in range(n):
        if (geneStds[gene] != 0):
            cellStatesNorm[gene, :] /= geneStds[gene]
            #plt.hist(cellStatesNorm[gene, :])
            #plt.show()
        else:
            cellStatesNorm[gene, :] = np.zeros(shape=cells)

    correlationMatrix = np.dot(cellStatesNorm,
                               np.transpose(cellStatesNorm)) / (cells - 1)

    distanceMatrix = (np.ones((n, n)) - correlationMatrix) / 2

    return distanceMatrix


def plotGenesOverPseudotime(n, pseudotimes, cellStates):
    sortedIndexes = [
        i[0]
        for i in sorted(enumerate(pseudotimes.tolist()), key=lambda x: x[1])
    ]
    if __outputDir is not None:
        for gene in range(n):
            fig = plt.figure()
            ax = plt.gca()
            plt.title("Gene " + str(gene) + " over Pseudotime")
            ax.scatter(pseudotimes[sortedIndexes],
                       cellStates[gene, sortedIndexes])
            plt.savefig(
                os.path.join(__outputDir, 'Gene_' + str(gene) +
                             '_over_pseudotime.tiff'))
            plt.close()


def crossValidateData(n,
                      cells,
                      cellStates,
                      pseudotimes,
                      cellTypesRecord,
                      W_real,
                      geneClusters,
                      numPartitions=5,
                      alphaMin=1e-4,
                      alphaMax=1e-3,
                      alphaMultStep=0.01,
                      pseudotimeThresh=0.1):
    meanSquaredErrorThreshold = 3

    minMeanSquaredErrorArray = np.zeros(numPartitions)
    alphaArray = np.zeros(numPartitions)
    sensitivityArray = np.zeros(numPartitions)
    specificityArray = np.zeros(numPartitions)
    precisionArray = np.zeros(numPartitions)
    W_array = np.zeros(shape=(n, n, numPartitions))

    totalMinMeanSquaredError = float("inf")
    totalMinMeanSquaredErrorIndex = -1

    if __verbose:
        print("Cluster Membership: ", Counter(geneClusters))

    cellIndices = [x for x in range(cells)]
    random.shuffle(cellIndices)
    cellsPerPartition = int(cells / numPartitions)
    for partition in range(numPartitions):
        cvStart = partition * cellsPerPartition
        cvEnd = (partition + 1) * cellsPerPartition
        if (partition + 1) == numPartitions:
            cvEnd = cells
        cvIndices = cellIndices[cvStart:cvEnd]
        dataIndices = list(set(cellIndices) - set(cvIndices))

        dataCellStates = cellStates[:, dataIndices]
        dataPseudotimes = pseudotimes[dataIndices]
        dataCellTypesRecord = cellTypesRecord[dataIndices]

        cvCellStates = cellStates[:, cvIndices]
        cvPseudotimes = pseudotimes[cvIndices]
        cvCellTypesRecord = cellTypesRecord[cvIndices]

        alpha = alphaMax

        minMeanSquaredError = float("inf")
        #minMeanSquaredErrorDiff = float("inf")
        optAlpha = -1
        optW = -1
        optW_random = -1
        optSensitivity = -1
        optSpecificity = -1
        optPrecision = -1
        meanSquaredError = 0

        while (alpha > alphaMin and meanSquaredError <
               meanSquaredErrorThreshold * minMeanSquaredError):
            if __verbose:
                print("Analysis: ")
                print("Alpha: ", alpha)
            W = findWMatrix(n, cells, W_real, dataCellStates,
                            dataCellTypesRecord, pseudotimeThresh,
                            dataPseudotimes, alpha)
            if __verbose:
                print("Total Magnitude of W predicted: ", np.linalg.norm(W))
                print("Total Magnitude of W_real: ", np.linalg.norm(W_real))

            sensitivity, specificity, precision = compareToRealW(
                n, W, W_real, geneClusters)
            meanSquaredError = getCVCost(pseudotimeThresh, W, cvCellStates,
                                         cvPseudotimes, cvCellTypesRecord)
            if __verbose:
                print("Mean Squared Error: ", meanSquaredError)
                print("Comparison to Random W:")

            W_random = generateRandomW(W)
            compareToRealW(n, W_random, W_real, geneClusters)
            randomMeanSquaredError = getCVCost(pseudotimeThresh, W_random,
                                               cvCellStates, cvPseudotimes,
                                               cvCellTypesRecord)
            if __verbose:
                print("Random W Mean Squared Error: ", randomMeanSquaredError)

            if (meanSquaredError < minMeanSquaredError):
                minMeanSquaredError = meanSquaredError
                #minMeanSquaredErrorDiff = (randomMeanSquaredError - meanSquaredError)
                optAlpha = alpha
                optSensitivity = sensitivity
                optSpecificity = specificity
                optPrecision = precision
                optW = W
                optW_random = W_random

            alpha *= alphaMultStep

        if __verbose:
            print("Optimal Alpha: ", optAlpha)
            print("Optimal MSE: ", minMeanSquaredError)
        minMeanSquaredErrorArray[partition] = minMeanSquaredError
        alphaArray[partition] = optAlpha
        sensitivityArray[partition] = optSensitivity
        specificityArray[partition] = optSpecificity
        precisionArray[partition] = optPrecision
        W_array[:, :, partition] = optW

        if (minMeanSquaredError < totalMinMeanSquaredError):
            minMeanSquaredError = totalMinMeanSquaredError
            totalMinMeanSquaredErrorIndex = partition

        #W_cluster = getWClusterPred(n, optW, geneClusters)
        #W_real_cluster = getWClusterReal(n, W_real, geneClusters)
        #W_random_cluster = getWClusterPred(n, optW_random, geneClusters)

        plotROCCurve(n, optW, W_real, geneClusters)
        plotROCCurve(n, optW_random, W_real, geneClusters)

    if __verbose:
        print("Summary Statistics: ")
        print("Alpha Mean and Sd: ", np.mean(alphaArray), ", ",
              np.std(alphaArray))
        print("MSE Mean and Sd: ", np.mean(minMeanSquaredErrorArray), ", ",
              np.std(minMeanSquaredErrorArray))
        print("Sensitivity Mean and Sd: ", np.mean(sensitivityArray), ", ",
              np.std(sensitivityArray))
        print("Specificity Mean and Sd: ", np.mean(specificityArray), ", ",
              np.std(specificityArray))
        print("Precision Mean and Sd: ", np.mean(precisionArray), ", ",
              np.std(precisionArray))
        print("W Sd: ", np.mean(np.std(W_array, axis=2)))

    finalW = W_array[:, :, totalMinMeanSquaredErrorIndex]

    return finalW


def findWMatrix(n, cells, W_real, cellStatesNew, cellTypesRecordNew,
                pseudotimeThresh, pseudotimesNew, alpha):
    #T = np.load("T.npy")#
    #T = np.zeros(shape=(n,n))#np.random.rand(n,n)#
    #M = np.load("M.npy")#
    #M = np.zeros(shape=(n,n))#np.random.rand(n,n)#
    W = np.zeros(shape=(n, n))
    cellStatesNew = np.reshape(cellStatesNew, cellStatesNew.shape + (1, ))
    pseudotimesNew = np.reshape(pseudotimesNew, pseudotimesNew.shape + (1, ))
    cellTypesRecordNew = np.reshape(cellTypesRecordNew,
                                    cellTypesRecordNew.shape + (1, ))

    # try:
    # 	raise Exception('Starting Anew')
    # 	cellStatesTotal = np.load("cellStatesInfSaved.npy")
    # 	pseudotimesTotal = np.load("pseudotimesInfSaved.npy")
    # 	cellTypesRecordTotal = np.load("cellTypesRecordTotal.npy")
    # 	cellStatesTotal = np.concatenate((cellStatesTotal, cellStatesNew), axis=2)
    # 	pseudotimesTotal = np.concatenate((pseudotimesTotal, pseudotimesNew), axis=1)
    # 	cellTypesRecordTotal = np.concatenate((cellTypesRecordTotal, cellTypesRecordNew), axis=1)
    # except:
    # 	print("Exception!")
    cellStatesTotal = cellStatesNew
    pseudotimesTotal = pseudotimesNew
    cellTypesRecordTotal = cellTypesRecordNew

    numDatasets = cellStatesTotal.shape[2]
    totalNumCellPairs = 0

    if __outputDir is not None:
        np.save(
            os.path.join(__outputDir, "cellStatesInfSaved"), cellStatesTotal)
        np.save(
            os.path.join(__outputDir, "pseudotimesInfSaved"), pseudotimesTotal)

    #print("Finding Cell Pairs")

    cellPairsTotal = [[] for i in range(numDatasets)]
    for datasetIndex in range(numDatasets):
        cellStates = cellStatesTotal[:, :, datasetIndex]
        pseudotimes = pseudotimesTotal[:, datasetIndex]
        cellTypesRecord = cellTypesRecordTotal[:, datasetIndex]
        cellPairs = [
        ]  #Find these using pseudotime - can be made more efficient
        for cellType in range(max(cellTypesRecord) + 1):
            cellsOfType = np.where(cellTypesRecord == cellType)
            for cell1 in np.nditer(cellsOfType):
                pseudotime1 = pseudotimes[cell1]
                for cell2 in np.nditer(cellsOfType):
                    pseudotime2 = pseudotimes[cell2]
                    if (math.fabs(pseudotime1 - pseudotime2) <=
                            pseudotimeThresh and pseudotime1 > pseudotime2):
                        if (np.any(cellStates[:, cell1] - cellStates[:, cell2]
                                   != 0)):
                            cellPairs.append([cell1, cell2])
        cellPairsTotal[datasetIndex] = cellPairs
        totalNumCellPairs += len(cellPairs)
    if __outputDir is not None:
        np.save(os.path.join(__outputDir, "cellPairsTotal"), cellPairsTotal)
    if __verbose:
        print("Total Num Cell Pairs: ", totalNumCellPairs)

    #print("Finding Tau")
    tau = np.zeros(shape=(n))
    #negRangeBound = np.amin(cellStates, axis = 1)
    #posRangeBound = np.amax(cellStates, axis = 1)
    #print("Negative Range Bound:", negRangeBound)
    #print("Positive Range Bound:", posRangeBound)
    # for gene in range(n):
    # 	tauEstCells = np.zeros(shape=len(cellPairs))
    # 	for cellPairIndex in range(len(cellPairs)):
    # 		cellPair = cellPairs[cellPairIndex]
    # 		cellIndex1 = cellPair[0]
    # 		cellIndex2 = cellPair[1]
    # 		Sc1 = cellStates[:, cellIndex1]
    # 		Sc2 = cellStates[:, cellIndex2]
    # 		derivative = (Sc1[gene] - Sc2[gene])/(pseudotimes[cellIndex1] - pseudotimes[cellIndex2])
    # 		if(derivative != 0):
    # 			if(derivative > 0):
    # 				tauEstCells[cellPairIndex] = (1 - Sc2[gene])/derivative
    # 			elif(derivative < 0):
    # 				tauEstCells[cellPairIndex] = (-1 - Sc2[gene])/derivative
    # 		else:
    # 			tauEstCells[cellPairIndex] = float("inf")
    # 	print(scp.mquantiles(tauEstCells))
    # 	tau[gene] = np.min(tauEstCells)
    # print(tau)
    tau = pseudotimeThresh * np.ones(shape=(
        n))  #Doing this right now because don't need tau
    if __outputDir is not None:
        np.save(os.path.join(__outputDir, "tau"), tau)

    #print("Finding W")

    expProdsTotal = [[] for i in range(numDatasets)]
    expProdMeansTotal = [[] for i in range(numDatasets)]
    for datasetIndex in range(numDatasets):
        cellPairs = cellPairsTotal[datasetIndex]
        cellStates = cellStatesTotal[:, :, datasetIndex]
        pseudotimes = pseudotimesTotal[:, datasetIndex]
        expProds = np.zeros(shape=(n, len(cellPairs)))
        #print("Num Cell Pairs: ", len(cellPairs))
        for cellPairIndex in range(len(cellPairs)):
            cellPair = cellPairs[cellPairIndex]
            cellIndex1 = cellPair[0]
            cellIndex2 = cellPair[1]
            Sc1 = cellStates[:, cellIndex1]
            Sc2 = cellStates[:, cellIndex2]
            derivative = (Sc1 - Sc2) / (
                pseudotimes[cellIndex1] - pseudotimes[cellIndex2])
            #print("Pseudotime Differences: ", pseudotimes[cellIndex1] - pseudotimes[cellIndex2] + pseudotimeOffset)
            #print("Derivative: ", derivative)
            #print("Pre Arctanh: ", np.multiply(tau, derivative) + Sc2)
            product = np.multiply(tau, derivative) + Sc2
            while (np.max(np.absolute(product)) >= 1):
                product = 0.999 * product
            expProds[:, cellPairIndex] = np.arctanh(product)
            #print("ExpProd: ", np.arctanh(np.multiply(tau, derivative) + Sc2))
            #print("ActualProd: ", np.dot(W_real, Sc2))
        expProdsTotal[datasetIndex] = expProds
        expProdMeansTotal[datasetIndex] = np.nanmean(expProds, axis=1)

    geneMeansTotal = np.mean(cellStatesTotal, axis=1)

    W = np.zeros(shape=(n, n))
    for gene in range(n):
        data = np.zeros(shape=(totalNumCellPairs, n))
        y = np.zeros(shape=(totalNumCellPairs))
        cellPairIndexOffset = 0
        for datasetIndex in range(numDatasets):
            cellPairs = cellPairsTotal[datasetIndex]
            cellStates = cellStatesTotal[:, :, datasetIndex]
            expProdMeans = expProdMeansTotal[datasetIndex]
            geneMeans = geneMeansTotal[:, datasetIndex]
            expProds = expProdsTotal[datasetIndex]
            for cellPairIndex in range(len(cellPairs)):
                cellPair = cellPairs[cellPairIndex]
                cellIndex1 = cellPair[0]
                cellIndex2 = cellPair[1]
                Sc1 = cellStates[:, cellIndex1]
                Sc2 = cellStates[:, cellIndex2]
                expProd = expProds[:, cellPairIndex]
                y[cellPairIndexOffset
                  + cellPairIndex] = expProd[gene] - expProdMeans[gene]
                data[cellPairIndexOffset + cellPairIndex, :] = Sc2 - geneMeans
        cellPairIndexOffset += len(cellPairs)
        reg_result = lasso_regression(
            data, y,
            alpha)  #neural_network(n, data, y, alpha)#linear_reg(data, y)#
        #print("Gene ", gene, " RSS: ", reg_result[0])
        #print("Expected Cost: ", np.sum((y - np.dot(data, np.transpose(W_real[gene, :])))**2))
        #print("Base Cost:", np.sum(y**2))
        #print("W[gene]: ", reg_result[2:])
        #print("W_real[gene]: ", W_real[gene, :])
        W[gene, :] = reg_result[2:]
        #print("W[gene] connections: ", np.nonzero(W[gene, :]))
        #print("W_real[gene] connections: ", np.nonzero(W_real[gene, :]))

    #print("Predicted W Matrix : ", W)
    #print("Real W Matrix: ", W_real)

    # totalCost = 0
    # totalExpCost = 0
    # for cellPairIndex in range(len(cellPairs)):
    # 	cellPair = cellPairs[cellPairIndex]
    # 	cellIndex1 = cellPair[0]
    # 	cellIndex2 = cellPair[1]
    # 	Sc1 = cellStates[:, cellIndex1]
    # 	Sc2 = cellStates[:, cellIndex2]
    # 	expProd = expProds[:, cellPairIndex]
    # 	predictedSc1 = np.dot(W, Sc2)
    # 	expSc1 = np.dot(W_real, Sc2)
    # 	totalCost = totalCost + np.sum((expProd - predictedSc1)**2)
    # 	totalExpCost = totalExpCost + np.sum((expProd - expSc1)**2)
    # print("Total Cost: ", totalCost)
    # print("Total Expected Cost: ", totalExpCost)

    if __outputDir is not None:
        np.save(os.path.join(__outputDir, "W"), W)
    return (W)


def lasso_regression(data, y_vals, alpha, models_to_plot={}):
    #Fit the model
    lassoreg = Lasso(alpha=alpha, max_iter=1e5, fit_intercept=False)
    #print("Data: ", data)
    #print("y_vals: ", y_vals)
    lassoreg.fit(data, y_vals)
    y_pred = lassoreg.predict(data)

    #Check if a plot is to be made for the entered alpha
    if alpha in models_to_plot and __outputDir is not None:
        fig = plt.figure()
        ax = plt.gca()
        ax.plot(data, y_pred)
        ax.plot(data, y_vals, '.')
        plt.title('Plot for alpha: %.3g' % alpha)
        plt.savefig(os.path.join(__outputDir, 'alpha_%.3g' % (alpha, )))
        plt.close()

    #Return the result in pre-defined format
    rss = sum((y_pred - y_vals)**2)
    ret = [rss]
    ret.extend([lassoreg.intercept_])
    ret.extend(lassoreg.coef_)
    return ret


def compareToRealW(n, W, W_real, geneClusters):
    geneClustersSet = list(set(geneClusters))
    W_diff_matrix = W - W_real
    W_diff = np.linalg.norm(W_diff_matrix)

    #print("Analysis: ")
    #print("Difference in W's: ", W_diff)
    #print("Proportion of Difference: ", W_diff/(np.absolute(np.linalg.norm(W_real) + np.linalg.norm(W))))
    #print("Total Magnitude of W predicted: ", np.linalg.norm(W))
    #print("Total Magnitude of W_real: ", np.linalg.norm(W_real))

    #plt.title("W Predicted vs W Real")
    #plt.scatter(W.flatten(), W_real.flatten())
    #plt.show()

    #print("Spearman rho:", spearmanr(W.flatten(), W_real.flatten())[0])

    #print("Cluster Membership: ", Counter(geneClusters))

    sensitivity = -1
    specificity = -1
    precision = -1

    truePositives = 0
    trueNegatives = 0
    falsePositives = 0
    falseNegatives = 0

    W_cluster = getWClusterPred(n, W, geneClusters)
    W_real_cluster = getWClusterReal(n, W_real, geneClusters)

    for i in range(n):
        for clusterNum in geneClustersSet:
            if (W_cluster[i, clusterNum] == 0):
                predictedVal = 0
            else:
                predictedVal = 1
            trueVal = W_real_cluster[i, clusterNum]

            if (predictedVal == 0 and trueVal == 0):
                trueNegatives += 1
            elif (predictedVal == 1 and trueVal == 1):
                truePositives += 1
                print(i, clusterNum)
            elif (predictedVal == 0 and trueVal == 1):
                falseNegatives += 1
            elif (predictedVal == 1 and trueVal == 0):
                falsePositives += 1

    if __verbose:
        print("TP: ", truePositives, ". FP: ", falsePositives, ". TN: ",
              trueNegatives, ". FN: ", falseNegatives, ".")

    try:
        sensitivity = float(truePositives) / (truePositives + falseNegatives)
        specificity = float(trueNegatives) / (trueNegatives + falsePositives)
        precision = float(truePositives) / (truePositives + falsePositives)
        negativePredictiveValue = float(trueNegatives) / (
            trueNegatives + falseNegatives)
        accuracy = float(truePositives + trueNegatives) / (
            truePositives + falseNegatives + trueNegatives + falsePositives)

        if __verbose:
            print("Sensitivity/Recall: ", sensitivity, ", Specificity: ",
                  specificity)
            print("Precision: ", precision, ", Negative Predictive Value: ",
                  negativePredictiveValue)
            print("Accuracy: ", accuracy)
    except:
        print("Error in calculating statistics")

    return sensitivity, specificity, precision


def plotROCCurve(n, W, W_real, geneClusters):
    #W_flattened = np.abs(W.flatten())
    #W_real_flattened = np.abs(W_real.flatten())
    #W_real_flattened_binary = np.ones(shape=W_real_flattened.shape)
    #W_real_flattened_binary[np.where(W_real_flattened == 0)] = 0
    #W_real_flattened_binary = W_real_flattened
    tpr, fpr = roc_metrics(n, geneClusters, W_real, W)
    AUC = metrics.auc(fpr, tpr)
    if __verbose:
        print("ROC AUC: ", AUC)
        print("FPR: ", fpr)
        print("TPR: ", tpr)
    if __outputDir is not None:
        fig = plt.figure()
        ax = plt.gca()
        ax.plot(fpr, tpr)
        plt.title("ROC Curve")
        plt.savefig(os.path.join(__outputDir, 'ROC_Curve.tiff'))
        plt.close()


def roc_metrics(n, geneClusters, W_real, W):
    geneClustersSet = list(set(geneClusters))
    tprList = [1]
    fprList = [1]

    sensitivity = -1
    specificity = -1
    precision = -1

    truePositives = 0
    trueNegatives = 0
    falsePositives = 0
    falseNegatives = 0

    W_flattened = np.abs(getWClusterPred(n, W, geneClusters)).flatten()
    W_real_flattened_binary = np.abs(
        getWClusterReal(n, W_real, geneClusters)).flatten()

    steps = 10
    stepSize = (max(W_flattened) - min(W_flattened)) / steps
    for step in range(steps + 1):
        currentThreshold = step * stepSize
        for i in range(len(W_flattened)):
            if (W_flattened[i] > currentThreshold):
                predictedVal = 1
            else:
                predictedVal = 0
            trueVal = W_real_flattened_binary[i]
            if (predictedVal == 0 and trueVal == 0):
                trueNegatives += 1
            elif (predictedVal == 1 and trueVal == 1):
                truePositives += 1
            elif (predictedVal == 0 and trueVal == 1):
                falseNegatives += 1
            elif (predictedVal == 1 and trueVal == 0):
                falsePositives += 1

        print("TP: ", truePositives, ". FP: ", falsePositives, ". TN: ",
              trueNegatives, ". FN: ", falseNegatives, ".")

        try:
            sensitivity = float(truePositives) / (
                truePositives + falseNegatives)
            specificity = float(trueNegatives) / (
                trueNegatives + falsePositives)
            precision = float(truePositives) / (truePositives + falsePositives)
            negativePredictiveValue = float(trueNegatives) / (
                trueNegatives + falseNegatives)
            accuracy = float(truePositives + trueNegatives) / (
                truePositives + falseNegatives + trueNegatives + falsePositives
            )

            #print("Sensitivity/Recall: ", sensitivity, ", Specificity: ", specificity)
            #print("Precision: ", precision, ", Negative Predictive Value: ", negativePredictiveValue)
            #print("Accuracy: ", accuracy)

            tprList.append(sensitivity)
            fprList.append(1 - specificity)
        except:
            print("Error in calculating statistics")

    tprList.append(0)
    fprList.append(0)

    tprList = list(reversed(tprList))
    fprList = list(reversed(fprList))

    return (tprList, fprList)


def saveW(W):
    if __outputDir is not None:
        np.save(os.path.join(__outputDir, "W_Output"), W)


def getWClusterPred(n, W, geneClusters):
    geneClustersSet = list(set(geneClusters))

    W_cluster = np.zeros(shape=(n, len(geneClustersSet)))
    for i in range(n):
        for clusterNum in geneClustersSet:
            clusterPredictedRegSum = 0
            for j in range(n):
                if (geneClusters[j] == clusterNum and W[i, j] != 0):
                    clusterPredictedRegSum += np.abs(W[i, j])
            W_cluster[i, clusterNum] = clusterPredictedRegSum
    return (W_cluster)


def getWClusterReal(n, W_real, geneClusters):
    geneClustersSet = list(set(geneClusters))

    W_real_cluster = np.zeros(shape=(n, len(geneClustersSet)))
    for i in range(n):
        for clusterNum in geneClustersSet:
            for j in range(n):
                if (geneClusters[j] == clusterNum and W_real[i, j] != 0):
                    W_real_cluster[i, clusterNum] = 1
    return (W_real_cluster)


def getCVCost(pseudotimeThresh, W, cvCellStates, cvPseudotimes,
              cvCellTypesRecord):
    totalNumCellPairs = 0
    meanSquaredError = 0
    tau = pseudotimeThresh
    cellPairs = []  #Find these using pseudotime - can be made more efficient

    for cellType in range(max(cvCellTypesRecord) + 1):
        cellsOfType = np.where(cvCellTypesRecord == cellType)
        for cell1 in np.nditer(cellsOfType):
            pseudotime1 = cvPseudotimes[cell1]
            for cell2 in np.nditer(cellsOfType):
                pseudotime2 = cvPseudotimes[cell2]
                if (math.fabs(pseudotime1 - pseudotime2) <= pseudotimeThresh
                        and pseudotime1 > pseudotime2):
                    if (np.any(cvCellStates[:, cell1] - cvCellStates[:, cell2]
                               != 0)):
                        cellPairs.append([cell1, cell2])
    totalNumCellPairs += len(cellPairs)

    for cellPairIndex in range(len(cellPairs)):
        cellPair = cellPairs[cellPairIndex]
        cellIndex1 = cellPair[0]
        cellIndex2 = cellPair[1]
        Sc1 = cvCellStates[:, cellIndex1]
        Sc2 = cvCellStates[:, cellIndex2]
        derivative = (Sc1 - Sc2) / (
            cvPseudotimes[cellIndex1] - cvPseudotimes[cellIndex2])
        actualValue = derivative * tau + Sc2
        predictedValue = np.tanh(np.dot(W, Sc2))
        errorSquared = np.mean((actualValue - predictedValue)**2)
        meanSquaredError += errorSquared
    meanSquaredError /= len(cellPairs)
    return (meanSquaredError)


def generateRandomW(W):
    W_random = np.copy(W)
    np.random.shuffle(W_random)
    W_random = np.transpose(W)
    np.random.shuffle(W_random)
    W_random = np.transpose(W)
    return (W_random)


def plotWGraphs(n, W, W_real):

    names = ['Inferred', 'Real']

    for iW, Wi in enumerate([W, W_real]):
        geneGraph = nx.DiGraph()
        for i in range(n):
            for j in range(n):
                if (Wi[i, j] != 0):
                    geneGraph.add_edge(j, i)
        nx.draw_random(geneGraph, with_labels=True)
        plt.savefig(os.path.join(__outputDir, 'W_Graph_%s' % (names[iW], )))
        plt.close()
