Matri-seq Example
=================

Import the `RNAscrutiny` package and initialize with `init` function. Here, we set the output directory to a new directory with the current timestamp, and we set the verbosity of the package to `False`.

.. ipython::

   In [0]: import RNAscrutiny.MatriSeq as ms

   In [0]: import datetime

   In [0]: ms.init(outputDir=datetime.datetime.now().strftime('%d_%m_%Y-%H:%M:%S'), verbose=False)

Use the `generateCorrlationMatrix` function to generate a gene correlation matrix. The `networkStructure` can be one of `SimulatedPowerLaw`, `TRRUST`, or `Random`. `n` denotes the number of genes.

.. ipython::

   In [0]: networkStructure = 'SimulatedPowerLaw'

   In [0]: n = 700

   In [0]: gc, tf = ms.generateGeneCorrelationMatrix(n=n, networkStructure=networkStructure)


Represent the cell development as a tree with the `formCellTypeTree` function.

.. ipython::

   In [0]: cells = 200

   In [0]: cellTypeTree, projMatrix, constMatrix = ms.formCellTypeTree(n=n, cells=cells, cellTypeParents=[-1], cellTypeNumCells=[200], cellTypeConstProps=[0.0], transcriptionFactors=tf, cellTypeMeans=[0])


Determine the minmum value of alpha that allows convergence of the cell states.

.. ipython::

   In [0]: initialStates = ms.generateInitialStates(n, cells)

   In [0]: alpha = ms.findOptimalAlpha(n, gc, cellTypeTree, initialStates, numToSample=50)

   In [0]: gc /= alpha


Recurse over cell type tree to develop cells.

.. ipython::

   In [0]: currentStates = initialStates

   In [0]: cellTypesRecord = np.zeros(shape=cells, dtype="int64")

   In [0]: totalCellIterations = np.zeros(shape=cells, dtype="int64")

   In [0]: ms.findAllNextStates(gc, cellTypeTree.head(), cellTypeTree, currentStates, cellTypesRecord, totalCellIterations)

   In [0]: finalCellStates = currentStates


Output plots to the output directory set above.

.. ipython::

   In [0]: timeStep = 0.01

   In [0]: pseudotimes = timeStep * totalCellIterations

   In [0]: ms.analyzeDataBeforeTransformation(n, cells, finalCellStates, alpha, pseudotimes, networkStructure=networkStructure)


Get the final cell states

.. ipython::

   In [0]: finalCellStates, cellSizeFactors = ms.transformData(n, cells, finalCellStates)


Output plots to the output directory set above:

.. ipython::

   In [0]: ms.analyzeSCRNAseqData(n, cells, cellTypesRecord, finalCellStates)
