scRutiNy Manual
===============

Introduction
############
	
Single-cell RNA-seq data allows us to address cell subpopulations within a tissue sample and gene expression heterogeneity. However, the ability of scRNA-seq to reveal the gene expression of individual cells also makes it more difficult to analyze than bulk RNA-seq data, as cell development, differentiation and replication can confound the results obtained from differential expression analyses. scRNA-seq is prone to more noise than bulk RNA-seq measurement, which are averaged over many cells. We present here a program ``MatriSeq`` for simulating single-cell RNA-seq data as a dynamical system based on a gene-gene interaction network that models the development of each cell over time. It has the ability to generate different cell types through simulation of differential gene regulation in the context of the underlying dynamical system.

Generative models in general learn the joint probability distribution of the data and the model, while discriminative models learn conditional probability distributions. While discriminative models are better suited to classification tasks, generative models allow us to model the distributions of individual classes, allowing representations of the relations in the data. As single-cell data may have temporal non-stationarity and related soft boundaries between developing cells and differentiating cell types, generative models may be able to better extract biological function from this data.

We attempt inference of gene regulatory networks from scRNA-seq data in the context of our generative model through the program ``RegNetInference``. We study how well such networks can be reverse-engineered from the data, and quantify the increase in accuracy with increasing numbers of cells. We apply our method to real data and compare the performance of our model with published methodologies.

``scRutiNy`` is comprised of two modules: :class:`RNAscrutiny.MatriSeq` and :class:`RNAscrutiny.RegNetInference`. ``MatriSeq`` and ``RegNetInference`` are standalone programs; however, the input from ``MatriSeq`` can be read into ``RegNetInference``, or ``RegNetInference`` can read real data input and operate on its own.

Requirements and Setting Up
***************************

1. Necessary Libraries for ``MatriSeq`` and ``RegNetInference``

   1. ``MatriSeq`` libraries:
      
      - ``scipy``
      - ``scikit-learn``
      - ``rpy2``
      - ``numpy``
      - ``matplotlib``
      - ``tree`` (custom library)

   2. `RegNetInference` libraries:

      - ``scipy``
      - ``scikit-learn``
      - ``rpy2``
      - ``numpy``
      - ``matplotlib``
      - ``seaborn``
      - ``networkx``

2. Setting up the working directories

   1. Data and graphs produced by ``MatriSeq`` by default are output to the directory ``scRutiNy/MatriSeq`` but the directory can be specified differently by the user with the ``outputDir`` parameter passed to ``MatriSeq.init``.
   2. ``RegNetInference`` reads the data output from ``MatriSeq``; it expects the data to be in the directory ``scRutiNy/MatriSeq`` but this can be specified differently by the user the ``inputDir`` parameter passed to ``RegNetInference.init``. Data and graphs produced by ``RegNetInference`` by default are output to the directory ``scRutiNy/RegNetInference`` but the directory can be specified differently by the user with the ``outputDir`` parameter passed to ``RegNetInference.init``.

``MatriSeq``
############

1. Selecting the parameters for the gene coupling matrix

   - Selecting the gene coupling matrix size: Select the number of genes (``n``) to be simulated. The gene coupling matrix :math:`W` will describe how the gene expression of a cell at a time point affects the gene expression of the same cell at the next time point. The resulting :math:`W` matrix generated will be an :math:`n\times n` matrix. Default recommended value of ``n`` is 3,000.

   - Selecting the network structure: There are three options for ``networkStructure`` for the gene coupling matrix that can be selected:

     - ``SimulatedPowerLaw`` - Simulated power law network with ``powerLawExponent`` as parameter
     - ``TRRUST`` - Real network from TRRUST database of transcription factor interactions
     - ``Random`` - Random network with ``networkSparsity`` as a parameter

     The "Simulated Power Law" network structure is a network which follows a power law distribution in terms of the indegree and outdregree of the nodes. The ``powerLawExponent`` parameter controls what the exponent used for the power law distribution is -- in order for the distribution to be a power law, this exponent is usually between 2 and 3.

     The "TRRUST" network structure is a network based on literature-curated transcription factor interactions in the human regulatory network. The network is sourced from http://www.grnpedia.org/trrust/ and can be generated from the raw file on the website using a script included in the program. This regulatory network is good for simulating the realistic structure of a human regulatory network. It also follows a power law distribution and can be approximated with the right exponent used in the "Simulated Power Law" distribution. *Note:* The number of genes (``n``) will be set by ``MatriSeq`` if the TRRUST network is used, so the user input will be overridden.

     The "Random" network structure is a network in which connections are randomly assigned. The ``networkSparsity`` parameter controls what fraction of possible connections are actually set as connections. This regulatory network is good for studying how the network sparsity affects the outcomes of the simulation.

     Default recommended value of ``networkStructure`` is ``TRRUST``.

2. Controlling the magnitude of the gene coupling matrix

   - Selecting the minimum :math:`alpha`: The parameter ``alpha`` controls the size of the couplings of the gene coupling matrix :math:`W`. Originally, the values of the entries of :math:`W` range from :math:`(-1, 1)`, where a -1 represents a strong inhibitory regulatory relationship, while a 1 represents a strong activation relationship. However, these couplings can be scaled, because the original :math:`W` matrix is not guaranteed to cause the convergence of cells. The ``alpha`` parameter is a scaling parameter such that :math:`W` is multiplied by :math:`1/\alpha` in order to get the new value of ``W``. The ``initialAlpha`` parameter can be set by the user, and ``MatriSeq`` will test to see if the ``initialAlpha`` provided will create a ``W`` that causes convergence of cell types. If not, ``MatriSeq`` will increase ``alpha`` until it reaches a value that scales ``W`` to a value range that causes convergence of cell types. More on this process is described in "Matri-seq Step 5 - Selecting the parameters to control the smoothness of gene trajectories and the selection of alpha". 

   - Selecting the ``alpha`` to be lower will cause the gene trajectories to be less uniform and to take on more nonlinear trajectories. Default recommended value of ``initialAlpha`` is 1.

   - Selecting the maximum value for the :math:`\beta` distribution parameters: For each row :math:`i` of :math:`W`, which represents the incoming regulatory relationships that the gene :math:`i` has from other genes, a :math:`\beta` distribution is generated from an :math:`\alpha` and a :math:`\beta` parameter. The couplings for that gene are selected from the :math:`\beta` distribution for the gene. The :math:`alpha` and the :math:`beta` parameter for each gene are randomly selected from a uniform distribution (1, ``maxAlphaBeta``), so ``maxAlphaBeta`` can be selected to control the range of :math:`\alpha` and :math:`\beta` parameters to be used in the generation of gene coupling values. Default recommended value of ``maxAlphaBeta`` is 10

   - Selecting whether to normalize rows of W by the number of connections: The rows of ``W`` can be normalized by the number of connections each row has -- the couplings for row :math:`i` will all be divided by the number of connections that row :math:`i` has. This can be done to reduce the high variance that genes with a lot of regulators may demonstrate. Default recommended value of ``normalizeWRows`` is ``False``.

3. Selecting the cell development structure

   - Selecting whether to use the same progenitor cell for all cells: Each cell in `MatriSeq` is "developed" over pseudotime in order to produce the single-cell development effects in the data. As a result, each cell must start with some initial state. There are two options: all cells can start with the same initial state, or all cells can start with different initial states, and the boolean ``sameInitialCell`` controls this. If all cells start with the same initial state, this can be used to model a dataset where a stem cell divides and differentiates over time into many different cell types. If all cells start with different initial states, then this can be model a situation where all the cells have already differentiated or are differentiating along different, unrelated tracks. Default Recommended value of ``sameInitialCell`` is ``True``.

   - Selecting whether to capture cells during development or once development has finished: Cells during the development process undergo many changes, and some even progress through many progenitor cell types before they stop developing and stay as a certain cell type. Some single-cell experiments capture cells once they have all stopped developing and have reached a steady state. However, other single-cell experiments capture single cells at different times in the cells' development trajectories, so some cells have just begun developing, while others are almost finished. Programs like "Monocle" attempt to order these developing cells along their trajectories and assign them "pseudotimes" based on the estimated developmental time of capture. ``MatriSeq`` has the ability to simulate this aspect by stopping the cells after a certain number of iterations prior to the optimal number of developmental iterations.
     If the boolean ``developmentalMode`` is set to ``True``, then when a cell is developing into its final cell type, the number of iterations that it will develop for will be drawn from a uniform distribution (0, ``optimalIterations``). This allows for pseudotime inference algorithms to draw pseudotime values from the generated data, and ``MatriSeq`` will also provide a normalized pseudotime value for each cell based on the number of iterations so that this value can be compared against those inferred by pseudotime inference programs. 
     However, if the boolean ``developmentalMode`` is set to ``False``, then when a cell is developing into its final cell type, the number of iterations that it will develop for will be drawn from a poisson distribution with the mean equal to ``optimalIterations``. This allows there to be some variation in the number of iterations each cell develops for, but not as much variation as when ``developmentalMode`` is set to ``True``. Default recommended value of ``developmentalMode`` is ``False``.

   - Creating the cell type hierarchy: The cell type hierarchy is at the core of ``MatriSeq``. Cell types can have progenitors and descendants, and different cell types can have different attributes. Below we describe the different aspects of this cell type hierarchy:

     1) Selecting the number of cells per cell type: ``MatriSeq`` can support having different cell types in one simulation. The list ``cellTypeNumCells`` represents how many cells each cell type will end up with. For example, if the value of the entry with index 0 is 400 in ``cellTypeNumCells``, then that means that cell type 0 will end up with 400 cells, while if the value of the entry with index 2 is 200, then that means that cell type 2 will end up with 200 cells. Default Recommended value of ``cellTypeNumCells`` is ``[500]``.

     2) Selecting parent-child relationship between cell types: All cell types in ``MatriSeq`` must have a parent, even if that parent is the empty cell type represented by the index -1. The list ``cellTypeParents`` describes what the parent of each cell type is. Cell types are ordered and correspond to their index in the list -- for example, cell type 0 corresponds to the 0 index in the list, and cell type 1 corresponds to the 1 index in the list. If the value in ``cellTypeParents`` for index 1 is 0, then that means that cell type 0 is the parent of cell type 1. However, if the value for index 0 is -1, that means that the cell type 0 has no parent, so its parent is the empty cell type. Multiple cell types can have the empty cell type as their parent, which allows for independent cell type lineages.
        All cells will develop through their parent cell types first, starting with the parent with the empty cell type as its parent. However, the specified number of cells will end up as each cell type. Default recommended value of ``cellTypeParents`` is ``[-1]``.

     3) Selecting the parameters for the amount of constant genes per cell type: Each cell type in ``MatriSeq`` has a projection vector and a constant vector that cause differential gene expression between cell types. The projection vector sets certain genes to 0, and then the constant vector sets those genes to take on a value of 1 or -1, both of which have an equal probability of being selected. These genes are meant to represent input transcription factors, which are either being highly expressed (value 1) or lowly expressed (value -1) due to cell-type specific transcription factor expression. As genes affect each other in development, these expression values of input genes affect other genes, and this leads to cell type differentiation. Genes that are defined as transcription factors can only be genes that have a non-zero outdegree, so a list of genes with nonzero outdegrees is generated to be used.
        For each cell type, a proportion of the genes with nonzero outdegrees will be selected to be set as cell-type specific transcription factors, and this proportion can be selected for each cell type. The array cellTypeConstProps describes this proportion of genes for each cell type. Additionally, each cell type inherits the projection and constant vectors from its parent cell type, as children cell types generally share cell-type specific transcription factor expression with their parent cell type, and the resultant projection and constant vectors for a cell type is the union of the selected genes of the parent and cell-type specific projection and constant vectors. Default recommended value of ``cellTypeConstProps`` is ``[0.05]``.

     4) Selecting the cell type means: Different cell types can have different average levels of gene expression, as some cell types may have very active gene expression, while other cell types may have very inactive gene expression. ``MatriSeq`` will allow this feature of cell types to be replicated through selection of a gene expression mean per cell type; however, this feature is not available in the current version of Matri-seq. Unlike the projection/constant vectors, cell means are not inherited. Default recommended value of ``cellTypeMeans`` is ``None``.


4. Selecting the parameters to control the smoothness of gene trajectories and the selection of alpha and the optimal number of iterations
   Controlling gene trajectory smoothness:

   - Selecting how much to update each cell during each development step: In each iteration of development, each cell's state is updated as by the following equations:

     - ``nextState = tanh(W*(initialState - geneMeanLevels) + normal(0, noiseDisp, n)``

     - ``derivative = (nextState - initialState) / tau``

     - ``currentState = initialState + timeStep*derivative``

     - ``currentState = proj*currentState + const + geneMeanLevels``

     Basically, the "next" state is defined to be the state of the cell after the next application of the :math:`W` matrix to the cell state. This "next" state reflects a time step of unit 1 later. However, the gene trajectory can be made smoother by using this "next" step to calculate a derivative in the gene levels at the cell's current time point, and then to approximate the cell's state at a much closer time point. This is an application of the Euler method for approximating the value of a function at a nearby time point. The variable ``timeStep`` controls how far apart in time the cell's current state and its approximated next state are. The final update step applies the cell-type specific projection and constant vectors to the cell. If the timeStep variable is reduced, then the gene trajectory will be smoother, and if it is increased then the gene trajectory will more rough, and if it is set to 1 then the update will be equivalent to simply applying the :math:`W` matrix to the cell state at each iteration. *Note*: ``geneMeanLevels`` is subtracted out of the original computation and added back in at the end so that the application of the :math:`W` matrix can reflect changes in each gene from its baseline level. However, in the current version of the program ``geneMeanLevels`` is set to an array of 0's, since ``geneMeans`` are not able to be specified by the user yet in the current version. Default recommended value of timeStep is 0.01.

   - Selecting the amount of noise to add to each cell at each development step: In the above cell state update equations, the cell's "next" state is calculated by applying the :math:`W` matrix and also adding in a Gaussian random variable. This Gaussian random variable is a noise variable that represents how at each iteration, a cell's state fluctuates due to forces outside of gene regulation. The amount of fluctuation is determined by the variable ``noiseDisp``, which controls the standard deviation of the Gaussian random variable, which is centered at 0. The ``timeStep`` and ``noiseDisp`` are independent because the ``noiseDisp`` is applied to the "next" state, which is always a time step of 1 unit away. As a result, the ``noiseDisp`` should be set relative to the S value bounds of :math:`(-1, 1)`, as any noise fluctuations will be bounded such that a cell's state can not be greater than 1 or less than -1. Default recommended value of ``noiseDisp`` is 0.05.

     Controlling the selection of ``alpha`` and the optimal number of iterations: In ``MatriSeq``, a value of ``alpha`` is selected by starting with ``initialAlpha`` and testing to see if the :math:`W` matrix scaled by it allows for the '” number of iterations is determined, as different cell types have different projection and constant vectors, so some make take longer than others to converge. The process for determining the optimal number of iterations is the same as that for determining the optimal ``alpha``, except for this step the ``alpha`` is already fixed, but the same step that is run for each iteration of ``alpha`` is run.

     1) Selecting when to consider a cell as having converged: At each testing of the ``W`` produced by ``alpha``, a subset of cells are selected and developed over a certain number of iterations. "Convergence" is defined by when the average derivative of the change in a cell's state over a specified number of iterations is below a certain threshold ``convergenceThreshold``. As illustrated by the equations below:

        - ``derivative = (currentState – initialState)/timeStep``

        - ``if mean(derivative) < convergenceThreshold``, then convergence

        The convergence threshold can control how far apart the cell states of cells of the same cell type are from each other, as a higher ``convergenceThreshold`` will allow for cells to be more variable right before they have "converged", while a lower ``convergenceThreshold`` will allow for cells to be less variable right before they have "converged". Also, a lower ``convergenceThreshold`` should allow to a cell to converge in fewer iterations, so the program should run faster. Default recommended value for ``convergenceThreshold`` is 0.10.

     2) Selecting how large to make each ``alpha`` step: When an ``alpha`` value fails to produce a ``W`` matrix that causes convergence within a specified number of iterations, then the ``alpha`` value is increased by ``alphaStep``, and the new value of alpha is tested. Decreasing ``alphaStep`` will allow for a more accurate ``alpha`` to be found, but ``MatriSeq`` will take longer to run. *Note:* Do not decrease ``alphaStep`` below 0.01, as the program is set to dynamically adjust the ``alphaStep`` to 0.01 once the ``alpha`` is near the optimal ``alpha``. Default recommended value of ``alphaStep`` is 0.02.

     3) Selecting the maximum number of iterations until convergence: Potentially, any value of ``alpha`` could cause convergence if the program was given a long enough time to run. However, as time is limited, an upper bound on the number of iterations can be set through ``maxNumIterations``. Thus, if the cell does not converge with the current ``alpha`` after ``maxNumIterations`` iterations, then the ``alpha`` is incremented by ``alphaStep``. If the user wants time series data, then ``maxNumIterations`` can be increased to give more time points, since the cells should undergo more iterations before development since the ``maxNumIterations`` has been increased, so a smaller ``alpha`` can be found that produces a ``W`` that causes cell convergence. Default recommended value of ``maxNumIterations`` is 500.

     4) Selecting how many cells to sample: In order to determine convergence, ``cellsToSample`` cells are randomly selected from the pool of all cells, and convergence tests are performed on them. In order for "convergence" to have occurred, all cells in the sample will need to have "converged". As a result, increasing ``cellsToSample`` will lead to a more accurate alpha and optimal number of iterations being found, but it will also cause Matri-seq to take longer. Default recommended value of ``cellsToSample`` is 50.

     5) Selecting how many iterations to average over to determine convergence: When a cell has converged, the mean of its derivative over ``convDistAvgNum`` occurred iterations must be below ``convergenceThreshold``. Not only must the current derivative be below ``convergenceThreshold``, but also the previous ``convDistAvgNum`` iterations. This is done so that the derivative does not fluctuate to a low value and trigger "convergence", but rather the cell must be consistently stable to have converged. Default recommended value of ``convDistAvgNum`` is 20.

5. Selecting the parameters to create realistic characteristics of scRNA-seq data: The data outputted by the first four steps is in the range :math:`(-1, 1)`, and the distribution of average gene expression per gene and the distribution of average gene expression per cell are both centered around 0. As scRNA-seq data does not exhibit these characteristics, since it cannot be negative and average gene expression is not 0, ``MatriSeq`` must apply some transformations to make the data more realistic. Many of the transformations are adapted from ``Splatter``.

   - Selecting the base of the exponent used to transform the data: The original ``S`` matrix created by the first three steps contains values in the range :math:`(-1, 1)`, where -1 represents a gene with very low expression, and a 1 represents a gene with very high expression. However, real gene expression is of the range :math:`(0, \inf)`. This aspect can be simulated by exponentiating the data, meaning that each value in the matrix is replaced by an exponent base to the power of the previous value. This exponent base can be selected as a power of e -- the parameter ``expScale`` controls the exponent base so that the base is :math:`\exp` (``expscale``). Now, the new range of values that  the data can take is (:math:`\exp` (-``expScale``), :math:`e` (``expScale``)). The exponent scale controls the variance of gene expression, but it does not control the mean of gene expression, since the gene means are overridden by the gene mean distribution selected in the next step. Default recommended value of ``expScale`` is 4.

   - Selecting the distribution of gene means:

     1) Selecting whether to use a Gamma or Exponentially-Modified Gaussian distributed gene mean: The average expression per gene can be modeled using a Gamma distribution, as is the practice in Splatter. Currently, only a Gamma distributed mean can be used, as the exponentially-modified Gaussian distribution will be available in future versions, so ``gammaDistMean`` must be set to ``True``. Both distributions give the distribution of the logarithm of the gene means a left skew. The gene means are scaled to match that of the selected distribution. Default recommended value of ``gammaDistMean`` is ``True``.

     2) Selecting the Gamma distributed gene mean shape and rate: The Gamma distribution is governed by the shape and rate parameters. Both of these parameters can be selected to make the distribution of gene means match that of a real dataset. The shape parameter controls the skewness of the distribution, while the rate can be used to adjust the mean. Default recommended value of ``shape`` is 0.35 and ``rate`` is 0.20.

     3) Selecting the Exponentially-Modified Gaussian distributed gene mean lambda and gene scale: This option is not currently available in the current version.

   - Selecting the dispersion of the cell radius: After the output of the previous two steps, the distribution of the average gene expression per cell is narrow, because on average each cell should have a similar mean gene expression. However, real scRNA-seq data has a wide distribution of average gene expression per cell, and this is because cells have different sizes. The expression values represented mRNA concentrations, which means that they are independent of the cell size. However, in this step the concentrations are scaled based on the cell size. It is assumed that the distribution of cell radius is normally distributed, with a mean of 1 and a standard deviation of ``cellRadiusDisp``, and each cell is assigned a radius.The mean of the cell radius distribution is 1 so that the mean gene expression should not be significantly affected across all genes. Each cell's volume is then calculated by cubing the cell radius, and then the mRNA count is found by multiplying the mRNA concentration by the cell volume. Increasing the ``cellRadiusDisp`` naturally increases the width of the distribution of the average gene expression per cell. Default recommended value of ``cellRadiusDisp`` is 0.15.

``RegNetInference``
###################

1. Selecting the maximum pseudotime distance between cell pairs

   - Select the maximum pseudotime distance: ``RegNetInference`` finds the gene correlation matrix by analyzing the differences in gene expression between one cell at one time point and another cell at another time point. The program then tries to find the best :math:`W` (gene correlation) matrix that explains the changes in gene expression over time. In ``RegNetInference``, the input pseudotimes must be scaled so that the minimum pseudotime is 0 and the maximum pseudotime is 1.
     If ``RegNetInference`` were to compute the :math:`W` matrix that fits best for all possible pairs of cells, the number of cell pairs would be too large, and the error between a cell pair with a pseudotime near 0 and another cell pair with a pseudotime near 1 would be so large that it would be difficult to make the :math:`W` matrix able to successfully model such a large change in pseudotime. As a result, the ``pseudotimeThresh`` parameter can be set to determine the maximum pseudotime difference between two cells in order for that pair to be considered in the calculation of :math:`W`. Default recommended value of ``pseudotimeThresh`` is 0.10.

2. Select the sliding window size used for averaging over cells

   - Select the sliding window size: The cells are ordered by pseudotime in ``RegNetInference`` so that for each gene, a trajectory over pseudotime can be found and plotted. However, this trajectory is often very rough, as there is a lot of error in scRNA-seq experiments, and the Poisson sampling of reads makes values discrete. As a result, ``RegNetInference`` provides a smoothing function that sets each cell's expression profile to the average of the previous and following cells. The parameter ``windowSize`` controls how many cells are averaged over to produce a cell's state value -- for each cell roughly ``windowSize``/2 cells before it and ``windowSize``/2 cells after it (in pseudotime) are averaged into the cell's state. Default recommended value of ``windowSize`` is 30.

