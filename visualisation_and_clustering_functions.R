# #cran
# install.packages(c("ggplot2", "Matrix", "dbscan", "pheatmap", "fpc", 
#"dynamicTreeCut", "factoextra", "digest", "RColorBrewer", "doParallel"))
# #bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("BiocParallel", "scran", "scater", "monocle", 
#"SingleCellExperiment", "KEGGREST"))

suppressMessages(library(BiocParallel, warn.conflicts = F))
suppressMessages(library(scran, warn.conflicts = F))
suppressMessages(library(ggplot2, warn.conflicts = F))
suppressMessages(library(scater, warn.conflicts = F))
suppressMessages(library(Matrix, warn.conflicts = F))
suppressMessages(library(monocle, warn.conflicts = F))
suppressMessages(library(SingleCellExperiment, warn.conflicts = F))
suppressMessages(library(dbscan, warn.conflicts = F))
suppressMessages(library(pheatmap, warn.conflicts = F))
suppressMessages(library(fpc, warn.conflicts = F))
suppressMessages(library(dynamicTreeCut, warn.conflicts = F))
suppressMessages(library(factoextra, warn.conflicts = F))
suppressMessages(library(digest, warn.conflicts = F))
suppressMessages(library(KEGGREST, warn.conflicts = F))
suppressMessages(library(RColorBrewer, warn.conflicts = F))
suppressMessages(require(doParallel))
suppressMessages(require(matrixStats))
suppressMessages(library(AnnotationDbi, warn.conflicts = F))
suppressMessages(library(dplyr, warn.conflicts = F))
suppressMessages(library(biomaRt, warn.conflicts = F))
suppressMessages(library(org.Mm.eg.db, warn.conflicts = F))

checkTSNE <- function(dataDirectory, experimentName, PCs=c(4, 6, 8, 10, 20, 40, 50), 
                      perplexities=c(30, 40)){
  
  SAMP <- experimentName
  
  eq <- c()
  
  for(i in 1:14){
    
    coordinatesName <- paste0(experimentName, '_tsne_coordinates_', i, "_",
                              PCs[i %% length(PCs)], "PCs_",
                              perplexities[((i-1) %/% length(PCs)) + 1], "perp")
  
    
    TSNEres_1 <- read.delim(file.path(paste0(dataDirectory, "/tsnes"), 
                                      paste0(coordinatesName, ".tsv")),
                            stringsAsFactors = FALSE)
    
    
    
    TSNEres_2 <- read.delim(file.path(paste0(dataDirectory, "/tsnes_first_run"), 
                                    paste0(SAMP,"_tsne_coordinates_",i,".tsv")),
                            stringsAsFactors = FALSE)
    
    eq <- c(eq, all.equal(TSNEres_1, TSNEres_2))
    
  }
  
  return(all(eq == TRUE))
}

### this function calculates PCA and then tSNE with PCs and perplexities ###
### it returns a list of pSNE = PCA+tSNE results ###
### to get XY coordinates call psne_res[1,i][[1]] ###
### i = [1:length(PCs)*length(perplexities)] is a number of iteration ###

getTSNEresults <- function(expressionMatrix, cores=1,
                           PCs=c(4, 6, 8, 10, 20, 40, 50), 
                           perplexities=c(30, 40), randomSeed=42){
    PCAData <- prcomp(t(expressionMatrix))$x
    myCluster <- makeCluster(cores, # number of cores to use
                             type = "PSOCK") # type of cluster
    registerDoParallel(myCluster)
    tSNECoordinates <- foreach(PCA=rep(PCs, length(perplexities)), 
                               perp=rep(perplexities, each=length(PCs)), 
                               .combine='cbind') %dopar% {
                                   library(SingleCellExperiment)
                        scater::plotTSNE(SingleCellExperiment(assays=list(
                            logcounts=t(PCAData[,1:PCA]))),
                        scale_features=FALSE, perplexity=perp, 
                        rand_seed=randomSeed, theme_size=13, return_SCESet=FALSE)
                        }
    stopCluster(myCluster)
    message(paste("Calculated", length(PCs)*length(perplexities), 
"2D-tSNE plots."))
    return(tSNECoordinates)
}

testClustering <- function(sceObject, dataDirectory, experimentName, 
                           dbscanEpsilon=1.4,
                           minPts=5,
                           perplexities = c(30), PCs = c(4),
                           randomSeed = 42,
                           width=7, height=7, onefile=FALSE, #pdf 
                           family, title, fonts, version,
                           paper, encoding, bg, fg, pointsize, 
                           pagecentre, colormodel,
                           useDingbats, useKerning, 
                           fillOddEven, compress){
  
  initialisePath(dataDirectory)
  dir.create(file.path(dataDirectory, "test_clustering"), showWarnings = F)
  
  message("Generating TSNE.")
  #1. Generating 2D tSNE plots
  tSNE <- getTSNEresults(expr = exprs(sceObject), cores=1, 
                         perplexities = perplexities, 
                         PCs = PCs,
                         randomSeed = randomSeed)
  
  message("Saving results.")
  #picture checker
  pdf(file.path(dataDirectory, "test_clustering", "test_tSNE.pdf"),
      width=width, height=height, onefile=onefile, # not changed by default
      family=family, title=title, fonts=fonts, version=version,
      paper=paper, encoding=encoding, bg=bg, fg=fg, pointsize=pointsize, 
      pagecentre=pagecentre, colormodel=colormodel,
      useDingbats=useDingbats, useKerning=useKerning, 
      fillOddEven=fillOddEven, compress=compress)
  print(tSNE)
  dev.off()
  
  
  #2. Clustering with dbscan
  # choosing for the best epsilon
  pdf(file.path(dataDirectory, "test_clustering", "distance_graph.pdf"),
      width=width, height=height, onefile=onefile, # not changed by default
      family=family, title=title, fonts=fonts, version=version,
      paper=paper, encoding=encoding, bg=bg, fg=fg, pointsize=pointsize, 
      pagecentre=pagecentre, colormodel=colormodel,
      useDingbats=useDingbats, useKerning=useKerning, 
      fillOddEven=fillOddEven, compress=compress)
  plotDistanceGraphWithEpsilon(tSNE$data, epsilon=dbscanEpsilon, 
                               minNeighbours = minPts)
  dev.off()
  
  pdf(file.path(dataDirectory, "test_clustering", "test_clustering.pdf"),
      width=width, height=height, onefile=onefile, # not changed by default
      family=family, title=title, fonts=fonts, version=version,
      paper=paper, encoding=encoding, bg=bg, fg=fg, pointsize=pointsize, 
      pagecentre=pagecentre, colormodel=colormodel,
      useDingbats=useDingbats, useKerning=useKerning, 
      fillOddEven=fillOddEven, compress=compress)
  print(plotTestClustering(tSNE$data, epsilon=dbscanEpsilon,
                           minNeighbours = minPts))
  dev.off()
  
  message("Pictures of test clustering were exported.")
  return(list(tSNE, 
              plotDistanceGraphWithEpsilon(tSNE$data, epsilon=dbscanEpsilon, 
                                                 minNeighbours = minPts),
              plotTestClustering(tSNE$data, epsilon=dbscanEpsilon,
                                 minNeighbours = minPts)))
  
}

choosePalette <- function(colorPalette, clustersNumber){
  # chooses the palette with the number of color fitted to the
  # number of clusters. Depending on number of clusters there are different default,
  # but you always can pick your own.
  
  colorPalette26 <- c( "yellow", "darkgoldenrod1", "coral1", "deeppink", 
                       "indianred", "coral4", "darkblue", "darkmagenta", 
                       "darkcyan", "mediumorchid", "plum2", "gray73", 
                       "cadetblue3", "khaki",
                       "darkolivegreen1", "lightgreen", "limegreen", 
                       "darkolivegreen4", "green", "#CC79A7", "violetred3", 
                       "brown3", "darkslategray1", "gray51", "slateblue2", 
                       "blue")
  
  pickDefaultPalette <- function(clustersNumber, colorPalette26){
    if(clustersNumber < 13) return(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", 
                                     "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                     "#CAB2D6", "#6A3D9A", "#FFFF99", 
                                     "#B15928")[1:clustersNumber])
    return(rep(colorPalette26,
              round(clustersNumber/length(colorPalette26))+1)[1:clustersNumber])
  }
  
  if(length(colorPalette) == 1){
    if(colorPalette == "default"){
      return(pickDefaultPalette(clustersNumber, colorPalette26))
    }
  }
  
  if(clustersNumber > length(colorPalette)){
    message("The number of clusters is greater than the number of colors. 
Using default CONCLUS palette instead.")
    return(pickDefaultPalette(clustersNumber, colorPalette26))
  }
  
  return(colorPalette[1:clustersNumber])
}

chooseStatePalette <- function(statesNumber){

  colorPalette <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
                      "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
                      "#CCEBC5", "#FFED6F")
  return(rep(colorPalette,
             round(statesNumber/length(colorPalette)) + 1)[1:statesNumber])
}

runCONCLUS <- function(sceObject, dataDirectory, experimentName, 
                       colorPalette="default",
                       statePalette="default",
                       clusteringMethod="ward.D2",
                       epsilon=c(1.3, 1.4, 1.5), minPoints=c(3, 4), k=0, 
                       PCs=c(4, 6, 8, 10, 20, 40, 50),
                       perplexities=c(30,40), 
                       randomSeed = 42,
                       deepSplit=4, preClustered = F,
                       orderClusters = FALSE,
                       cores=14,
                       plotPDFcellSim = TRUE,
                       deleteOutliers = TRUE){
  
  initialisePath(dataDirectory)
  
  # Generating 2D tSNE plots
  tSNEResults <- generateTSNECoordinates(sceObject, dataDirectory, 
                                          experimentName, PCs=PCs, 
                                          perplexities=perplexities,
                                          randomSeed = randomSeed)

  if(preClustered){
    # Running dbscan
    message("Running dbscan using ", cores, " cores.")
    dbscanResults <- runDBSCAN(tSNEResults, sceObject, dataDirectory, 
                               experimentName, epsilon=epsilon, 
                               minPoints=minPoints,
                               cores=cores)
    
    # assigning cells to clusters
    message("Calculating cells similarity matrix.")
    cellsSimilarityMatrix <- clusterCells(dbscanResults, sceObject, clusterNumber=k, 
                                      deepSplit=deepSplit, cores=cores,
                                      clusteringMethod=clusteringMethod)[[2]]
    sceObjectFiltered <- sceObject
  } else {
    # Running clustering
    clusteringResults <- runClustering(tSNEResults, sceObject, dataDirectory, 
                                       experimentName,
                                       epsilon=epsilon, minPoints=minPoints, 
                                       k=k, deepSplit=deepSplit,
                                       cores=cores,
                                       clusteringMethod=clusteringMethod,
                                       deleteOutliers = deleteOutliers,
                                       PCs=PCs,
                                       perplexities=perplexities, 
                                       randomSeed = randomSeed)
    sceObjectFiltered <- clusteringResults[[1]]
    cellsSimilarityMatrix <- clusteringResults[[2]]
  }
  
  print(table(colData(sceObjectFiltered)$clusters, 
              dnn=list("Cells distribuion by clusters")))
  
  clustersNumber <- length(unique(colData(sceObjectFiltered)$clusters))
  colorPalette <- choosePalette(colorPalette, clustersNumber)

  # Plotting cluster stablility and 2D visualisations
  plotCellSimilarity(sceObjectFiltered, cellsSimilarityMatrix, dataDirectory,
                     experimentName, colorPalette, 
                     orderClusters = orderClusters, 
                     statePalette = statePalette, 
                     clusteringMethod = clusteringMethod,
                     plotPDF = plotPDFcellSim)
  
  plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
                    PCs=PCs, perplexities=perplexities, colorPalette,
                    columnName = "clusters")
  plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
                    PCs=PCs, perplexities=perplexities, colorPalette,
                    columnName = "noColor")
  if(any(colnames(colData(sceObjectFiltered)) %in% "state")){
      plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
                        PCs=PCs, perplexities=perplexities, statePalette,
                        columnName = "state")
  }
  
  # Calculating cluster similarity and marker genes
  clustersSimilarityMatrix <- 
      calculateClustersSimilarity(cellsSimilarityMatrix, 
          sceObject = sceObjectFiltered,
          clusteringMethod = clusteringMethod)[[1]]
  
  plotClustersSimilarity(clustersSimilarityMatrix, 
                         sceObject = sceObjectFiltered,
                         dataDirectory = dataDirectory, 
                         experimentName = experimentName, 
                         colorPalette = colorPalette,
                         statePalette = statePalette,
                         clusteringMethod = clusteringMethod)
  
  rankGenes(sceObjectFiltered, clustersSimilarityMatrix, dataDirectory, 
            experimentName)
  
  # Exporting key matrices 
  exportMatrix(cellsSimilarityMatrix, dataDirectory, experimentName, 
               "cellsSimilarityMatrix")
  exportMatrix(clustersSimilarityMatrix, dataDirectory, experimentName, 
               "clustersSimilarityMatrix")
  
  return(sceObjectFiltered)
}

exportMatrix <- function(matrix, dataDirectory, experimentName, name){
  fileName <- paste0(experimentName, "_", name, ".csv")
  write.table(matrix, file=file.path(dataDirectory, "output_tables", fileName), 
              sep = ",")
}

initialisePath <- function(dataDirectory){
  # creates directories for further writing of results
  # names of directories are hardcoded
  # no idea is it good or bad
  
  graphsDirectory <- "pictures"
  markerGenesDirectory <- "marker_genes"
  tSNEDirectory <- "tsnes"
  outputDataDirectory <- "output_tables"
  tSNEPicturesDirectory <- "tSNE_pictures"
  
  
  dir.create(dataDirectory, showWarnings=F)
  dir.create(file.path(dataDirectory, graphsDirectory), showWarnings=F)
  dir.create(file.path(dataDirectory, graphsDirectory, tSNEPicturesDirectory), 
             showWarnings=F)
  dir.create(file.path(dataDirectory, markerGenesDirectory), showWarnings=F)
  dir.create(file.path(dataDirectory, tSNEDirectory), showWarnings=F)
  dir.create(file.path(dataDirectory, outputDataDirectory), showWarnings=F)
  
}

generateTSNECoordinates <- function(sceObject, dataDirectory, experimentName,
                                    randomSeed=42, cores=14, 
                                    PCs=c(4, 6, 8, 10, 20, 40, 50), 
                                    perplexities=c(30,40)){
  # generates several tSNE coordinates based on 
  # different perplexity and number of PCs
  # writes coordinates in special folder
  
  
  tSNEDirectory <- "tsnes"
  message(paste0("Running TSNEs using ", cores, " cores."))
  TSNEres <- getTSNEresults(exprs(sceObject), cores=cores, PCs=PCs, 
                            perplexities=perplexities, randomSeed=randomSeed)
  
  PCA <- rep(PCs, length(perplexities))
  perp <- rep(perplexities, each=length(PCs))
  
  outputDir <- file.path(dataDirectory, tSNEDirectory)
  filesList <- list.files(outputDir, pattern = "_tsne_coordinates_")
  deletedFiles <- sapply(filesList, function(fileName) deleteFile(fileName, 
                                                                  outputDir))
  for (i in 1:(length(PCs)*length(perplexities))){
    write.table(TSNEres[1, i][[1]],
                file=file.path(dataDirectory, tSNEDirectory,
                    paste0(experimentName,'_tsne_coordinates_', i, "_" ,
                           PCA[i], "PCs_", perp[i], "perp.tsv")), 
                quote=FALSE, sep='\t')
  }
  rm(tSNEDirectory, PCA, perp)
  return(TSNEres)
}

checkTSNEPicture <- function(tSNEResults, sceObject){
  # gives the unclustered picture of tSNE 
  # to be sure that normalization step was
  # successful
  
  tSNECoords <- tSNEResults[[1]]
  tSNECoords <- tSNECoords[rownames(tSNECoords) %in% 
                               colData(sceObject)$cellName, ]
  
  ggplot(tSNECoords, aes_string(x=names(tSNECoords)[1], 
                                y=names(tSNECoords)[2])) + 
    geom_point(size=I(1))
}

plotDistanceGraph <- function(tSNEResults, minNeighbours=5, tSNEIndex=3){
  # plots kNN distance graph 
  # for choosing right values 
  # of epsilon for further DBSCAN analysis
  
  
  tSNECoords <- tSNEResults[1,tSNEIndex][[1]]
  dbscan::kNNdistplot(tSNECoords, k=minNeighbours)
}

plotDistanceGraphWithEpsilon <- function(tSNEData, minNeighbours=5, 
                                         epsilon=1.2){
  # similar function as plotDistanceGraph,
  # but with already known epsilon value
  #
  
  dbscan::kNNdistplot(tSNEData, k=minNeighbours)
  abline(h=epsilon, lty=2)
  
}

plotTestClustering <- function(tSNEData, minNeighbours=5, 
                               epsilon=1.2){
  # plots test DBSCAN on one of the pictures
  # for being ensured that clustering will 
  # probably work successful
  
  dbscanResults <- fpc::dbscan(tSNEData, eps=epsilon, MinPts=minNeighbours)
  fviz_cluster(dbscanResults, tSNEData, ellipse=TRUE, geom="point", 
               legend="bottom") 
}

### This function calculates dbscan for all pSNE from TSNEtables with all 
### combinations of paramenters from epsilon and minPoints
### it does not set random seed. It allows to vary this parameter automatically
### it returns a matrix where columns are iterations
### number of iterations is equal to
### ncol(TSNEtables)*length(epsilon)*length(epsilon)

mkDbscan <- function(TSNEtables, cores = 14, epsilon = c(1.2, 1.5, 1.8),
                    minPoints = c(15, 20)){
    myCluster <- makeCluster(cores, # number of cores to use
                             type = "PSOCK") # type of cluster
    registerDoParallel(myCluster)
    dbscanResults <- foreach(i=rep(rep(1:ncol(TSNEtables),
                                   each=length(minPoints)),
                                   length(epsilon)),
                             eps=rep(epsilon,
                                 each=ncol(TSNEtables)*length(minPoints)), 
                             MinPts=rep(minPoints, 
                                    ncol(TSNEtables)*length(epsilon)), 
                         .combine='cbind') %dopar% {
                             fpc::dbscan(TSNEtables[1,i][[1]], eps=eps, 
                                         MinPts=MinPts)$cluster
                         }
    stopCluster(myCluster)
    return(dbscanResults)
}

runDBSCAN <- function(tSNEResults, sceObject, dataDirectory, experimentName, 
                      cores=14, epsilon=c(1.3, 1.4, 1.5), minPoints=c(3, 4)){
  # consensus clustering by DBSCAN
  # returning matrix of DBSCAN results
  
  outputDataDirectory <- "output_tables"
  # taking only cells from the sceObject
  for(i in 1:ncol(tSNEResults)){
      tmp <- tSNEResults[1,i][[1]]
      tSNEResults[1,i][[1]] <- tmp[colnames(sceObject),]
  }
  dbscanResults <- mkDbscan(tSNEResults, cores = cores, epsilon = epsilon, 
                             minPoints = minPoints) 
  dbscanResults <- t(dbscanResults)
  colnames(dbscanResults) <- colData(sceObject)$cellName
  write.table(dbscanResults, file.path(dataDirectory, outputDataDirectory,
              paste0(experimentName, "_dbscan_results.tsv")), quote = FALSE, 
              sep = "\t")
  rm(tmp)
  return(dbscanResults)
}

### This function calculates how many time a cell were not assigned to 
### any clusters by dbscan. ###
### it returns a data frame ###
mkOutlierScoreDf <- function(mat){
  
    outlierScoreDf <- as.data.frame(colnames(mat))
    colnames(outlierScoreDf) <- "cellName"
    outlierScoreDf <- dplyr::mutate(outlierScoreDf, outlierScore=NA)
    
    for(i in 1:ncol(mat)){
        vec <- mat[,i]
        outlierScoreDf$outlierScore[i] <- length(vec[vec == 0])
    }
    
    outlierScoreDf$outlierScorePer <- outlierScoreDf$outlierScore / nrow(mat)
    return(outlierScoreDf)
}

excludeOutliers <- function(dbscanMatrix, sceObject, threshold=0.3){
  # exclude outliers based on DBSCAN clustering
  # outliers are the cells which cannot be assigned
  # to any final cluster

  outlierInfo <- mkOutlierScoreDf(dbscanMatrix)
  colData <- colData(sceObject)[colData(sceObject)$cellName
                                %in% colnames(dbscanMatrix),]

  if(is.vector(colData)){
    colData <- DataFrame(cellName=colData, row.names=colData)
  }

  # Case if coldata has the outliers scores already
  colData$outlierScorePer <- NULL
  colData$outlierScore <- NULL

  colData <- merge(colData, outlierInfo,
                 by.x="cellName", by.y="cellName",
                 all.x=TRUE, all.y=TRUE, sort=FALSE)
  rownames(colData) <- colData$cellName

  numberOfCellsBefore <- dim(colData)[1]
  sceObject <- sceObject[,colData$outlierScorePer < threshold]
  colData <- colData[colData$outlierScorePer < threshold,]
  numberOfCellsAfter <- dim(colData)[1]

  dbscanMatrix <- dbscanMatrix[,outlierInfo$outlierScorePer < threshold]
  colData(sceObject) <- colData

  message(paste(numberOfCellsBefore - numberOfCellsAfter,
              "outliers were excluded from the SingleCellExperiment object.\n"))

  return(list(sceObject, dbscanMatrix))
}

mkSimMat <- function(mat, cores=14){
    
    myCluster <- makeCluster(cores, # number of cores to use
                             type = "PSOCK") # type of cluster
    registerDoParallel(myCluster)
    
    simMats <- foreach(i=1:nrow(mat)) %dopar% {
        simMat <- matrix(0, ncol=ncol(mat), nrow=ncol(mat))
        colnames(simMat) <- colnames(mat)
        rownames(simMat) <- colnames(mat)
        
        vec <- unique(mat[i,])
        for(j in vec[vec!=0]){
            selCol <- colnames(mat)[mat[i,] == j]
            simMat[rownames(simMat) %in% selCol, 
                    colnames(simMat) %in% selCol] <-
                simMat[rownames(simMat) %in% selCol, 
                        colnames(simMat) %in% selCol] + 1
        }
        rm(cl, selCol, i, j)
        return(simMat)
    }
    stopCluster(myCluster)
    
    simMat <- matrix(0, ncol=ncol(mat), nrow=ncol(mat))
    colnames(simMat) <- colnames(mat)
    rownames(simMat) <- colnames(mat)
    
    for(i in 1:nrow(mat)){
        simMat <- simMat + simMats[[i]]
    }
    
    rm(simMats)
    simMat <- simMat / nrow(mat)
    stopifnot(isSymmetric(simMat))
    
    return(simMat)
}

clusterCells <- function(dbscanMatrix, sceObject, clusterNumber=0, 
                         deepSplit, cores=14, 
                         clusteringMethod = "ward.D2") {
  # using hierarchical clustering for getting the number
  # of clusters we need from DBSCAN results
  # returning sceObject and cells correlation matrix
  
  cellsSimilarityMatrix <- mkSimMat(dbscanMatrix, cores=cores)
  
  distanceMatrix <- as.dist(sqrt((1-cellsSimilarityMatrix)/2))
  clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
  
  if(clusterNumber == 0){
    message(paste0("Assigning cells to clusters. DeepSplit = ", deepSplit))
    clusters <- unname(cutreeDynamic(clusteringTree, 
                                     distM=as.matrix(distanceMatrix), 
                                     verbose=0, deepSplit = deepSplit))  
  } else {
    message(paste0("Assigning cells to ", clusterNumber, " clusters."))
    clusters <- cutree(clusteringTree, k=clusterNumber)
  }
  
  colData(sceObject)$clusters <- factor(clusters)
  
  return(list(sceObject, cellsSimilarityMatrix))
}

plotCellSimilarity <- function(sceObject, cellsSimilarityMatrix, dataDirectory, 
                               experimentName, colorPalette="default",
                               statePalette="default", clusteringMethod="ward.D2",
                               orderClusters = FALSE,
                               plotPDF = TRUE,
                               returnPlot = FALSE,
                               width=7, height=6, onefile=FALSE, #pdf 
                               family, title, fonts, version,
                               paper, encoding, bg, fg, pointsize, 
                               pagecentre, colormodel,
                               useDingbats, useKerning, 
                               fillOddEven, compress,
                               color = colorRampPalette(rev( #pheatmap
                                   brewer.pal(n = 7, 
                                              name = "RdYlBu")))(100), 
                               kmeans_k = NA, breaks = NA, 
                               border_color = "grey60",
                               cellwidth = NA, cellheight = NA, 
                               scale = "none", 
                               cluster_rows = FALSE, #not default
                               cluster_cols = FALSE, #not default
                               clustering_distance_rows = "euclidean",
                               clustering_distance_cols = "euclidean", 
                               clustering_method = "complete",
                               clustering_callback = identity2, 
                               cutree_rows = NA, cutree_cols = NA,
                               treeheight_row = ifelse((
                                   class(cluster_rows) == "hclust") || 
                                       cluster_rows, 50, 0), 
                               treeheight_col = ifelse((
                                   class(cluster_cols) == "hclust") ||
                                       cluster_cols, 50, 0), 
                               legend = TRUE, 
                               legend_breaks = NA,
                               legend_labels = NA, 
                               annotation_row = NA, 
                               annotation_col = NA,
                               annotation = NA, 
                               annotation_colors = NA, 
                               annotation_legend = TRUE,
                               annotation_names_row = TRUE, 
                               annotation_names_col = TRUE,
                               drop_levels = TRUE, 
                               show_rownames = FALSE, #not default
                               show_colnames = FALSE, #not default
                               fontsize = 7.5, #not default (10)
                               fontsize_row = 0.03, 
                               fontsize_col = fontsize,
                               display_numbers = F, 
                               number_format = "%.2f", 
                               number_color = "grey30",
                               fontsize_number = 0.8 * fontsize, 
                               gaps_row = NULL, gaps_col = NULL,
                               labels_row = NULL, labels_col = NULL, 
                               filename = NA, 
                               widthHeatmap = NA, heightHeatmap = NA, #edited
                               silent = FALSE,
                               main = paste0("Cells similarity matrix ", 
                                             ncol(cellsSimilarityMatrix),
                                             " columns, ", 
                                             nrow(cellsSimilarityMatrix), 
                                             " rows."),
                               widthPNG = 500, heightPNG = 480 #png
                               ){
  # plots cells correlation matrix gained form 
  # clusterCells() function as the result of DBSCAN
  
  graphsDirectory <- "pictures"
  colData <- colData(sceObject)
  clustersNumber <- length(unique(colData$clusters))
  
  if(orderClusters == "name"){
      # Ordering expressionMatrixrix
      newOrder <- unname(unlist(sapply(levels(colData$clusters),
                                   function(cluster)
                                       orderCellsInCluster(cluster, 
                                           colData, 
                                           expressionMatrix,
                                           clusteringMethod=clusteringMethod))))
      
      cellsSimilarityMatrix <- cellsSimilarityMatrix[newOrder, newOrder]
      cluster_cols <- FALSE
      cluster_rows <- FALSE
  }else if(orderClusters == FALSE){
      distanceMatrix <- as.dist(sqrt((1-cellsSimilarityMatrix)/2))
      clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
      cluster_cols <- clusteringTree
      cluster_rows <- clusteringTree
  }else{
      message("Unknown option of orderClusters. Options are 'FALSE' or 'name'. 
Using the default version 'FALSE'.")
      distanceMatrix <- as.dist(sqrt((1-cellsSimilarityMatrix)/2))
      clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
      cluster_cols <- clusteringTree
      cluster_rows <- clusteringTree
  }
  
  annotationColors <- generateAnnotationColors(colData, colorPalette,
                                               statePalette)
  columnsToPlot <- switch(is.null(colData$state) + 1, c("clusters", "state"), 
                          c("clusters"))
  
  if(plotPDF){
      pdf(file=file.path(dataDirectory, graphsDirectory, paste(experimentName,
                    "cells_correlation", clustersNumber,"clusters.pdf", sep="_")),
          width=width, height=height, onefile=onefile, # not changed by default
          family=family, title=title, fonts=fonts, version=version,
          paper=paper, encoding=encoding, bg=bg, fg=fg, pointsize=pointsize, 
          pagecentre=pagecentre, colormodel=colormodel,
          useDingbats=useDingbats, useKerning=useKerning, 
          fillOddEven=fillOddEven, compress=compress)
  }else{
      message("Plot type is not pdf. Saving in png.")
      png(filename=file.path(dataDirectory, graphsDirectory, 
                             paste(experimentName,
                               "cells_correlation", clustersNumber,
                               "clusters.png", sep="_")),
          width = widthPNG, height = heightPNG, type = "cairo")
  }
  pheatmap(cellsSimilarityMatrix,
           show_colnames=show_colnames,
           show_rownames=show_rownames,
           annotation_col=as.data.frame(colData[columnsToPlot]),
           annotation_colors=annotationColors,
           fontsize_row=fontsize_row,
           cluster_cols=cluster_cols,
           cluster_rows=cluster_rows,
           fontsize=fontsize,
           main = main,
           color = color, kmeans_k = kmeans_k, # not changed by default
           breaks = breaks, 
           border_color = border_color,
           cellwidth = cellwidth, cellheight = cellheight, scale = scale, 
           clustering_distance_rows = clustering_distance_rows,
           clustering_distance_cols = clustering_distance_cols,
           clustering_method = clustering_method,
           clustering_callback = clustering_callback, 
           treeheight_row = treeheight_row, 
           treeheight_col = treeheight_col, legend = legend, 
           legend_breaks = legend_breaks,
           legend_labels = legend_labels, annotation_row = annotation_row,
           annotation = annotation, annotation_legend = annotation_legend,
           annotation_names_row = annotation_names_row, 
           annotation_names_col = annotation_names_col,
           drop_levels = drop_levels,
           fontsize_col = fontsize_col,
           display_numbers = display_numbers, 
           number_format = number_format, 
           number_color = number_color,
           fontsize_number = fontsize_numbere, gaps_row = gaps_row, 
           gaps_col = gaps_col, labels_row = labels_row, 
           labels_col = labels_col, filename = filename, width = widthHeatmap,
           height = heightHeatmap, silent = silent)
  dev.off()
  
  if(returnPlot){
      return(pheatmap(cellsSimilarityMatrix,
               show_colnames=show_colnames,
               show_rownames = show_rownames,
               annotation_col=as.data.frame(colData[columnsToPlot]),
               annotation_colors=annotationColors,
               fontsize_row=fontsize_row,
               cluster_cols=cluster_cols,
               cluster_rows=cluster_rows,
               fontsize=fontsize,
               main = main,
               color = color, kmeans_k = kmeans_k, # not changed by default
               breaks = breaks, 
               border_color = border_color,
               cellwidth = cellwidth, cellheight = cellheight, scale = scale, 
               clustering_distance_rows = clustering_distance_rows,
               clustering_distance_cols = clustering_distance_cols,
               clustering_method = clustering_method,
               clustering_callback = clustering_callback, 
               treeheight_row = treeheight_row, 
               treeheight_col = treeheight_col, legend = legend, 
               legend_breaks = legend_breaks,
               legend_labels = legend_labels, annotation_row = annotation_row,
               annotation = annotation, annotation_legend = annotation_legend,
               annotation_names_row = annotation_names_row, 
               annotation_names_col = annotation_names_col,
               drop_levels = drop_levels,
               fontsize_col = fontsize_col,
               display_numbers = display_numbers, 
               number_format = number_format, 
               number_color = number_color,
               fontsize_number = fontsize_numbere, gaps_row = gaps_row, 
               gaps_col = gaps_col, labels_row = labels_row, 
               labels_col = labels_col, filename = filename, width = widthHeatmap,
               height = heightHeatmap, silent = silent))
  }
}

plotClusteredTSNE <- function(sceObject, dataDirectory, experimentName,
                              colorPalette = "default",
                              PCs=c(4, 6, 8, 10, 20, 40, 50),
                              perplexities=c(30, 40),
                              columnName="clusters",
                              returnPlot = FALSE,
                              width=6, height=5, onefile=FALSE, #pdf 
                              family, title, fonts, version,
                              paper, encoding, bg, fg, pointsize, 
                              pagecentre, colormodel,
                              useDingbats, useKerning, 
                              fillOddEven, compress){
  # plots picture based on tSNE coordinates from
  # generateTSNECoordinates() and clustering results
  # from clusterCells() or runClustering() 
  
  tSNEDirectory <- "tsnes"
  graphsDirectory <- "pictures"
  graphsTSNEDirectory <- "tSNE_pictures"
  
  ### Plot all precalculated pSNEs to show your clusters ###
  
  if(columnName == "noColor"){
      numberElements <- NULL
  }else{
      numberElements <- length(unique(colData(sceObject)[,columnName]))
      colorPalette <- choosePalette(colorPalette, numberElements)
  }
  
  outputDir <- file.path(dataDirectory, graphsDirectory, graphsTSNEDirectory, 
                        paste("tSNE", numberElements, columnName, sep="_"))
  dir.create(outputDir, showWarnings = F)
  
  filesList <- list.files(outputDir, pattern = "_tsne_coordinates_")
  deletedFiles <- sapply(filesList, function(fileName) deleteFile(fileName, 
                                                                  outputDir))
  
  PCA <- rep(PCs, length(perplexities))
  perp <- rep(perplexities, each=length(PCs))
  
  tSNEplots <- rep(list(NA),(length(PCs)*length(perplexities)))
  
  for(i in 1:(length(PCs)*length(perplexities))){
    
    coordinatesName <- paste0(experimentName, '_tsne_coordinates_', i, "_",
                              PCA[i], "PCs_", perp[i], "perp")
      
    TSNEres <- read.delim(file.path(dataDirectory, tSNEDirectory, 
                                    paste0(coordinatesName, ".tsv")),
                          stringsAsFactors = FALSE)

    TSNEres <- TSNEres[rownames(TSNEres) %in% colData(sceObject)$cellName, ]
    
    if(columnName != "noColor"){
        TSNEres[columnName] <- factor(colData(sceObject)[,columnName])
    }
    
    pdf(file.path(dataDirectory, graphsDirectory, graphsTSNEDirectory,
                  paste("tSNE", numberElements, columnName, sep="_"),
                  paste0(coordinatesName, ".pdf")),
        width=width, height=height, onefile=onefile, # not changed by default
        family=family, title=title, fonts=fonts, version=version,
        paper=paper, encoding=encoding, bg=bg, fg=fg, pointsize=pointsize, 
        pagecentre=pagecentre, colormodel=colormodel,
        useDingbats=useDingbats, useKerning=useKerning, 
        fillOddEven=fillOddEven, compress=compress)
    if(columnName == "noColor"){
        tmp <- ggplot(TSNEres, aes_string(x=names(TSNEres)[1], 
                                          y=names(TSNEres)[2])) + 
            geom_point(size=I(1)) + theme_bw()
    }else{
        tmp <- ggplot(TSNEres, aes_string(x=names(TSNEres)[1], 
                                          y=names(TSNEres)[2], 
                                          color=columnName)) + 
            geom_point(size=I(1)) +
            scale_color_manual(values=colorPalette) + theme_bw()
    }
    print(tmp)
    dev.off()
    tSNEplots[[i]] <- tmp
  }
  if(returnPlot){
      return(tSNEplots)
  }
  rm(PCA, perp)
}

### This function returns a matrix with "protocells" representing clusters ###
### values show how much two "protocells" are similar ###
### 1 if clusters are very similar, 0 if very different ###

mkSimMed <- function(simMat, clusters){
    
    clusMed <- matrix(ncol=length(unique(clusters)), nrow=nrow(simMat))
    clusterNames <- levels(clusters)
    
    for(i in 1:ncol(clusMed)){
        clusMed[,i] <- rowMedians(simMat[,clusters == clusterNames[i]])
    }
    
    clusMed <- t(clusMed)
    
    simMed <- matrix(ncol=length(unique(clusters)), 
                     nrow=length(unique(clusters)))
    
    for(i in 1:ncol(simMed)){
        simMed[,i] <- rowMedians(clusMed[,clusters == clusterNames[i]])
    }
    
    # colnames(simMed) = 1:length(unique(clusters))
    # rownames(simMed) = 1:length(unique(clusters))
    
    colnames(simMed) <- clusterNames
    rownames(simMed) <- clusterNames
    
    return(simMed)
}

calculateClustersSimilarity <- function(cellsSimilarityMatrix, sceObject, 
                                        clusteringMethod){
  
    # Calculating cluster similarity for plotting picture 
    # and ranking genes result is the square matrix with
    # dimension equal to number of cluster. Numbers in matrix
    # are similarity between cluster.
    
    clusters <- colData(sceObject)$clusters
    clustersNumber <- length(unique(clusters))
    clustersNames <- levels(clusters)
    
    # Calculating matrix
    clustersSimilarityMatrix <- 
      mkSimMed(simMat=as.matrix(cellsSimilarityMatrix), clusters=clusters)
    
    # Plotting matrix
    distanceMatrix <- as.dist(sqrt((1-clustersSimilarityMatrix)/2))
    clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
    
    clustersSimOrdered <- data.frame(clusterNames = clustersNames, 
                                     clusterIndexes = 1:clustersNumber)
    rownames(clustersSimOrdered) <- clustersSimOrdered$clusterIndexes
    clustersSimOrdered <- clustersSimOrdered[clusteringTree$order,]

    return(list(clustersSimilarityMatrix, clustersSimOrdered$clusterNames))

}

plotClustersSimilarity <- function(clustersSimilarityMatrix, sceObject,
                                   dataDirectory,
                                   experimentName, colorPalette, 
                                   statePalette,
                                   clusteringMethod,
                                   returnPlot = FALSE,
                                   width=7, height=5.5, onefile=FALSE, #pdf 
                                   family, title, fonts, version,
                                   paper, encoding, bg, fg, pointsize, 
                                   pagecentre, colormodel,
                                   useDingbats, useKerning, 
                                   fillOddEven, compress,
                                   color = colorRampPalette(rev( #pheatmap
                                       brewer.pal(n = 7, 
                                                  name = "RdYlBu")))(100), 
                                   kmeans_k = NA, breaks = NA, 
                                   border_color = "grey60",
                                   cellwidth = NA, cellheight = NA, 
                                   scale = "none", cluster_rows = TRUE,
                                   cluster_cols = TRUE, 
                                   clustering_distance_rows = "euclidean",
                                   clustering_distance_cols = "euclidean", 
                                   clustering_method = "complete",
                                   clustering_callback = identity2, 
                                   cutree_rows = NA, cutree_cols = NA,
                                   treeheight_row = ifelse((
                                       class(cluster_rows) == "hclust") || 
                                           cluster_rows, 50, 0), 
                                   treeheight_col = ifelse((
                                       class(cluster_cols) == "hclust") ||
                                           cluster_cols, 50, 0), 
                                   legend = TRUE, 
                                   legend_breaks = NA,
                                   legend_labels = NA, 
                                   annotation_row = NA, 
                                   annotation_col = NA,
                                   annotation = NA, 
                                   annotation_colors = NA, 
                                   annotation_legend = TRUE,
                                   annotation_names_row = TRUE, 
                                   annotation_names_col = TRUE,
                                   drop_levels = TRUE, 
                                   show_rownames = TRUE, 
                                   show_colnames = TRUE,
                                   fontsize = 7.5, #not default (10)
                                   fontsize_row = fontsize, 
                                   fontsize_col = fontsize,
                                   display_numbers = F, 
                                   number_format = "%.2f", 
                                   number_color = "grey30",
                                   fontsize_number = 0.8 * fontsize, 
                                   gaps_row = NULL, gaps_col = NULL,
                                   labels_row = NULL, labels_col = NULL, 
                                   filename = NA, widthHeatmap = NA,
                                   heightHeatmap = NA, silent = FALSE,
                                   main = 
                                       paste0("Clusters similarity matrix ", 
                                              ncol(clustersSimilarityMatrix),
                                              " columns, ", nrow(clustersSimilarityMatrix), 
                                              " rows.")){
    clusters <- colData(sceObject)$clusters
    clustersNumber <- length(unique(clusters))
    clustersNames <- levels(clusters)
    dataDirectory <- dataDirectory
    experimentName <- experimentName
    graphsDirectory <- "pictures"
    
    # Plotting matrix
    distanceMatrix <- as.dist(sqrt((1-clustersSimilarityMatrix)/2))
    clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
    
    colDataSimilarity <- data.frame(clusters = clustersNames)
    rownames(colDataSimilarity) <- colDataSimilarity$clusters
    
    annotationColors <- generateAnnotationColors(colData(sceObject), 
                                                 colorPalette,
                                                 statePalette)
    message("\nSaving a heatmap with the clusters similarity matrix.")
    pdf(file.path(dataDirectory, graphsDirectory,
                  paste(experimentName,"clusters_similarity", clustersNumber, 
                        "clusters.pdf", sep="_")),
        width=width, height=height, onefile=onefile, # not changed by default
        family=family, title=title, fonts=fonts, version=version,
        paper=paper, encoding=encoding, bg=bg, fg=fg, pointsize=pointsize, 
        pagecentre=pagecentre, colormodel=colormodel,
        useDingbats=useDingbats, useKerning=useKerning, 
        fillOddEven=fillOddEven, compress=compress)
    pheatmap(clustersSimilarityMatrix,
             show_colnames=show_colnames,
             show_rownames=show_rownames,
             annotation_col=colDataSimilarity, 
             annotation_colors=annotationColors,
             cluster_cols=clusteringTree,
             cluster_rows=clusteringTree,
             fontsize=fontsize,
             main = main,
             color = color, kmeans_k = kmeans_k, # not changed by default
             breaks = breaks, 
             border_color = border_color,
             cellwidth = cellwidth, cellheight = cellheight, scale = scale, 
             clustering_distance_rows = clustering_distance_rows,
             clustering_distance_cols = clustering_distance_cols,
             clustering_method = clustering_method,
             clustering_callback = clustering_callback, 
             treeheight_row = treeheight_row, 
             treeheight_col = treeheight_col, legend = legend, 
             legend_breaks = legend_breaks,
             legend_labels = legend_labels, annotation_row = annotation_row,
             annotation = annotation, annotation_legend = annotation_legend,
             annotation_names_row = annotation_names_row, 
             annotation_names_col = annotation_names_col,
             drop_levels = drop_levels,
             fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             display_numbers = display_numbers, 
             number_format = number_format, 
             number_color = number_color,
             fontsize_number = fontsize_numbere, gaps_row = gaps_row, 
             gaps_col = gaps_col, labels_row = labels_row, 
             labels_col = labels_col, filename = filename, width = widthHeatmap,
             height = heightHeatmap, silent = silent)
    dev.off()
    
    if(returnPlot){
        pheatmap(clustersSimilarityMatrix,
                 show_colnames=show_colnames,
                 show_rownames=show_rownames,
                 annotation_col=colDataSimilarity, 
                 annotation_colors=annotationColors,
                 cluster_cols=clusteringTree,
                 cluster_rows=clusteringTree,
                 fontsize=fontsize,
                 main = main,
                 color = color, kmeans_k = kmeans_k, # not changed by default
                 breaks = breaks, 
                 border_color = border_color,
                 cellwidth = cellwidth, cellheight = cellheight, scale = scale, 
                 clustering_distance_rows = clustering_distance_rows,
                 clustering_distance_cols = clustering_distance_cols,
                 clustering_method = clustering_method,
                 clustering_callback = clustering_callback, 
                 treeheight_row = treeheight_row, 
                 treeheight_col = treeheight_col, legend = legend, 
                 legend_breaks = legend_breaks,
                 legend_labels = legend_labels, annotation_row = annotation_row,
                 annotation = annotation, annotation_legend = annotation_legend,
                 annotation_names_row = annotation_names_row, 
                 annotation_names_col = annotation_names_col,
                 drop_levels = drop_levels,
                 fontsize_row = fontsize_row, fontsize_col = fontsize_col,
                 display_numbers = display_numbers, 
                 number_format = number_format, 
                 number_color = number_color,
                 fontsize_number = fontsize_numbere, gaps_row = gaps_row, 
                 gaps_col = gaps_col, labels_row = labels_row, 
                 labels_col = labels_col, filename = filename, width = widthHeatmap,
                 height = heightHeatmap, silent = silent)
    }
}

# rankGenesInternal() saves n files in your outputDir, n=number of groups
deleteFile <- function(fileName, outputDir){
    file.remove(file.path(outputDir, fileName))
}

rankGenesInternal <- function(expr, colData, column, simMed, outputDir, 
                             experimentName){
  
    message("Ranking marker genes for each cluster.")
  
    filesList <- list.files(outputDir, pattern = "_genes.tsv")
    deletedFiles <- sapply(filesList, function(fileName) deleteFile(fileName, 
                                                                    outputDir))
  
    stopifnot(all(colnames(expr) == rownames(colData)))
    groups <- unique(colData[,column])
    simMed = simMed + 0.05
    for(i in 1:length(groups)){
        message(paste("Working on the file", i))
        tTestPval <- data.frame(Gene=rownames(expr))
        otherGroups <- groups[groups!=groups[i]]
        
        for(k in 1:length(otherGroups)){
            tTestPval[,paste0("vs_", otherGroups[k])] <- NA
            x <- expr[,colData[,c(column)] == groups[i]]
            y <- expr[,colData[,c(column)] == otherGroups[k]]
            t <- (rowMeans(x) - rowMeans(y)) / 
                sqrt(apply(expr, 1, var)*(1/ncol(x) + 1/ncol(y)))
            df <- ncol(x) + ncol(y) - 2
            tTestPval[, paste0("vs_", otherGroups[k])] <- 
                pt(t, df, lower.tail=FALSE)
        }
        
        tTestFDR <- data.frame(Gene=tTestPval$Gene)
        for(l in 1:length(otherGroups)){
            tTestFDR[,paste0("vs_", otherGroups[l])] <- 
                p.adjust(tTestPval[,paste0("vs_", otherGroups[l])], 
                         method="fdr")
        }
        
        submat <- as.matrix(tTestFDR[,2:(length(otherGroups)+1)])
        tTestFDR$mean_log10_fdr <- rowMeans(log10(submat+1e-300))
        tTestFDR$n_05 <- apply(submat, 1, function(x) 
            length(x[!is.na(x) & x < 0.05]))
        
        weights <- simMed[i, otherGroups]
        tTestFDR$score <- apply(submat, 1, function(x)
            sum(-log10(x+1e-300) * weights) / ncol(submat))
        
        tTestFDR <- tTestFDR[order(tTestFDR$score,decreasing=TRUE),]
        
        write.table(tTestFDR, file.path(outputDir,
                                        paste0(experimentName, "_cluster_", 
                                               groups[i], "_genes.tsv")), 
                    col.names=TRUE, row.names=FALSE, quote=FALSE, 
                    sep="\t")
    }
    
    rm(tTestFDR, tTestPval, i, k, l, x, y, t, df, otherGroups, submat, 
       weights, groups)
}

rankGenes <- function(sceObject, clustersSimilarityMatrix, dataDirectory, 
                      experimentName, column="clusters"){
  # Generates marker genes for each cluster. Creates tables in the 
  # marker genes directory, each table represents its own cluster.
  
  markerGenesDirectory <- "marker_genes"
  rankGenesInternal(exprs(sceObject), colData(sceObject), column, 
             clustersSimilarityMatrix,
             file.path(dataDirectory, markerGenesDirectory), experimentName)
  
}

### This function reads results of rankGenes() from outputDir ###
### and select equal number of nTop genes ###
### it returns a vector of gene names ###
### gene names are not necessarily unique ###

getMarkerGenes <- function(dataDirectory, sceObject, genesNumber=14, 
                           experimentName){
    
    markerGenesDirectory <- "marker_genes"
    numberOfClusters <- length(unique(colData(sceObject)$clusters))
    dir = file.path(dataDirectory, markerGenesDirectory)
    nTop = genesNumber
    clusters = unique(colData(sceObject)$clusters)
    
    markersClusters <- as.data.frame(matrix(ncol = 2,
                                            nrow = nTop*numberOfClusters))
    colnames(markersClusters) = c("geneName", "clusters")
    
    (fns <- list.files(dir, pattern = "_genes.tsv", full.names = FALSE))
    if(length(fns) != numberOfClusters){
        message(paste("Something wrong with number of files. 
It is supposed to be equal to number of clusters:", numberOfClusters))
        message(paste("Returning the marker genes from 
first", clusters, "clusters."))
        runUntil = numberOfClusters
    }else{
        runUntil = length(fns)
    }
    for(i in 1:runUntil){
        tmpAll <- read.delim(file.path(dir, paste(experimentName, 
                                                  "cluster", clusters[i], "genes.tsv", sep="_")), 
                             stringsAsFactors = FALSE)
        markersClusters$clusters[(nTop*(i-1)+1):(nTop*i)] <- as.character(clusters[i])
        markersClusters$geneName[(nTop*(i-1)+1):(nTop*i)] <- tmpAll$Gene[1:nTop]
    }
    markersClusters <- markersClusters[!duplicated(markersClusters$geneName),]
    return(markersClusters)
}

orderCellsInCluster <- function(cluster, colData, mtx, 
                                clusteringMethod="ward.D2"){
  # Order cells according to clustering results
  # Uses for ordering matrix to further plot it with pheatmap()
  
  cells <- colData[colData$clusters == cluster, ]$cellName
  if(length(cells) > 2){
      tree <- hclust(dist(t(mtx[, cells])), method=clusteringMethod)
      return(cells[tree$order])
  }else{
      return(cells)
  }
  
}

orderGenesInCluster <- function(cluster, markersClusters, mtx,
                                clusteringMethod="ward.D2"){
    # Order cells according to clustering results
    # Uses for ordering matrix to further plot it with pheatmap()

    genes <- markersClusters[markersClusters$clusters == cluster, ]$geneName
    if(length(genes) > 2){
        tree <- hclust(dist(mtx[genes, ]), method=clusteringMethod)
        return(genes[tree$order])
    }else{
        return(genes)
    }
    
    
}

generateAnnotationColors <- function(colData, colorPaletteParameter,
                                     statePalette){
  
  clusters <- levels(colData$clusters)
  states <- unique(colData$state)
  clusterNumber <- length(unique(colData$clusters))
  
  colorAnnotationClusters <- choosePalette(colorPaletteParameter, clusterNumber)
  #colorAnnotationState <- chooseStatePalette(length(states))
  colorAnnotationState <- choosePalette(statePalette, length(states))
  names(colorAnnotationState) <- states
  names(colorAnnotationClusters) <- clusters
  
  return(list(state=colorAnnotationState, clusters=colorAnnotationClusters))
}

plotCellHeatmap <- function(markersClusters, sceObject, dataDirectory, 
                            experimentName,
                            fileName, meanCentered=TRUE, colorPalette="default", 
                            statePalette="default", clusteringMethod="ward.D2",
                            orderClusters = FALSE, #FALSE, TRUE, name, similarity
                            orderGenes = FALSE, # FALSE, TRUE (will be ordered the same as clusters)
                            returnPlot = FALSE,
                            width=10, height=8.5, onefile=FALSE, #pdf 
                            family, title, fonts, version,
                            paper, encoding, bg, fg, pointsize, 
                            pagecentre, colormodel,
                            useDingbats, useKerning, 
                            fillOddEven, compress,
                            color = colorRampPalette(c("#023b84","#4b97fc", 
                                            "#c9d9ef","#FEE395", #not default
                                            "#F4794E", "#D73027",#pheatmap
                                            "#a31008","#7a0f09"))(100), 
                            kmeans_k = NA, breaks = NA, 
                            border_color = "grey60",
                            cellwidth = NA, cellheight = NA, 
                            scale = "none", cluster_rows = TRUE,
                            cluster_cols = FALSE, # not original default
                            clustering_distance_rows = "euclidean",
                            clustering_distance_cols = "euclidean", 
                            clustering_method = "complete",
                            clustering_callback = identity2, 
                            cutree_rows = NA, cutree_cols = NA,
                            treeheight_row = ifelse((
                                class(cluster_rows) == "hclust") || 
                                    cluster_rows, 50, 0), 
                            treeheight_col = ifelse((
                                class(cluster_cols) == "hclust") ||
                                    cluster_cols, 50, 0), 
                            legend = TRUE, 
                            legend_breaks = NA,
                            legend_labels = NA, 
                            annotation_row = NA, 
                            annotation_col = NA,
                            annotation = NA, 
                            annotation_colors = NA, 
                            annotation_legend = TRUE,
                            annotation_names_row = TRUE, 
                            annotation_names_col = TRUE,
                            drop_levels = TRUE, 
                            show_rownames = TRUE, 
                            show_colnames = FALSE, # not original default (T)
                            main = NA,
                            fontsize = 7.5, # not original default (10)
                            fontsize_row = 8, # not original default
                            fontsize_col = fontsize,
                            display_numbers = F, 
                            number_format = "%.2f", 
                            number_color = "grey30",
                            fontsize_number = 0.8 * fontsize, 
                            gaps_row = NULL, gaps_col = NULL,
                            labels_row = NULL, labels_col = NULL, 
                            filename = NA, widthHeatmap = NA,
                            heightHeatmap = NA, silent = FALSE){
  # plots correlation between cells between clusters
  
  colData <- colData(sceObject)
  expressionMatrix <- exprs(sceObject)[rownames(exprs(sceObject)) %in% 
                                           markersClusters$geneName,]
  
  if(meanCentered == TRUE){
    meanRows <- rowSums(expressionMatrix)/ncol(expressionMatrix)
    expressionMatrix <- expressionMatrix - meanRows
  }
  
  if(orderClusters == FALSE & orderGenes == TRUE){
      message("Genes cannot be ordered without clusters.
              Returning heatmap with orderCluster = similarity")
      orderClusters <- "similarity"
  }
  if(orderClusters == "name"){
      # Ordering expressionMatrixrix
      newOrder <- unname(unlist(sapply(levels(colData$clusters),
                                   function(cluster)
                                         orderCellsInCluster(cluster, 
                                         colData, 
                                         expressionMatrix,
                                         clusteringMethod=clusteringMethod))))
      
      expressionMatrix <- expressionMatrix[, newOrder]
      cluster_cols <- FALSE
      if(orderGenes){
          newOrder <- unname(unlist(sapply(levels(colData$clusters),
                                   function(cluster)
                                          orderGenesInCluster(cluster, 
                                          markersClusters, 
                                          expressionMatrix,
                                          clusteringMethod=clusteringMethod))))
          expressionMatrix <- expressionMatrix[newOrder, ]
          cluster_rows <- FALSE
      }
  }else if((orderClusters == TRUE) | (orderClusters == "similarity")){
      cellsSimilarityMatrix <- read.delim(file.path(dataDirectory, 
                                              "output_tables",
                                              paste0(experimentName, 
                                              "_cellsSimilarityMatrix.csv")),
                                          stringsAsFactors = FALSE, 
                                          header = TRUE, 
                                          sep = ",")
      
      clustersSimOrdered <- calculateClustersSimilarity(cellsSimilarityMatrix, 
                                                        sceObject, 
                                                        clusteringMethod)[[2]]
      
      newOrder <- unname(unlist(sapply(clustersSimOrdered,
                                       function(cluster)
                                           orderCellsInCluster(cluster, 
                                           colData, 
                                           expressionMatrix,
                                           clusteringMethod=clusteringMethod))))
      
      expressionMatrix <- expressionMatrix[, newOrder]
      cluster_cols <- FALSE
      if(orderGenes){
          newOrder <- unname(unlist(sapply(clustersSimOrdered,
                                       function(cluster)
                                          orderGenesInCluster(cluster, 
                                          markersClusters, 
                                          expressionMatrix,
                                          clusteringMethod=clusteringMethod))))
          expressionMatrix <- expressionMatrix[newOrder, ]
          cluster_rows <- FALSE
      }
  }else if(orderClusters == FALSE){
      distanceMatrix <- dist(t(expressionMatrix))
      cluster_cols <- hclust(distanceMatrix, method="ward.D2")
  }else{
      message("Unknown option of orderClusters. Options are 'TRUE' ('similarity'), 
'FALSE', 'name'. Using the default version 'FALSE'.")
      distanceMatrix <- dist(t(expressionMatrix))
      cluster_cols <- hclust(distanceMatrix, method="ward.D2")
  }
  
  if(!orderGenes){
      cluster_rows <- hclust(dist(expressionMatrix), method="ward.D2")
  }
  
  annotationColors <- generateAnnotationColors(colData, colorPalette, 
                                               statePalette)
  columnsToPlot <- switch(is.null(colData$state) + 1, c("clusters", "state"),
                          c("clusters"))
  
  if(is.null(colData$clusters)){
      annCol <- switch(is.null(colData$state) + 1, 
                       as.data.frame(colData["state"]), NA)
      annColors <- switch(is.null(colData$state) + 1, 
                          annotationColors[names(annotationColors) == "state"],
                          NA)
  }else{
      annCol <- as.data.frame(colData[columnsToPlot])
      annColors <- annotationColors
  }
  
  pdf(file.path(dataDirectory, "pictures", 
                paste0(experimentName, "_", fileName, ".pdf")),
      width=width, height=height, onefile=onefile, # not changed by default
      family=family, title=title, fonts=fonts, version=version,
      paper=paper, encoding=encoding, bg=bg, fg=fg, pointsize=pointsize, 
      pagecentre=pagecentre, colormodel=colormodel,
      useDingbats=useDingbats, useKerning=useKerning, 
      fillOddEven=fillOddEven, compress=compress)
  pheatmap(expressionMatrix,
           show_colnames=show_colnames,
           show_rownames=show_rownames,
           annotation_col=annCol,
           annotation_colors=annColors,
           fontsize_row=fontsize_row,
           cluster_cols=cluster_cols,
           main=main,
           cluster_rows=cluster_rows,
           color=color,
           fontsize=fontsize, kmeans_k = kmeans_k, # not changed by default
           breaks = breaks, 
           border_color = border_color,
           cellwidth = cellwidth, cellheight = cellheight, scale = scale, 
           clustering_distance_rows = clustering_distance_rows,
           clustering_distance_cols = clustering_distance_cols,
           clustering_method = clustering_method,
           clustering_callback = clustering_callback, 
           treeheight_row = treeheight_row, 
           treeheight_col = treeheight_col, legend = legend, 
           legend_breaks = legend_breaks,
           legend_labels = legend_labels, 
           annotation = annotation,
           annotation_row = annotation_row,
           annotation_legend = annotation_legend,
           annotation_names_row = annotation_names_row, 
           annotation_names_col = annotation_names_col,
           drop_levels = drop_levels,
           fontsize_col = fontsize_col,
           display_numbers = display_numbers, 
           number_format = number_format, 
           number_color = number_color,
           fontsize_number = fontsize_numbere, gaps_row = gaps_row, 
           gaps_col = gaps_col, labels_row = labels_row, 
           labels_col = labels_col, filename = filename, width = widthHeatmap,
           height = heightHeatmap, silent = silent)
  dev.off()
  
  exportMatrix(expressionMatrix, dataDirectory, experimentName, fileName)
  
  if(returnPlot){
      return(pheatmap(expressionMatrix,
                      show_colnames=show_colnames,
                      show_rownames=show_rownames,
                      annotation_col=annCol,
                      annotation_colors=annColors,
                      fontsize_row=fontsize_row,
                      cluster_cols=cluster_cols,
                      main=main,
                      cluster_rows=cluster_rows,
                      color=color,
                      fontsize=fontsize, kmeans_k = kmeans_k, # not changed by default
                      breaks = breaks, 
                      border_color = border_color,
                      cellwidth = cellwidth, cellheight = cellheight, scale = scale, 
                      clustering_distance_rows = clustering_distance_rows,
                      clustering_distance_cols = clustering_distance_cols,
                      clustering_method = clustering_method,
                      clustering_callback = clustering_callback, 
                      treeheight_row = treeheight_row, 
                      treeheight_col = treeheight_col, legend = legend, 
                      legend_breaks = legend_breaks,
                      legend_labels = legend_labels, 
                      annotation = annotation,
                      annotation_row = annotation_row,
                      annotation_legend = annotation_legend,
                      annotation_names_row = annotation_names_row, 
                      annotation_names_col = annotation_names_col,
                      drop_levels = drop_levels,
                      fontsize_col = fontsize_col,
                      display_numbers = display_numbers, 
                      number_format = number_format, 
                      number_color = number_color,
                      fontsize_number = fontsize_numbere, gaps_row = gaps_row, 
                      gaps_col = gaps_col, labels_row = labels_row, 
                      labels_col = labels_col, filename = filename, width = widthHeatmap,
                      height = heightHeatmap, silent = silent))
  }
}

exportData <- function(sceObject, dataDirectory, experimentName){
  # exports all the data from workflow, including .RData files
  
  outputDataDirectory <- "output_tables"
  
  ################ EXPORT MATRIX, COLDATA, ROWDATA, MAT, ALL OBJECTS
  write.table(exprs(sceObject), file=file.path(dataDirectory, 
            outputDataDirectory, paste0(experimentName, "_", 
                                        "expression_matrix.tsv")), sep="\t", 
            row.names = TRUE, quote = FALSE, col.names = TRUE)
  write.table(colData(sceObject), file=file.path(dataDirectory, 
            outputDataDirectory, paste0(experimentName, "_", 
                                        "colData.tsv")), sep="\t", 
            row.names = TRUE, quote = FALSE, col.names = TRUE)
  write.table(rowData(sceObject), file=file.path(dataDirectory, 
            outputDataDirectory, paste0(experimentName, "_", 
                                        "rowData.tsv")), sep="\t", 
            row.names = TRUE, quote = FALSE, col.names = TRUE)
  save.image(file=file.path(dataDirectory, outputDataDirectory, 
                        paste0(experimentName, "_", "full_workspace.RData")))
  
}

runClustering <- function(tSNEResults, # for deleteOutliers = FALSE
                          sceObject, dataDirectory, 
                          experimentName, epsilon=c(1.3, 1.4, 1.5), 
                          minPoints=c(3, 4), k=0, deepSplit=4,
                          clusteringMethod = "ward.D2",
                          cores=14,
                          deleteOutliers = TRUE,
                          PCs=c(4, 6, 8, 10, 20, 40, 50), 
                          perplexities=c(30,40), # for deleteOutliers = TRUE
                          randomSeed = 42){
  # Aggregates all the clustering parts. Takes tSNE coordinates and gives
  # results of final clustering: sceObject and cell correlation matrix
  
  # doing dbscan
  message("Running dbscan using ", cores, " cores.")
  dbscanResults <- runDBSCAN(tSNEResults, sceObject, dataDirectory, 
                             experimentName, epsilon=epsilon, 
                             minPoints=minPoints, cores=cores)
  if(deleteOutliers){
      # filtering by cluster outliers
      message("Excluding clustering outliers.")
      filteredResults <- excludeOutliers(dbscanResults, sceObject)
      sceObjectFiltered <- filteredResults[[1]]
      #dbscanResultsFiltered <- filteredResults[[2]]
      message("Getting TSNE coordinates for the filtered sceObject.")
      tSNEResultsFiltered <- generateTSNECoordinates(sceObjectFiltered, 
                                                     dataDirectory, 
                                                     experimentName, PCs=PCs, 
                                                     perplexities=perplexities,
                                                     randomSeed=randomSeed)
      
      dbscanResultsFiltered <- runDBSCAN(tSNEResultsFiltered, sceObjectFiltered, 
                                         dataDirectory, 
                                         experimentName, epsilon=epsilon, 
                                         minPoints=minPoints, cores=cores)
  }else{
      sceObjectFiltered = sceObject
      dbscanResultsFiltered = dbscanResults
  }
  
  # assigning cells to cluster
  clusteringResults <- clusterCells(dbscanResultsFiltered, sceObjectFiltered, 
                                    clusterNumber=k, deepSplit=deepSplit,
                                    clusteringMethod=clusteringMethod,
                                    cores=cores)
  sceObjectFiltered <- clusteringResults[[1]]
  cellsSimilarityMatrix <- clusteringResults[[2]]
  
  return(list(sceObjectFiltered, cellsSimilarityMatrix))
}

getKEGGGenes <- function(pathwayID, sceObject, species="mmu"){
  # Extracting genes from KEGG database by pathway ID. Returns only the genes
  # that found in expression matrix
  
  keggOutput <- keggGet(paste0(species, pathwayID))[[1]]$GENE
  keggOutput <- keggOutput[seq(2, length(keggOutput), by = 2)]
  geneList <- unique(unname(sapply(keggOutput, 
                                   function(x) strsplit(x, ";")[[1]][1])))
  filteredGeneList <- geneList[geneList %in% rownames(sceObject)]
  
  message(paste0("Taking ", length(filteredGeneList), " genes out of ", 
             length(geneList), " from KEGG pathway: ", species, pathwayID, 
             " ", keggGet(paste0(species, pathwayID))[[1]]$NAME))
  
  return(filteredGeneList)

}

plotGeneExpression <- function(geneName, experimentName, dataDirectory, 
                               sceObject, tSNEpicture=1,
                               returnPlot = FALSE,
                               width=6, height=5, onefile=FALSE, #pdf 
                               family, title, fonts, version,
                               paper, encoding, bg, fg, pointsize, 
                               pagecentre, colormodel,
                               useDingbats, useKerning, 
                               fillOddEven, compress){
  # Plots expression of the gene on tSNE coordinates.
  
  # filename with tsne is hardcoded
  
  experimentName <- experimentName
  dataDirectory <- dataDirectory
  tSNEDirectory <- "tsnes"
  graphsDirectory <- "pictures"
  
  ### Plot all precalculated pSNEs to show your clusters ###
  
  clustersNumber <- length(unique(colData(sceObject)$clusters))
  
  coordsName <- list.files(file.path(dataDirectory, tSNEDirectory),
                          pattern = paste0(experimentName,'_tsne_coordinates_',
                                           tSNEpicture, "_"))
  
  tSNECoords <- read.delim(file.path(dataDirectory, tSNEDirectory, coordsName),
                        stringsAsFactors=FALSE)
  
  tSNECoords <- tSNECoords[colData(sceObject)$cellName, ]

  if(!geneName %in% rownames(exprs(sceObject))){
    print("Gene is not found in expression matrix")
  }
  
  stopifnot(all(rownames(tSNECoords) == colnames(sceObject)))
  tSNECoords$expression <- exprs(sceObject)[geneName, ]
  
  pdf(file.path(dataDirectory, graphsDirectory, paste0(paste(experimentName, 
                    "tSNE", clustersNumber, "clusters" , geneName,
                    "tSNEpicture", tSNEpicture, 
                    sep="_"), ".pdf")), 
      width=width, height=height, onefile=onefile, # not changed by default
      family=family, title=title, fonts=fonts, version=version,
      paper=paper, encoding=encoding, bg=bg, fg=fg, pointsize=pointsize, 
      pagecentre=pagecentre, colormodel=colormodel,
      useDingbats=useDingbats, useKerning=useKerning, 
      fillOddEven=fillOddEven, compress=compress)
  tmp <- ggplot(tSNECoords, aes(x=tSNECoords[,1], 
                                y=tSNECoords[,2], color=expression)) + 
    geom_point(size=I(1)) + theme_bw() +
    scale_colour_gradientn(colours=alpha(colorRampPalette(c("grey","red",
                                                        "#7a0f09", 
                                                        "black"))(100), 0.8)) +
      ggtitle(geneName)
  #brewer.pal(9, "OrRd")[0:9]
  print(tmp)
  
  dev.off()
  
  if(returnPlot){
      return(tmp)
  }
}

exportClusteringResults <- function(sceObject, dataDirectory,
                                    experimentName, fileName){
  # Exports clustering results into the table. Rows are cell names.
  
  tableData <- DataFrame(clusters = colData(sceObject)$clusters, 
                         row.names = colData(sceObject)$cellName) 
  write.table(tableData, 
              file = file.path(dataDirectory, "output_tables",
                               paste0(experimentName,"_", fileName)), 
              sep = "\t", quote = FALSE)
}

addClusteringManually <- function(fileName, sceObject, dataDirectory, 
                                  experimentName, columnName = "clusters"){
  
  # Replaces clustering results with clustering provided in the user table. 
  # Table must contain all the cells from experiment and nothing more.
  
  tableData <- read.table(file.path(dataDirectory, "output_tables", 
                                paste0(experimentName,"_", fileName)), sep="\t")
  
  if(all(rownames(colData(sceObject)) %in% rownames(tableData))){
    if(ncol(tableData) == 1){
      colData(sceObject)$clusters <- 
          factor(tableData[rownames(colData(sceObject)),])
    } else {
      colData(sceObject)$clusters <- 
          factor(tableData[rownames(colData(sceObject)),][,columnName])
    }
    
    return(sceObject)
  } else {
    message("Rownames in colData are not equal to rownames in table. 
Returning initial SCE object.")
    return(sceObject)
  }
}

# This function assumes that rownames(countMatrix) are either ENSEMBL IDs or 
# or SYMBOLs. It will return a rowData with the same rownames as in countMatrix
# but genes which are not ENSEMBL IDs or SYMBOLs will not receive the annotation.
# Both types are possible in one matrix but not in one gene.
# Genes like Lrrc8a_ENSMUSG00000007476 will not be annotated.
# Iformation about cell surface localization will be used in the Shine application.
# It is useful to find cell surface markers. But this part of the function 
# takes some time, so if you need this info, you can disable this option by 
# cellSurface = FALSE to speed up the calculations.
annotateGenes <- function(countMatrix, species = "mmu",
                          genomeAnnot, ensemblPattern, rowData = NULL){
    if(missing(species) & (missing(genomeAnnot) | missing(ensemblPattern))){
        message("Species is either not selected or not equal to 'mmu' or 'human'.
Please, select among default species or use options genomeAnnot and 
ensemblPattern. See example in the help page.")
        return(NULL)
    }else if(!missing(genomeAnnot) & !missing(ensemblPattern)){
        species = "manual"
    }else if(species == "mmu"){
        suppressMessages(library(org.Mm.eg.db, warn.conflicts = F))
        genomeAnnot <- org.Mm.eg.db
        ensemblPattern <- "ENSMUSG"
    }else if(species == "human"){
        suppressMessages(library(org.Hs.eg.db, warn.conflicts = F))
        genomeAnnot <- org.Hs.eg.db
        ensemblPattern <- "ENSG"
    }
    ensemblGenes <- rownames(countMatrix)[grep(ensemblPattern, 
                                               rownames(countMatrix))]
    symbolGenes <- rownames(countMatrix)[!grepl(ensemblPattern, 
                                                rownames(countMatrix))]
    message(paste0("Annotating ",length(ensemblGenes), " genes containing ", 
                   ensemblPattern, " pattern."))
    # if none of ENSEMBL genes are in database we fill their rowData rows with NA
    if(length(intersect(AnnotationDbi::keys(genomeAnnot, 
                                            keytype = "ENSEMBL"), 
                        ensemblGenes)) == 0){
        rowdataEnsembl <- data.frame(ENSEMBL = ensemblGenes,
                                     SYMBOL = NA,
                                     GENENAME = NA)
    }else{
        rowdataEnsembl <- AnnotationDbi::select(genomeAnnot, keys=ensemblGenes, 
                                                keytype="ENSEMBL", 
                                                columns=c("SYMBOL", 
                                                          "GENENAME"), 
                                                multiVals="first")
        rowdataEnsembl <- rowdataEnsembl[!duplicated(rowdataEnsembl$ENSEMBL),]
    }
    rowdataEnsembl$nameInCountMatrix <- ensemblGenes
    message("Annotating rest ", length(symbolGenes), " genes 
considering them as SYMBOLs.")
    rowdataSymbol <- AnnotationDbi::select(genomeAnnot, keys=symbolGenes, 
                                           keytype="SYMBOL", 
                                           columns=c("ENSEMBL", 
                                                     "GENENAME"), 
                                           multiVals="first")
    rowdataSymbol <- rowdataSymbol[!duplicated(rowdataSymbol$SYMBOL),]
    rowdataSymbol$nameInCountMatrix <- symbolGenes
    rowdata <- base::rbind(rowdataSymbol, rowdataEnsembl)
    rm(rowdataEnsembl, rowdataSymbol)
    
    # sometimes several ensembls give one symbol
    (mult_sym <- rowdata$SYMBOL[!is.na(rowdata$SYMBOL) & 
                                    duplicated(rowdata$SYMBOL)])
    
    # we decided don't combine such ensembls, but leave them unique with
    #"symbol_ensembl" ###
    (rowdata$SYMBOL[rowdata$SYMBOL %in% mult_sym] <- 
            paste(rowdata$SYMBOL[rowdata$SYMBOL %in% mult_sym],
                  rowdata$ENSEMBL[rowdata$SYMBOL %in% mult_sym],
                  sep = "_"))
    rm(mult_sym)
    
    ensembl <- useMart(biomart = "ensembl", dataset="mmusculus_gene_ensembl")
    message("Retrieving information about genes from biomaRt. 
            It can take up to five minutes, depends on Internet connection.")
    res <- getBM(attributes=c("ensembl_gene_id", "go_id", "name_1006",
                              "chromosome_name", "gene_biotype"), 
                 mart=ensembl)
    tmp <- res[!duplicated(res$ensembl_gene_id),]
    rowdata <- merge(rowdata, tmp[c("ensembl_gene_id",
                                    "chromosome_name", "gene_biotype")],
                     by.x = "ENSEMBL", by.y = "ensembl_gene_id",
                     all.x = TRUE, all.y = FALSE, sort = FALSE)
    rowdata_GO <- merge(rowdata, res[c("ensembl_gene_id",
                                       "go_id", "name_1006")],
                        by.x = "ENSEMBL", by.y = "ensembl_gene_id",
                        all.x = TRUE, all.y = FALSE, sort = FALSE)
    rowdata_GO <- rowdata_GO[!is.na(rowdata_GO$name_1006) &
                                 ((rowdata_GO$name_1006 == "cell surface") |
            (rowdata_GO$name_1006=="cell surface receptor signaling pathway")),]
    rowdata_GO$name_1006[duplicated(rowdata_GO$ENSEMBL)] <- 
        "cell surface receptor signaling pathway"
    rowdata_GO <- rowdata_GO[!duplicated(rowdata_GO$ENSEMBL),]
    rowdata <- merge(rowdata, rowdata_GO[c("nameInCountMatrix", 
                                           "go_id", "name_1006")],
                     by.x = "nameInCountMatrix", by.y = "nameInCountMatrix",
                     all.x = TRUE, all.y = TRUE, sort = FALSE)
    rm(tmp, ensembl, res, rowdata_GO)
    
    if(!is.null(rowData)){
        rowData$nameInCountMatrix <- rownames(rowData)
        rowdata <- merge(rowData, rowdata,
                         by.x = "nameInCountMatrix", by.y = "nameInCountMatrix",
                         all.x = TRUE, all.y = TRUE, sort = FALSE)
    }
    
    rownames(rowdata) <- rowdata$nameInCountMatrix
    rowdata <- rowdata[rownames(countMatrix),]
    stopifnot(all(rownames(rowdata) == rownames(countMatrix)))
    
    return(rowdata)
}

filterCells <- function(countMatrix, colData, genesSumThr = 100,
                        MoreBetter = c("genesNum", "sumCodPer", "genesSum"),
                        MoreWorse = c("sumMtPer")){
    message("Running filterCells.")
    countMatrix <- countMatrix[,colSums(countMatrix) > genesSumThr]
    colData <- colData[colnames(countMatrix),]
    mb <- MoreBetter
    mw <- MoreWorse
    
    reportTable <- data.frame(matrix(NA, ncol = length(mb)+length(mw), 
                                     nrow = nrow(colData)))
    colnames(reportTable) <- c(mb, mw)
    reportTable <- cbind(cellName = colData$cellName, reportTable)
    rownames(reportTable) <- reportTable$cellName
    
    stopifnot(all(colData$cellName==reportTable$cellName))
    
    for(j in 1:length(mb)){
        quan <- quantile(colData[,colnames(colData) == mb[j]])
        threshold <- 2.5*quan[2] - 1.5*quan[4]
        if(threshold < 0){
            threshold <- (quan[1]+quan[2]) / 2
        }
        reportTable[colData[,colnames(colData)==mb[j]] >= as.numeric(threshold),
                    colnames(reportTable)==mb[j]] = 1
        reportTable[colData[,colnames(colData)==mb[j]] < as.numeric(threshold),
                    colnames(reportTable)==mb[j]] = 0
    }
    for(j in 1:length(mw)){
        quan <- quantile(colData[,colnames(colData) == mw[j]])
        threshold <- 2.5*quan[4] - 1.5*quan[2]
        if(threshold > quan[5]){
            threshold <- (quan[3]+quan[4]) / 2
        }
        reportTable[colData[,colnames(colData)==mw[j]] <= as.numeric(threshold),
                    colnames(reportTable)==mw[j]] = 1
        reportTable[colData[,colnames(colData)==mw[j]] > as.numeric(threshold),
                    colnames(reportTable)==mw[j]] = 0
    }
    
    ### add columns with filtering score and verdict ###
    reportTable <- dplyr::mutate(reportTable, score = NA)
    reportTable$score <- rowSums(reportTable[,colnames(reportTable) %in% 
                                                 c(mb,mw)])
    reportTable <- dplyr::mutate(reportTable, filterPassed = NA)
    reportTable$filterPassed[reportTable$score >= length(mb)+length(mw)] <- 1
    reportTable$filterPassed[reportTable$score < length(mb)+length(mw)] <- 0
    
    ### add filtering verdict to colData ###
    colData <- dplyr::mutate(colData, filterPassed = NA)
    colData$filterPassed[colData$cellName %in% 
                    reportTable$cellName[reportTable$filterPassed == 1]] <- 1
    colData$filterPassed[colData$cellName %in% 
                    reportTable$cellName[reportTable$filterPassed == 0]] <- 0
    
    reportTable <- reportTable[order(reportTable$score, decreasing  =  FALSE), ]
    
    rm(threshold, j, mb, mw)
    
    rownames(colData) <- colData$cellName
    colData <- colData[colnames(countMatrix), ]
    stopifnot(all(rownames(colData) == colnames(countMatrix)))
    
    countMatrix <- countMatrix[,colData$filterPassed == 1]
    colData <- colData[colData$filterPassed == 1,]
    
    return(list(countMatrix, colData))
}

# from 2_mk_coldata_light.R
# I rewrote it with the new style
# This function creates colData or add columns mtGenes, genesNum, 
# codGenes, genesSum, codSum, mtPer, codPer, sumMtPer, sumCodPer to the 
# existing colData.
addCellsInfo <- function(countMatrix, rowData, colData = NULL){
    message("Adding cell info for cells filtering.")
    coldata <- data.frame(cellName = colnames(countMatrix),
                          stringsAsFactors = FALSE)
    
    ### add info about all genes in a cell ###
    coldata <- dplyr::mutate(coldata, genesNum = NA, genesSum = NA, oneUMI = NA)
    coldata$genesSum <- colSums(countMatrix)
    for(i in 1:ncol(countMatrix)){
        vec <- countMatrix[,i]
        coldata$genesNum[coldata$cellName == colnames(countMatrix)[i]] <- 
            length(vec[vec > 0])
        coldata$oneUMI[coldata$cellName == colnames(countMatrix)[i]] <- 
            length(vec[vec == 1])
    }
    rm(vec)
    coldata <- dplyr::mutate(coldata, 
                            oneUMIper = 100 * coldata$oneUMI / coldata$genesNum)
    
    ### add info about mitochondrial and protein-coding genes ###
    coldata <- dplyr::mutate(coldata, 
                             mtGenes = NA, mtSum = NA,
                             codGenes = NA, codSum = NA)
    for(i in 1:ncol(countMatrix)){
        mt <- countMatrix[rownames(countMatrix) %in% 
                            rowData$nameInCountMatrix[rowData$chromosome_name == 
                                                                        "MT"],i]
        coldata$mtGenes[coldata$cellName == colnames(countMatrix)[i]] <- 
            length(mt[mt > 0])
        coldata$mtSum[coldata$cellName == colnames(countMatrix)[i]]   <- sum(mt)
        
        cod <- countMatrix[rownames(countMatrix) %in% 
                    rowData$nameInCountMatrix[rowData$gene_biotype == 
                                                            "protein_coding"],i]
        coldata$codGenes[coldata$cellName == colnames(countMatrix)[i]] <- 
            length(cod[cod > 0])
        coldata$codSum[coldata$cellName == colnames(countMatrix)[i]] <- sum(cod)
    }
    
    coldata <- dplyr::mutate(coldata, 
                    mtPer      = 100 * coldata$mtGenes  / coldata$genesNum,
                    codPer     = 100 * coldata$codGenes / coldata$genesNum,
                    sumMtPer  = 100 * coldata$mtSum    / coldata$genesSum,
                    sumCodPer = 100 * coldata$codSum   / coldata$genesSum)
    rm(mt, cod)
    
    if(!is.null(colData)){
        colData$cellName <- rownames(colData)
        coldata <- merge(colData, coldata,
                         by.x = "cellName", by.y = "cellName",
                         all.x = FALSE, all.y = TRUE, sort = FALSE)
    }
    
    rownames(coldata) <- coldata$cellName
    coldata <- coldata[colnames(countMatrix),]
    stopifnot(all(rownames(coldata) == colnames(countMatrix)))
    
    return(coldata)
}

filterGenes <- function(countMatrix, rowData){
    # filter genes which are more than in 10 cells and less than (all-10) cells
    selRows <- ((rowSums(countMatrix[,] >= 1)) > 10)
    countMatrix <- countMatrix[selRows,]
    rowData <- rowData[rowData$nameInCountMatrix %in% rownames(countMatrix),]
    
    return(list(countMatrix, rowData))
}

normaliseCountMatrix <- function(countMatrix, 
                                 species,
                                 method="default", 
                                 sizes=c(20,40,60,80,100),
                                 rowData=NULL, 
                                 colData=NULL,
                                 alreadyCellFiltered = FALSE,
                                 runQuickCluster = TRUE){
    # Does normalisation of count matrix with.
    # There are 2 possible methods: "default" or "census"
    # The function returns SCE object with normalised count matrix
    if(method == "default"){
        rowData <- annotateGenes(countMatrix, species = species,
                                 rowData = rowData)
        colData <- addCellsInfo(countMatrix, rowData = rowData,
                                colData = colData)
        if(!alreadyCellFiltered){
            filterCellsResult <- filterCells(countMatrix, colData)
            countMatrix <- filterCellsResult[[1]]
            colData <- filterCellsResult[[2]]
            rm(filterCellsResult)
        }
        filterGenesResult <- filterGenes(countMatrix, rowData)
        countMatrix <- filterGenesResult[[1]]
        rowData <- filterGenesResult[[2]]
        rm(filterGenesResult)
        
        stopifnot(all(rownames(countMatrix) == rownames(rowData)))
        stopifnot(all(colnames(countMatrix) == rownames(colData)))
        
        sce <- 
            SingleCellExperiment(assays = list(counts = as.matrix(countMatrix)), 
                                 colData=colData, rowData=rowData)
        
        # normalization
        message("Running normalization. It can take a while, depends on the 
number of cells.")
        if(runQuickCluster){
            cl <- tryCatch(scran::quickCluster(sce), error=function(e) NULL)
        }else{
            cl <- NULL
        }
        
        # compute sizeFactors which will be used for normalization
        sceNorm <- scran::computeSumFactors(sce, sizes = sizes, clusters = cl)
        
        message("summary(sizeFactors(sceObject)):")
        print(summary(sizeFactors(sceNorm)))
        if(length(sizeFactors(sceNorm)[sizeFactors(sceNorm) <= 0]) > 0){
            message("Cells with negative sizeFactors will be deleted before the 
downstream analysis.")
        }
        sceNorm <- sceNorm[, sizeFactors(sceNorm) > 0]
        sceNorm <- scater::normalize(sceNorm)
        rm(sce)
        
        return(sceNorm)
        
    }else if(method == "census"){
        sceObject <- normalize_dataset(as.matrix(countMatrix))
        colData(sceObject)$cellName = rownames(colData(sceObject))
        return(sceObject)
    }else{
        print("Wrong method. Unmodified count matrix returned.")
        return(countMatrix)
    }
}


#Please, note that this dataset is provided only for understanding the concept of an optimal normalization
#and does not include any regressions of batch-corrections
#However, for an average dataset this is more than enough to get meaningful result
#you'll find similar or simpler approach in a vest majority of the public tutorials

#Try to adjust different values and understand filters impact on final expression
#Code to create QC plots you'll find here
#https://hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#normalization-practice-reads

#Test dataset does not contain mt genes and spikeIns, but has a lot of rear cell types
#Try to find filtering thresholds to eliminate some rear cell populations
#Take any arbitrary dataset form this database https://hemberg-lab.github.io/scRNA.seq.datasets/
#and have your final cleaned and normalized expression matrix for this dataset ready on the next Wednesday

normalize_dataset <- function(mtx)
{
    filter_na_inf <- function(mtx)
    {
        colsS = Matrix::colSums(mtx)
        rowsS = Matrix::rowSums(mtx)
        
        goodCols = !is.infinite(colsS)&!is.na(colsS)
        goodRows = !is.infinite(rowsS)&!is.na(rowsS)
        
        percC = (table(goodCols)/ncol(mtx))["TRUE"]
        percR = (table(goodRows)/nrow(mtx))["TRUE"]
        
        if(is.na(percR))
        {
            mtx = mtx[,goodCols]
        }else if(is.na(percC))
        {
            mtx = mtx[goodRows,]
        }else if(percC>percR)
        {
            mtx = mtx[,goodCols]
        }else{
            mtx = mtx[goodRows,]
        }
        return(mtx)
    }
    
    get_types <- function(arr)
    {
        return(sapply(arr, function(cell) { return(unlist(strsplit(cell, ".", 
                                                            fixed = T))[1])}))
    }
    
    #rds contains pre-created Single Cell Dataset
    #lets consider that you don`t have a "counts" slot, pre-normalized "exprs" only
    #download dataset: https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/deng-reads.rds
    
    
    mtx = filter_na_inf(mtx)
    
    #census de-normalization expression to read-counts
    #rel2abs article: https://www.ncbi.nlm.nih.gov/pubmed/28114287
    fd = as.matrix(rownames(mtx))
    rownames(fd) = fd
    colnames(fd)[1]<-"gene_short_name"
    
    pd = as.matrix(colnames(mtx))
    rownames(pd) = pd
    colnames(pd)[1]<-"cell_name"
    
    pd = new("AnnotatedDataFrame", data = as.data.frame(pd))
    fd = new("AnnotatedDataFrame", data = as.data.frame(fd))
    
    
    relative = newCellDataSet(mtx,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.1,
                              expressionFamily = tobit(Lower = 0.1))
    
    rpc_matrix = relative2abs(relative, t_estimate = estimate_t(exprs(relative)), 
                              method = "num_genes", cores = 8)
    #filter na/inf
    #census normalization often ends with  some percentage of NaNs for some cells.
    #keep in mind, that sometimes it can "NaNify"" away up to half of the dataset
    rpc_matrix = filter_na_inf(rpc_matrix)
    
    
    {
        #NOTE! This step is outside of the optimal pipeline and logic!
        #if you've lost too much cells, you can try to normalize lost part separately
        #and then merge them and then normalize using "normalization through pooling"
        #or perform a batch-effect regression? marking two datasets as batches
        #but this is just a dirty hack, never tell anyone that you used it
        
        #you can use it only if your "lost cells" seem to be evenly distributed 
        #in the low-dimensional representation of the dataset
        #optimal way to visually assess it is to use 
        #PCA(no less than 50 components) over TSNE (will tell about them in the next lecture)
        cells_res = setdiff(colnames(mtx), colnames(rpc_matrix))
        #I will not provide you an implementation, but you may try do it yourself
    }
    
    sce = SingleCellExperiment(assays = list(counts = rpc_matrix), 
                               colData = get_types(colnames(rpc_matrix)))
    rowData(sce)$feature_symbol = rownames(sce)
    sce = sce[!duplicated(rowData(sce)$feature_symbol), ]
    
    #Filter infrequent cells and genes
    lowerDetectionLimit_cell = 2
    lowerDetectionLimit_gene = 2
    numcells_sh = 2
    numgenes_sh = 2
    if(package.version("scater") == "1.8.4"){
        numcells = nexprs(sce, detection_limit = lowerDetectionLimit_cell, 
                          byrow = T)
        numgenes = nexprs(sce, detection_limit = lowerDetectionLimit_gene, 
                          byrow = F)
    }else{
        numcells = nexprs(sce, lowerDetectionLimit = lowerDetectionLimit_cell, 
                          byrow = T)
        numgenes = nexprs(sce, lowerDetectionLimit = lowerDetectionLimit_gene, 
                          byrow = F)
    }
    keep.gene = numcells >= numcells_sh
    keep.cell = numgenes >= numgenes_sh
    
    
    #cat("genes_left:", round(length(which(keep.gene))/nrow(sce)*100, 2), "%")
    #cat("cells_left:", round(length(which(keep.cell))/ncol(sce)*100, 2), "%")
    sce = sce[keep.gene, keep.cell]
    
    #Filter genes that has flat expression profile
    #Outside of the well-known best practices, good from my experiense
    gene_levels_sh = 3
    gene_levels_unique = apply(counts(sce), 1, FUN = unique)
    gene_levels_lengths = unlist(lapply(gene_levels_unique, length))
    genes_good = names(which(gene_levels_lengths>=gene_levels_sh))
    
    #cat("genes_left:", round(length(genes_good)/nrow(sce)*100, 2), "%")
    sce = sce[genes_good,]
    
    
    #MT QC
    # You can use databases to detect gene symbol annotation and automate gene selection step
    
    # detected_genome = detect_genome(rownames(sce))
    # if(is.null(detected_genome))
    # {
    #   return(NULL)
    # }
    # anno = get_genome_annotation_names_mapping(detected_genome$org)
    # #Filter unknown genes
    # is.mito =  rownames(sce) %in% anno$gene_symbol[which(anno$chr  ==  "MT")]
    
    #however here we are using simpler approach
    is.mito_offline = (grepl("^mt-", rownames(sce)) | grepl("^MT-", rownames(sce)))
    
    
    # sce = calculateQCMetrics(sce, exprs_values = "counts", 
    #                          feature_controls = list(mt = rownames(sce)[which(is.mito_offline)], 
    #                                                   ercc = rownames(sce)[which(isSpike(sce, "ercc"))]),
    #                          cell_controls = NULL, nmads = 3, pct_feature_controls_threshold = 80)
    
    if(package.version("scater") == "1.8.4"){
        sce = calculateQCMetrics(sce, exprs_values = "counts", 
                                 feature_controls = list(mt = 
                                        rownames(sce)[which(is.mito_offline)]), 
                                 cell_controls = NULL)
    }else{
        sce = calculateQCMetrics(sce, exprs_values = "counts", 
                                 feature_controls = list(mt = 
                                        rownames(sce)[which(is.mito_offline)]), 
                                 cell_controls = NULL, nmads = 3, 
                                 pct_feature_controls_threshold = 80)
    }
    
    final_drop = rep(F, ncol(sce))
    
    libsize.drop = isOutlier(sce$total_counts, nmads = 3, type = "low", log = T)
    final_drop = libsize.drop|final_drop
    
    feature.drop = isOutlier(sce$total_features, nmads = 3, type = "low", log = T)
    final_drop = feature.drop|final_drop
    
    mito.drop = isOutlier(sce$total_counts_mt, nmads = 3, type = "high", log = F)
    final_drop = mito.drop|final_drop
    
    # spike.drop = isOutlier(sce$total_counts_ercc, nmads = 3, type = "both", log = F)
    # final_drop = spike.drop|final_drop
    
    #cat("cells_left:", round((ncol(sce)-length(which(final_drop)))/ncol(sce)*100, 2), "%")
    sce = sce[,!final_drop]
    
    
    
    #Calculating size-factors for cells normalization based on our spike-ins separately from all genes
    # sce = computeSpikeFactors(sce, general.use = FALSE)
    
    #Calculating main size-factors by cells pooling
    #article about pooling normalization https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4848819/
    #I use to use formula like this to calculate norm levels
    norm_sizes = seq(min(20, max(2, round(ncol(sce)/4))), min(max(2, round(ncol(sce)/2)), 50), 1)
    # norm_sizes = seq(20, 100, 5)
    #but you can use default values for now
    sce = computeSumFactors(sce, sizes = norm_sizes, positive = F)
    
    #Compute normalized expression values using the size factors stored in the object
    sce <- normalise(sce)
    return(sce)
}

getGenesInfo <- function(genes, databaseDir, groupBy = "clusters",
                         orderGenes = "initial",
                         silent = FALSE, coresGenes = 20){
    
    # MGI
    getMGIentry <- function(MGIid){
        library('rvest')
        library(S4Vectors)
        data <- NULL
        if(!is.na(MGIid)){
            url <- paste0("http://www.informatics.jax.org/marker/",
                          MGIid)
            webpage <- read_html(url)
            data_html <- html_nodes(webpage,'#mpMarkerClip')
            data <- html_text(data_html)
            data <- gsub("\t", "", data, perl = TRUE)
            data <- gsub("\n", "", data, perl = TRUE)
            data <- gsub(";", ",", data)
            
            rm(url, data_html, webpage)
        }
        if(!S4Vectors::isEmpty(data)){
            return(data)
        }else{
            return(NA)
        }
    }
    
    #NCBI
    getNCBIentry <- function(NCBIid){
        library('rvest')
        library(S4Vectors)
        data <- NULL
        if(!is.na(NCBIid)){
            url <- paste0("https://www.ncbi.nlm.nih.gov/gene?cmd=Retrieve&dopt=full_report&list_uids=",
                          NCBIid)
            webpage <- read_html(url)
            data_html <- html_nodes(webpage,'#summaryDl dd:nth-child(20)')
            data <- html_text(data_html)
            data <- gsub(" See more", "", data)
            data <- gsub("\n          human\n          all\n", "", data)
            data <- gsub(";", ",", data)
            
            rm(url, data_html, webpage)
        }
        if(!S4Vectors::isEmpty(data)){
            return(data)
        }else{
            return(NA)
        }
    }
    
    # Uniprot
    getUniprotEntry <- function(UniprotID){
        library('rvest')
        library(S4Vectors)
        data <- NULL
        if(!is.na(UniprotID)){
            url <- paste0("https://www.uniprot.org/uniprot/",
                          UniprotID, "#function")
            webpage <- read_html(url)
            data_html <- html_nodes(webpage,'h2+ .annotation > span')
            data <- html_text(data_html)
            data <- gsub("* Publication.*", " Publications", data)
            data <- gsub('<p><a href="/manual/evidences#ECO:0000250">More...</a></p>', 
                         "", data)
            data <- gsub('\r\n .* <p>Manually curated information .* toiUniProtKB:.* (.*_.*).', 
                         "", data)
            data <- gsub('\r\n .* <p>Manual validated information.*', 
                         "", data)
            data <- gsub('UniRule annotation', "", data)
            data <- gsub(";", ",", data)
            
            rm(url, data_html, webpage)
        }
        if(!S4Vectors::isEmpty(data)){
            return(data)
        }else{
            return(NA)
        }
    }
    
    library(DataCombine)
    library(doParallel)
    #database <- read.delim(file.path(databaseDir, "Mmus_gene_database.tsv"),
    #                       stringsAsFactors = FALSE)
    database <- read.delim(file.path(databaseDir, "Mmus_gene_database_secretedMol.tsv"),
                           stringsAsFactors = FALSE)
    
    genesEnsembl <- genes[grepl("ENSMUSG", genes$geneName),]
    genesSymbol <- genes[!grepl("ENSMUSG", genes$geneName),]
    
    resultEnsembl <- merge(genesEnsembl, database,
                           by.x = "geneName", by.y = "Ensembl",
                           all.x = TRUE, all.y = FALSE, sort = FALSE)
    resultEnsembl$Ensembl <- resultEnsembl$geneName
    
    resultSymbol <- merge(genesSymbol, database,
                          by.x = "geneName", by.y = "Symbol",
                          all.x = TRUE, all.y = FALSE, sort = FALSE)
    resultSymbol$Symbol <- resultSymbol$geneName
    
    colnames(resultEnsembl)
    colnamesOrder = c("geneName", "clusters", "Ensembl", "Symbol", "Name",            
                      "Feature.Type", "MGI.Gene.Marker.ID", "Entrez.Gene.ID", 
                      "Uniprot.ID", "chromosome_name", "go_id", "name_1006")
    resultEnsembl <- resultEnsembl[,colnamesOrder]
    resultSymbol <- resultSymbol[,colnamesOrder]
    
    result <- rbind(resultSymbol, resultEnsembl)
    #colnames(result)[11:12] <- c("CellSurface.GOid", "CellSurface.GOname")
    
    rownames(result) <- result$geneName
    result <- result[genes$geneName,]
    
    if(orderGenes == "alphabetical"){
        result <- result[order(result$geneName),]
    }
    
    myCluster <- makeCluster(coresGenes, # number of cores to use
                             type = "PSOCK") # type of cluster
    registerDoParallel(myCluster)
    
    # MGI vector
    if(!silent){
        message("Collecting knockout phenotype information from MGI.")
    }
    
    MGI <- unname(unlist(foreach(MGIid=result$MGI.Gene.Marker.ID) %dopar% getMGIentry(MGIid)))
    
    #NCBI vector
    if(!silent){
        message("Retrieving info from NCBI.")
    }
    
    NCBI <- unname(unlist( foreach(NCBIid=result$Entrez.Gene.ID) %dopar% getNCBIentry(NCBIid) ))
    
    # Uniprot
    if(!silent){
        message("Getting summary from Uniprot.")
    }
    
    UniprotFunction <- unname(unlist( foreach(UniprotID=result$Uniprot.ID) %dopar% getUniprotEntry(UniprotID) ))
    
    stopCluster(myCluster)
    
    #MGI <- data.frame(MGI, stringsAsFactors = FALSE)
    #NCBI <- data.frame(NCBI, stringsAsFactors = FALSE)
    #UniprotFunction <- data.frame(UniprotFunction, stringsAsFactors = FALSE)
    
    result <- cbind(result, MGI, NCBI, UniprotFunction)
    
    result$MGI <- as.character(result$MGI)
    result$NCBI <- as.character(result$NCBI)
    result$UniprotFunction <- as.character(result$UniprotFunction)
    
    rownames(result) <- c(1:nrow(result))
    
    # inserting space for comments
    if(any(colnames(result) %in% groupBy) & 
       (orderGenes == "initial") &
       (length(unique(result$clusters)) > 1)){
        resultFinal <- result
        groupingTable <- table(resultFinal[,groupBy])
        groupingTable <- groupingTable[unique(resultFinal$clusters)]
        resultFinal <- InsertRow(resultFinal, c("For notes:", 
                                                rep("", (ncol(result))-1)), 
                                 RowNum = 1)
        
        RowNum <- groupingTable[1] + 1
        for(i in 1:(length(groupingTable)-1)){
            resultFinal <- InsertRow(resultFinal, rep("", ncol(result)), 
                                     RowNum = (RowNum + 1))
            resultFinal <- InsertRow(resultFinal, c("For notes:", 
                                                    rep("", (ncol(result))-1)), 
                                     RowNum = (RowNum + 2)) 
            RowNum <- RowNum + 2 + groupingTable[i+1]
        }
        result <- resultFinal
        rm(resultFinal)
    }
    
    rm(resultEnsembl, resultSymbol, database, colnamesOrder)
    
    #result <- result[,c("geneName", "clusters", "Name",             
    #                    "Feature.Type",  "CellSurface.GOid",
    #                    "CellSurface.GOname", "MGI", "NCBI",    
    #                    "UniprotFunction", "chromosome_name", 
    #                    "Symbol", "Ensembl", "MGI.Gene.Marker.ID", 
    #                    "Entrez.Gene.ID", "Uniprot.ID")]
    result <- result[,c("geneName", "clusters", "Name",             
                        "Feature.Type",  "go_id",
                        "name_1006", "MGI", "NCBI",    
                        "UniprotFunction", "chromosome_name", 
                        "Symbol", "Ensembl", "MGI.Gene.Marker.ID", 
                        "Entrez.Gene.ID", "Uniprot.ID")]
    return(result)
}

saveGenesInfo <- function(inputDir, pattern, outputDir, databaseDir, 
                          sep = ";", header = TRUE, 
                          startFromFile = 1, 
                          groupBy = "clusters", # getGenesInfo params
                          orderGenes = "initial",
                          silent = TRUE, coresGenes = 20){
    filesList <- list.files(inputDir, pattern = pattern)
    filePrefix <- do.call(rbind, strsplit(filesList, "[.]"))[,1]
    
    for(file in startFromFile:length(filesList)) {
        cat(file, "\n")
        genes <- read.delim(file.path(inputDir, filesList[file]), 
                            sep = ";", header = header,
                            stringsAsFactors = FALSE)
        
        result <- getGenesInfo(genes, databaseDir, 
                               groupBy = "clusters", 
                               orderGenes = "initial",
                               silent = FALSE, coresGenes = 20)
        
        message("Writing the output file number ", file, "\n")
        
        write.table(result, file = file.path(outputDir, 
                                             paste0(filePrefix[file], 
                                                    "_genesInfo.csv")),
                    quote = FALSE, sep = ";", row.names = FALSE)
    }
}

saveMarkersLists <- function(experimentName, dataDirectory, 
                             inputDir = file.path(dataDirectory, "marker_genes"), 
                             outputDir = file.path(dataDirectory, 
                                                   paste0("marker_genes/markers_lists")), 
                             pattern = "genes.tsv", Ntop = 100){
    dir.create(outputDir, showWarnings=F)
    fnames <- list.files(inputDir, pattern = pattern)
    for(i in 1:length(fnames)){
        tmp <- read.delim(file.path(inputDir, fnames[i]), stringsAsFactors = FALSE)
        markerList <- as.data.frame(tmp$Gene[1:Ntop])
        outputName <- gsub(paste0("_", pattern), "", fnames[i])
        clusterName <- gsub(paste0(experimentName, "_"), "", outputName)
        colnames(markerList) <- "geneName"
        markerList$clusters <- clusterName
        write.table(markerList, file.path(outputDir, 
                                          paste0(outputName, "_markers.csv")), 
                    row.names = FALSE, quote = FALSE, sep = ";")
    }
}