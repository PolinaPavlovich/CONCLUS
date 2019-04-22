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
suppressMessages(library(xlsx, warn.conflicts = F))
suppressMessages(library(grDevices, warn.conflicts = F))
suppressMessages(library(S4Vectors, warn.conflicts = F))
suppressMessages(library(Biobase, warn.conflicts = F))
suppressMessages(library(foreach, warn.conflicts = F))

### Internal function.
### this function calculates PCA and then tSNE with PCs and perplexities ###
### it returns a list of pSNE = PCA+tSNE results ###
### to get XY coordinates call psne_res[1,i][[1]] ###
### i = [1:length(PCs)*length(perplexities)] is a number of iteration ###

getTSNEresults <- function(expressionMatrix, cores=1,
                           PCs=c(4, 6, 8, 10, 20, 40, 50),
                           perplexities=c(30, 40), randomSeed=42){
    PCAData <- prcomp(t(expressionMatrix))$x
    myCluster <- parallel::makeCluster(cores, # number of cores to use
                             type = "PSOCK") # type of cluster
    doParallel::registerDoParallel(myCluster)
    tSNECoordinates <- foreach::foreach(PCA=rep(PCs, length(perplexities)),
                               perp=rep(perplexities, each=length(PCs)),
                               .combine='cbind') %dopar% {
                                   library(SingleCellExperiment)
                        scater::plotTSNE(SingleCellExperiment(assays=list(
                            logcounts=t(PCAData[,1:PCA]))),
                        scale_features=FALSE, perplexity=perp,
                        rand_seed=randomSeed, theme_size=13, return_SCESet=FALSE)
                        }
    parallel::stopCluster(myCluster)
    message(paste("Calculated", length(PCs)*length(perplexities),
"2D-tSNE plots."))
    return(tSNECoordinates)
}

# Do not export this function.
.testClustering <- function(sceObject, dataDirectory, experimentName,
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
  tSNE <- getTSNEresults(expr = Biobase::exprs(sceObject), cores=1,
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

#' To check one iteration of clustering before running full workflow CONCLUS.
#' 
#' This function generates a single clustering iteration of CONCLUS to check whether
#' chosen parameters for dbscan are suitable for your data.
#'
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @param dbscanEpsilon a parameter of fpc::dbscan() function.
#' @param minPts a parameter of fpc::dbscan() function.
#' @param PCs a vector of PCs for plotting.
#' @param perplexities vector of perplexities (t-SNE parameter).
#' @param randomSeed random seed for reproducibility.
#' @param width plot width.
#' @param height plot height.
#' @param ... other pdf() arguments.
#'
#' @return t-SNE results, a distance graph plot, a t-SNE plot colored by test clustering solution.
#' @export
testClustering <- function(sceObject, dataDirectory, experimentName,
                           dbscanEpsilon=1.4,
                           minPts=5,
                           perplexities = c(30), PCs = c(4),
                           randomSeed = 42,
                           width=7, height=7, ...){

  .testClustering(sceObject, dataDirectory, experimentName,
                  dbscanEpsilon=dbscanEpsilon,
                  minPts=minPts,
                  perplexities = perplexities, PCs = PCs,
                  randomSeed = randomSeed,
                  width=width, height=height, ...)
}

#' Choose palette for a plot.
#'
#' It is an internal function usually applied for choosing the palette for clusters.
#' Depending if the number of clusters is more than 12 or not, one of two built-in palettes will be applied.
#' If you give your vector of colors, the function will not change them.
#' If the number of clusters is more than 26, it will copy colors to get the needed length of the palette.
#'
#' @param colorPalette Either "default" or a vector of colors, for example c("yellow", "#CC79A7"). 
#' @param clustersNumber number of clusters in the output palette.
#'
#' @return Color palette with the number of colors equal to the clusterNumber parameter.
#' @export
choosePalette <- function(colorPalette, clustersNumber){

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

#' Run CONCLUS in one click
#'
#' This function performs core CONCLUS workflow. It generates PCA and t-SNE coordinates, 
#' runs DBSCAN, calculates similarity matrices of cells and clusters, assigns cells to clusters,
#' searches for positive markers for each cluster. The function saves plots and tables into dataDirectory.
#'
#' @param sceObject a SingleCellExperiment object with your data.
#' @param dataDirectory CONCLUS will create this directory if it doesn't exist and store there all output files.
#' @param experimentName most of output file names of CONCLUS are hardcoded.
#' experimentName will stay at the beginning of each output file name to
#' distinguish different runs easily.
#' @param colorPalette a vector of colors for clusters.
#' @param statePalette a vector of colors for states.
#' @param clusteringMethod a clustering methods passed to hclust() function.
#' @param epsilon a parameter of fpc::dbscan() function.
#' @param minPoints a parameter of fpc::dbscan() function.
#' @param k preferred number of clusters. Alternative to deepSplit. A parameter of cutree() function.
#' @param PCs a vector of first principal components.
#' For example, to take ranges 1:5 and 1:10 write c(5, 10).
#' @param perplexities a vector of perplexity for t-SNE.
#' @param randomSeed random seed for reproducibility.
#' @param deepSplit intuitive level of clustering depth. Options are 1, 2, 3, 4.
#' @param preClustered if TRUE, it will not change the column clusters after the run.
#' However, it will anyway run DBSCAN to calculate similarity matrices.
#' @param orderClusters can be either FALSE (default) of "name".
#' If "name", clusters in the similarity matrix of cells will be ordered by name.
#' @param cores maximum number of jobs that CONCLUS can run in parallel.
#' @param plotPDFcellSim if FALSE, the similarity matrix of cells will be saved in png format.
#' FALSE is recommended for count matrices with more than 2500 cells due to large pdf file size.
#' @param deleteOutliers whether cells which were often defined as outliers by dbscan must be deleted.
#' It will require recalculating of the similarity matrix of cells. Default is FALSE.
#' Usually those cells form a separate "outlier" cluster and can be easier distinguished and deleted later
#' if necessary.
#' @param tSNEalreadyGenerated if you already ran CONCLUS ones and have t-SNE coordinated saved
#' You can set TRUE to run the function faster since it will skip the generation of t-SNE coordinates and use the stored ones. 
#' Option TRUE requires t-SNE coordinates to be located in your 'dataDirectory/tsnes' directory.
#' @param tSNEresExp experimentName of t-SNE coordinates which you want to use.
#' This argument allows copying and pasting t-SNE coordinates between different CONCLUS runs without renaming the files.
#'
#' @keywords CONCLUS
#' @export
#' @return A SingleCellExperiment object.

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
                       deleteOutliers = TRUE,
                       tSNEalreadyGenerated = FALSE,
                       tSNEresExp = ""){

  initialisePath(dataDirectory)

  # Generating 2D tSNE plots
  if(!tSNEalreadyGenerated){
      tSNEResults <- generateTSNECoordinates(sceObject, dataDirectory,
                                             experimentName, PCs=PCs,
                                             perplexities=perplexities,
                                             randomSeed = randomSeed)
  }else{
      tSNEResults <- readRDS(file.path(dataDirectory, "output_tables",
                                     paste0(tSNEresExp,"_tSNEResults.rds")))
  }

  if(preClustered){
    # Running dbscan
    message("Running dbscan using ", cores, " cores.")
    dbscanResults <- runDBSCAN(tSNEResults, sceObject, dataDirectory,
                               experimentName, epsilon=epsilon,
                               minPoints=minPoints,
                               cores=cores)

    # assigning cells to clusters
    message("Calculating cells similarity matrix.")
    cellsSimilarityMatrix <- clusterCellsInternal(dbscanResults, sceObject, clusterNumber=k,
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

  print(table(SummarizedExperiment::colData(sceObjectFiltered)$clusters,
              dnn=list("Cells distribution by clusters")))

  clustersNumber <- length(unique(SummarizedExperiment::colData(sceObjectFiltered)$clusters))
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
                    columnName = "clusters",
                    tSNEresExp = tSNEresExp)
  plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
                    PCs=PCs, perplexities=perplexities, colorPalette,
                    columnName = "noColor",
                    tSNEresExp = tSNEresExp)
  if(any(colnames(SummarizedExperiment::colData(sceObjectFiltered)) %in% "state")){
      plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
                        PCs=PCs, perplexities=perplexities, statePalette,
                        columnName = "state",
                        tSNEresExp = tSNEresExp)
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

#' Export matrix to a file.
#' 
#' The function allows you to export a matrix to a .csv file with a hard-coded filename (according to experimentName) 
#' in the "dataDirectory/output_tables" directory for further analysis.
#'
#' @param matrix your matrix (e.g., expression matrix)
#' @param dataDirectory CONCLUS output directory for a given experiment (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear at the beginning of the filenames 
#' (supposed to be the same for one experiment during the workflow).
#' @param name name of the file. Will be placed after the experimentName header.
#' @export
exportMatrix <- function(matrix, dataDirectory, experimentName, name){

  fileName <- paste0(experimentName, "_", name, ".csv")
  write.table(matrix, file=file.path(dataDirectory, "output_tables", fileName),
              sep = ",")
}

#' Create all needed directories for CONCLUS output.
#'
#' @param dataDirectory output directory for a given CONCLUS run (supposed to be the same for one experiment during the workflow).
#' @export
initialisePath <- function(dataDirectory){
  # creates directories for further writing of results.
  # names of directories are hardcoded.
  # no idea if it is good or bad.

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

#' Generate and save t-SNE coordinates with selected parameters.
#'
#' The function generates several t-SNE coordinates based on given perplexity and ranges of PCs. 
#' Final number of t-SNE plots is length(PCs)*length(perplexities)
#' It writes coordinates in "dataDirectory/tsnes" subfolder.
#'
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @param randomSeed random seed for reproducibility.
#' @param cores maximum number of jobs that CONCLUS can run in parallel.
#' @param PCs a vector of first principal components.
#' For example, to take ranges 1:5 and 1:10 write c(5, 10).
#' @param perplexities a vector of perplexity (t-SNE parameter).
#'
#' @return An object with t-SNE results (coordinates for each plot).
#' @export
generateTSNECoordinates <- function(sceObject, dataDirectory, experimentName,
                                    randomSeed=42, cores=14,
                                    PCs=c(4, 6, 8, 10, 20, 40, 50),
                                    perplexities=c(30,40)){

  tSNEDirectory <- "tsnes"
  message(paste0("Running TSNEs using ", cores, " cores."))
  TSNEres <- getTSNEresults(Biobase::exprs(sceObject), cores=cores, PCs=PCs,
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
  saveRDS(TSNEres, file = file.path(dataDirectory, "output_tables",
                             paste0(experimentName,"_tSNEResults.rds")))
  rm(tSNEDirectory, PCA, perp)
  return(TSNEres)
}

checkTSNEPicture <- function(tSNEResults, sceObject){
  # gives the unclustered picture of tSNE
  # to be sure that normalization step was
  # successful

  tSNECoords <- tSNEResults[[1]]
  tSNECoords <- tSNECoords[rownames(tSNECoords) %in%
                               SummarizedExperiment::colData(sceObject)$cellName, ]

  ggplot2::ggplot(tSNECoords, aes_string(x=names(tSNECoords)[1],
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
  factoextra::fviz_cluster(dbscanResults, tSNEData, ellipse=TRUE, geom="point",
               legend="bottom")
}

### This function calculates dbscan for all t-SNE from TSNEtables with all
### combinations of paramenters from epsilon and minPoints
### it does not set random seed. It allows to vary this parameter automatically
### it returns a matrix where columns are iterations
### number of iterations is equal to
### ncol is (TSNEtables)*length(epsilon)*length(epsilon)

mkDbscan <- function(TSNEtables, cores = 14, epsilon = c(1.2, 1.5, 1.8),
                    minPoints = c(15, 20)){
    myCluster <- parallel::makeCluster(cores, # number of cores to use
                             type = "PSOCK") # type of cluster
    doParallel::registerDoParallel(myCluster)
    dbscanResults <- foreach::foreach(i=rep(rep(1:ncol(TSNEtables),
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
    parallel::stopCluster(myCluster)
    return(dbscanResults)
}

#' Run clustering iterations with selected parameters using DBSCAN.
#'
#' This function returns a matrix of clustering iterations of DBSCAN.
#'
#' @param tSNEResults results of conclus::generateTSNECoordinates() function.
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames 
#' (supposed to be the same for one experiment during the workflow).
#' @param cores maximum number of jobs that CONCLUS can run in parallel.
#' @param epsilon a fpc::dbscan() parameter.
#' @param minPoints a fpc::dbscan() parameter.
#'
#' @return A matrix of DBSCAN results.
#' @export
runDBSCAN <- function(tSNEResults, sceObject, dataDirectory, experimentName,
                      cores=14, epsilon=c(1.3, 1.4, 1.5), minPoints=c(3, 4)){

  outputDataDirectory <- "output_tables"
  # taking only cells from the sceObject
  for(i in 1:ncol(tSNEResults)){
      tmp <- tSNEResults[1,i][[1]]
      tSNEResults[1,i][[1]] <- tmp[colnames(sceObject),]
  }
  dbscanResults <- mkDbscan(tSNEResults, cores = cores, epsilon = epsilon,
                             minPoints = minPoints)
  dbscanResults <- t(dbscanResults)
  colnames(dbscanResults) <- SummarizedExperiment::colData(sceObject)$cellName
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
  colData <- SummarizedExperiment::colData(sceObject)[SummarizedExperiment::colData(sceObject)$cellName
                                %in% colnames(dbscanMatrix),]

  if(is.vector(colData)){
    colData <- S4Vectors::DataFrame(cellName=colData, row.names=colData)
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
  SummarizedExperiment::colData(sceObject) <- colData

  message(paste(numberOfCellsBefore - numberOfCellsAfter,
              "outliers were excluded from the SingleCellExperiment object.\n"))

  return(list(sceObject, dbscanMatrix))
}

mkSimMat <- function(mat, cores=14){

    myCluster <- parallel::makeCluster(cores, # number of cores to use
                             type = "PSOCK") # type of cluster
    doParallel::registerDoParallel(myCluster)

    simMats <- foreach::foreach(i=1:nrow(mat)) %dopar% {
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
    parallel::stopCluster(myCluster)

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

#' Cluster cells and get similarity matrix of cells.
#' 
#' The function returns consensus clusters by using hierarchical clustering on the similarity matrix of cells.
#' It provides two options: to specify an exact number of clusters (with clusterNumber parameter)
#' or to select the depth of splitting (deepSplit parameter).
#' 
#' @param dbscanMatrix an output matrix of conclus::runDBSCAN() function.
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param clusterNumber a parameter, specifying the exact number of cluster.
#' @param deepSplit a parameter, specifying how deep we will split the clustering tree. It takes integers from 1 to 4.
#' @param cores maximum number of jobs that CONCLUS can run in parallel.
#' @param clusteringMethod a clustering methods passed to hclust() function.
#'
#' @return A SingleCellExperiment object with modified/created "clusters" column in the colData, and cells similarity matrix.
#' @export
clusterCellsInternal <- function(dbscanMatrix, sceObject, clusterNumber=0,
                         deepSplit, cores=14,
                         clusteringMethod = "ward.D2") {
  # 

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

  SummarizedExperiment::colData(sceObject)$clusters <- factor(clusters)

  return(list(sceObject, cellsSimilarityMatrix))
}

# Do not export this function.
.plotCellSimilarity <- function(sceObject, cellsSimilarityMatrix, dataDirectory,
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
                                   RColorBrewer::brewer.pal(n = 7,
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
                               widthPNG = 800, heightPNG = 750 #png
                               ){
  # plots cells correlation matrix gained form
  # clusterCellsInternal() function as the result of DBSCAN

  graphsDirectory <- "pictures"
  colData <- SummarizedExperiment::colData(sceObject)
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
      #message("Plot type is not pdf. Saving in svg.")
      #svg(file.path(dataDirectory, graphsDirectory,
      #              paste(experimentName,
      #                    "cells_correlation", clustersNumber,
      #                    "clusters.svg", sep="_")),width=8,height=8)
  }
  pheatmap::pheatmap(cellsSimilarityMatrix,
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
      return(pheatmap::pheatmap(cellsSimilarityMatrix,
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

#' Save a cells similarity matrix.
#' 
#' This function plots similarity matrix as a heatmap, so one can see similarity between parts of different clusters.
#'
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @param cellsSimilarityMatrix an output matrix from the conclus::clusterCellsInternal() function.
#' @param colorPalette "default" or a vector of colors for the column "clusters" in the colData, for example c("yellow", "#CC79A7"). 
#' @param statePalette "default" or a vector of colors for the column "state" in the colData, for example c("yellow", "#CC79A7"). 
#' @param clusteringMethod a clustering methods passed to hclust() function.
#' @param orderClusters boolean, order clusters or not.
#' @param plotPDF if TRUE export to pdf, if FALSE export to png. 
#' FALSE is recommended for datasets with more than 2500 cells due to large pdf file size.
#' @param returnPlot boolean, return plot or not. Default if FALSE.
#' @param width plot width.
#' @param height plot height.
#' @param ... other parameters of pdf(), pheatmap() and png() functions.
#'
#' @return A ggplot object or nothing (depends on the returnPlot parameter).
#' It saves the pdf in "dataDirectory/pictures" folder.
#' @export
plotCellSimilarity <- function(sceObject, cellsSimilarityMatrix, dataDirectory,
                               experimentName, colorPalette="default",
                               statePalette="default", clusteringMethod="ward.D2",
                               orderClusters = FALSE,
                               plotPDF = TRUE,
                               returnPlot = FALSE,
                               width=7, height=6, ...){

  .plotCellSimilarity(sceObject, cellsSimilarityMatrix, dataDirectory,
                      experimentName, colorPalette=colorPalette,
                               statePalette=statePalette, clusteringMethod=clusteringMethod,
                               orderClusters = orderClusters,
                               plotPDF = plotPDF,
                               returnPlot = returnPlot,
                               width=width, height=height, ...)
}

# Do not export this function.
.plotClusteredTSNE <- function(sceObject, dataDirectory, experimentName,
                              tSNEresExp = "",
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
  # plots picture based on t-SNE coordinates from
  # generateTSNECoordinates() and clustering results
  # from clusterCellsInternal() or runClustering()

  tSNEDirectory <- "tsnes"
  graphsDirectory <- "pictures"
  graphsTSNEDirectory <- "tSNE_pictures"

  if(tSNEresExp == ""){
      tSNEresExp <- experimentName
  }

  ### Plot all precalculated pSNEs to show your clusters ###

  if(columnName == "noColor"){
      numberElements <- NULL
  }else{
      numberElements <- length(unique(SummarizedExperiment::colData(sceObject)[,columnName]))
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

    coordinatesName <- paste0(tSNEresExp, '_tsne_coordinates_', i, "_",
                              PCA[i], "PCs_", perp[i], "perp")

    #fns <- list.files(file.path(dataDirectory, tSNEDirectory), full.names = TRUE)

    TSNEres <- read.delim(file.path(dataDirectory, tSNEDirectory,
                                    paste0(coordinatesName, ".tsv")),
                          stringsAsFactors = FALSE)
    #TSNEres <- read.delim(fns[grepl(paste0(PCA[i], "PCs_", perp[i], "perp.tsv"), fns)],
    #                            stringsAsFactors = FALSE)

    TSNEres <- TSNEres[rownames(TSNEres) %in% SummarizedExperiment::colData(sceObject)$cellName, ]

    if(columnName != "noColor"){
        TSNEres[columnName] <- factor(SummarizedExperiment::colData(sceObject)[,columnName])
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
        tmp <- ggplot2::ggplot(TSNEres, aes_string(x=names(TSNEres)[1],
                                          y=names(TSNEres)[2])) +
            geom_point(size=I(1)) + theme_bw()
    }else{
        tmp <- ggplot2::ggplot(TSNEres, aes_string(x=names(TSNEres)[1],
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

#' Plot t-SNE. Addtionally, it can highlight clusters or states.
#'
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @param tSNEresExp if t-SNE coordinates were generated in a different CONCLUS run, you can use them without renaming the files.
#' Please copy tsnes folder from the source run to the current one and write that experimentName in the tSNEresExp argument.
#' @param colorPalette "default" or a vector of colors for the column "clusters" in the colData, for example c("yellow", "#CC79A7").
#' @param PCs vector of PCs (will be specified in filenames).
#' @param perplexities vector of perplexities (will be specified in filenames).
#' @param columnName name of the column to plot on t-SNE dimensions.
#' @param returnPlot boolean, return plot or not.
#' @param width plot width.
#' @param height plot height.
#' @param ... other arguments of the pdf() function.
#'
#' @return A ggplot object or nothing (depends on the returnPlot parameter).
#' @export
plotClusteredTSNE <- function(sceObject, dataDirectory, experimentName,
                              tSNEresExp = "",
                              colorPalette = "default",
                              PCs=c(4, 6, 8, 10, 20, 40, 50),
                              perplexities=c(30, 40),
                              columnName="clusters",
                              returnPlot = FALSE,
                              width=6, height=5, ...){

  .plotClusteredTSNE(sceObject, dataDirectory, experimentName,
                     tSNEresExp = tSNEresExp,
                     colorPalette = colorPalette,
                     PCs=PCs,
                     perplexities=perplexities,
                     columnName=columnName,
                     returnPlot = returnPlot,
                     width=width, height=height, ...)
}

### This function returns a matrix with "protocells" representing clusters ###
### values show how much two "protocells" are similar ###
### 1 if clusters are very similar, 0 if very different ###

mkSimMed <- function(simMat, clusters){

    clusMed <- matrix(ncol=length(unique(clusters)), nrow=nrow(simMat))
    clusterNames <- levels(clusters)

    for(i in 1:ncol(clusMed)){
        clusMed[,i] <- matrixStats::rowMedians(simMat[,clusters == clusterNames[i]])
    }

    clusMed <- t(clusMed)

    simMed <- matrix(ncol=length(unique(clusters)),
                     nrow=length(unique(clusters)))

    for(i in 1:ncol(simMed)){
        simMed[,i] <- matrixStats::rowMedians(clusMed[,clusters == clusterNames[i]])
    }

    # colnames(simMed) = 1:length(unique(clusters))
    # rownames(simMed) = 1:length(unique(clusters))

    colnames(simMed) <- clusterNames
    rownames(simMed) <- clusterNames

    return(simMed)
}

#' Having cells similarity, calculate clusters similarity.
#'
#' @param cellsSimilarityMatrix a similarity matrix, one of the results of conclus::clusterCellsInternal() function.
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param clusteringMethod a clustering methods passed to hclust() function.
#'
#' @return A list contating the cluster similarity matrix and cluster names (order).
#' @export
calculateClustersSimilarity <- function(cellsSimilarityMatrix, sceObject,
                                        clusteringMethod){

    # Calculating cluster similarity for plotting picture
    # and ranking genes result is the square matrix with
    # dimension equal to number of cluster. Numbers in matrix
    # are similarity between cluster.

    clusters <- SummarizedExperiment::colData(sceObject)$clusters
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

# Don't export this function
.plotClustersSimilarity <- function(clustersSimilarityMatrix, sceObject,
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
                                       RColorBrewer::brewer.pal(n = 7,
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

    clusters <- SummarizedExperiment::colData(sceObject)$clusters
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

    annotationColors <- generateAnnotationColors(SummarizedExperiment::colData(sceObject),
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
    pheatmap::pheatmap(clustersSimilarityMatrix,
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
        pheatmap::pheatmap(clustersSimilarityMatrix,
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

#' Save a similarity cluster matrix.
#'
#' @param clustersSimilarityMatrix a matrix, result of conclus::calculateClustersSimilarity() function.
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @param colorPalette "default" or a vector of colors for the column "clusters" in the colData, for example c("yellow", "#CC79A7").
#' @param statePalette "default" or a vector of colors for the column "state" in the colData, for example c("yellow", "#CC79A7").
#' @param clusteringMethod a clustering methods passed to hclust() function.
#' @param returnPlot boolean, return plot or not.
#' @param width plot width.
#' @param height plot height.
#' @param ... other parameters of pdf() and pheatmap() functions.
#'
#' @return A ggplot object or nothing (depends on returnPlot parameter). It saves the pdf in "dataDirectory/pictures" folder.
#' @export
plotClustersSimilarity <- function(clustersSimilarityMatrix, sceObject,
                                   dataDirectory,
                                   experimentName, colorPalette,
                                   statePalette,
                                   clusteringMethod,
                                   returnPlot = FALSE,
                                   width=7, height=5.5, ...) {

  .plotClustersSimilarity(clustersSimilarityMatrix, sceObject,
                          dataDirectory, 
                          experimentName, colorPalette,
                          statePalette,
                          clusteringMethod,
                          returnPlot = returnPlot,
                          width=width, height=height, ...)
}

#' @export
deleteFile <- function(fileName, outputDir){
    file.remove(file.path(outputDir, fileName))
}

# rankGenesInternal() saves n files in your outputDir, n=number of groups
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

#' Rank marker genes by statistical significance.
#'
#' This function searches marker genes for each cluster. It saves tables in the "dataDirectory/marker_genes" directory,
#' one table per cluster.
#' 
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param clustersSimilarityMatrix matrix, result of conclus::calculateClustersSimilarity() function.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @export
#' @param column name of the column with a clustering result.
rankGenes <- function(sceObject, clustersSimilarityMatrix, dataDirectory,
                      experimentName, column="clusters"){
  # 

  markerGenesDirectory <- "marker_genes"
  rankGenesInternal(Biobase::exprs(sceObject), SummarizedExperiment::colData(sceObject), column,
             clustersSimilarityMatrix,
             file.path(dataDirectory, markerGenesDirectory), experimentName)

}

#' Get top N marker genes from each cluster. 
#' 
#' This function reads results of conclus::rankGenes() from "dataDirectory/marker_genes" and selects top N markers for each cluster.
#' 
#' @param dataDirectory output directory for a run of CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param genesNumber top N number of genes to get from one cluster.
#' @param experimentName name of the experiment which appears in filenames (supposed to be the same for one experiment during the workflow).
#' @param removeDuplicates boolean, if duplicated genes must be deleted or not.
#'
#' @return A data frame where the first columns are marker genes ("geneName") and 
#' the second column is the groups ("clusters").
#' @export
getMarkerGenes <- function(dataDirectory, sceObject, genesNumber=14,
                           experimentName, removeDuplicates = TRUE){

    markerGenesDirectory <- "marker_genes"
    numberOfClusters <- length(unique(SummarizedExperiment::colData(sceObject)$clusters))
    dir = file.path(dataDirectory, markerGenesDirectory)
    nTop = genesNumber
    clusters = unique(SummarizedExperiment::colData(sceObject)$clusters)

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
    if(removeDuplicates){
        markersClusters <- markersClusters[!duplicated(markersClusters$geneName),]
    }
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

# Do not export this function
.plotCellHeatmap <- function(markersClusters, sceObject, dataDirectory,
                            experimentName,
                            fileName, meanCentered=TRUE, colorPalette="default",
                            statePalette="default", clusteringMethod="ward.D2",
                            orderClusters = FALSE, #FALSE, TRUE, name, similarity
                            orderGenes = FALSE, # FALSE, TRUE (will be ordered the same as clusters)
                            returnPlot = FALSE,
                            saveHeatmapTable = FALSE,
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

  colData <- SummarizedExperiment::colData(sceObject)
  expressionMatrix <- Biobase::exprs(sceObject)[rownames(Biobase::exprs(sceObject)) %in%
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
  pheatmap::pheatmap(expressionMatrix,
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

  if(saveHeatmapTable){
      exportMatrix(expressionMatrix, dataDirectory, experimentName, fileName)
  }

  if(returnPlot){
      return(pheatmap::pheatmap(expressionMatrix,
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

#' Save markers heatmap.
#' 
#' This function plots heatmap with marker genes on rows and clustered cells on columns. 
#'
#' @param markersClusters a data frame where the first column is "geneName" containing genes names from sceObject, 
#' and the second column is corresponding "clusters". All names from that column must come from the column "clusters" in the colData(sceObject).
#' The data frame can be obtained from conclus::getMarkerGenes() function or created manually.
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory of a given CONCLUS run (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which appears in filenames (supposed to be the same for one experiment during the workflow).
#' @param fileName name of the ouput file
#' @param meanCentered boolean, should mean centering be applied to the expression data or not.
#' @param colorPalette "default" or a vector of colors for the column "clusters" in the colData, for example c("yellow", "#CC79A7").
#' @param statePalette "default" or a vector of colors for the column "state" in the colData, for example c("yellow", "#CC79A7").
#' @param clusteringMethod a clustering methods passed to hclust() function.
#' @param orderClusters boolean, should the heatmap be structured by clusters.
#' @param orderGenes boolean, should the heatmap be structured by genes.
#' @param returnPlot boolean, whether to return a ggplot object with the plot or not.
#' @param saveHeatmapTable boolean, whether to save the expression matrix used for heatmap into a .csv file or not.
#' The file will be saved into 'dataDirectory/output_tables' with the same name as the .pdf plot.
#' @param width plot width.
#' @param height plot height.
#' @param ... other parameters from pdf() and pheatmap() functions.
#'
#' @return A ggplot object of the plot if needed. The function saves pdf in "dataDirectiry/pictures" folder.
#' @export
plotCellHeatmap <- function(markersClusters, sceObject, dataDirectory,
                            experimentName,
                            fileName, meanCentered=TRUE, colorPalette="default",
                            statePalette="default", clusteringMethod="ward.D2",
                            orderClusters = FALSE, #FALSE, TRUE, name, similarity
                            orderGenes = FALSE, # FALSE, TRUE (will be ordered the same as clusters)
                            returnPlot = FALSE,
                            saveHeatmapTable = FALSE,
                            width=10, height=8.5, ...){

  .plotCellHeatmap(markersClusters, sceObject, dataDirectory,
                   experimentName,
                   fileName, meanCentered=meanCentered, colorPalette=colorPalette,
                   statePalette=statePalette, clusteringMethod=clusteringMethod,
                   orderClusters = orderClusters, #FALSE, TRUE, name, similarity
                   orderGenes = orderGenes, # FALSE, TRUE (will be ordered the same as clusters)
                   returnPlot = returnPlot,
                   saveHeatmapTable = saveHeatmapTable,
                   width=width, height=height, ...)
}

exportData <- function(sceObject, dataDirectory, experimentName){
  # exports all the data from workflow, including .RData files

  outputDataDirectory <- "output_tables"

  ################ EXPORT MATRIX, COLDATA, ROWDATA, FULL WORKSPACE
  write.table(Biobase::exprs(sceObject), file=file.path(dataDirectory,
            outputDataDirectory, paste0(experimentName, "_",
                                        "expression_matrix.tsv")), sep="\t",
            row.names = TRUE, quote = FALSE, col.names = TRUE)
  write.table(SummarizedExperiment::colData(sceObject), file=file.path(dataDirectory,
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

#' DBSCAN clustering on t-SNE results.
#' 
#' This function provides consensus DBSCAN clustering based on the results of t-SNE. 
#' You can tune algorithm parameters in options to get the number of clusters you want.
#'
#' @param tSNEResults the result of conclus::generateTSNECoordinates() function.
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory of a given CONCLUS run (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which appears in filenames (supposed to be the same for one experiment during the workflow).
#' @param epsilon a parameter of fpc::dbscan() function.
#' @param minPoints a parameter of fpc::dbscan() function.
#' @param k preferred number of clusters. Alternative to deepSplit.
#' @param PCs a vector of first principal components.
#' For example, to take ranges 1:5 and 1:10 write c(5, 10).
#' @param perplexities a vector of perplexity for t-SNE.
#' @param randomSeed random seed for reproducibility.
#' @param deepSplit intuitive level of clustering depth. Options are 1, 2, 3, 4.
#' @param clusteringMethod a clustering methods passed to hclust() function.
#' @param cores maximum number of jobs that CONCLUS can run in parallel.
#' @param deleteOutliers Whether cells which were often defined as outliers by dbscan must be deleted.
#' It will require recalculating of the similarity matrix of cells. Default is FALSE.
#' Usually those cells appear in an "outlier" cluster and can be easier distinguished and deleted later
#' if necessary.
#' 
#'
#' @return A list containing filtered from outliers SingleCellExperiment object and cells similarity matrix.
#' @export
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
  # It combines all the clustering parts. Takes tSNE coordinates and gives
  # results of final clustering: sceObject and cell correlation matrix

  # run dbscan
  message("Running dbscan using ", cores, " cores.")
  dbscanResults <- runDBSCAN(tSNEResults, sceObject, dataDirectory,
                             experimentName, epsilon=epsilon,
                             minPoints=minPoints, cores=cores)
  if(deleteOutliers){
      # filter cluster outliers
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

  # assign cells to cluster
  clusteringResults <- clusterCellsInternal(dbscanResultsFiltered, sceObjectFiltered,
                                    clusterNumber=k, deepSplit=deepSplit,
                                    clusteringMethod=clusteringMethod,
                                    cores=cores)
  sceObjectFiltered <- clusteringResults[[1]]
  cellsSimilarityMatrix <- clusteringResults[[2]]

  return(list(sceObjectFiltered, cellsSimilarityMatrix))
}

# getKEGGGenes
# Extracting genes from KEGG database by pathway ID.
# Returns only the genes that found in expression matrix.

# Please do not export this function because it does not work now.
.getKEGGGenes <- function(pathwayID, sceObject, species="mmu"){
  #

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

# Do not export this function.
.plotGeneExpression <- function(geneName, experimentName, dataDirectory,
                               graphsDirectory = "pictures",
                               sceObject, tSNEpicture=1,
                               commentName = "", palette = c("grey","red",
                                                             "#7a0f09",
                                                             "black"),
                               returnPlot = FALSE,
                               savePlot = TRUE,
                               alpha = 1, limits = NA,
                               pointSize = 1,
                               width=6, height=5, onefile=FALSE, #pdf
                               family, title, fonts, version,
                               paper, encoding, bg, fg, pointsize,
                               pagecentre, colormodel,
                               useDingbats, useKerning,
                               fillOddEven, compress){

  experimentName <- experimentName
  dataDirectory <- dataDirectory
  tSNEDirectory <- "tsnes"

  ### Plot all precalculated t-SNEs to show your clusters ###

  clustersNumber <- length(unique(SummarizedExperiment::colData(sceObject)$clusters))

  coordsName <- list.files(file.path(dataDirectory, tSNEDirectory),
                          pattern = paste0(experimentName,'_tsne_coordinates_',
                                           tSNEpicture, "_"))

  tSNECoords <- read.delim(file.path(dataDirectory, tSNEDirectory, coordsName),
                        stringsAsFactors=FALSE)

  tSNECoords <- tSNECoords[SummarizedExperiment::colData(sceObject)$cellName, ]

  if(!geneName %in% rownames(Biobase::exprs(sceObject))){
    print("Gene is not found in expression matrix")
  }

  stopifnot(all(rownames(tSNECoords) == colnames(sceObject)))
  tSNECoords$expression <- Biobase::exprs(sceObject)[geneName, ]

  if(length(limits) == 1){
      limits <- c(min(tSNECoords$expression),
                  max(tSNECoords$expression))
  }

  if(savePlot){
      pdf(file.path(dataDirectory, graphsDirectory, paste0(paste(experimentName,
                                                                 "tSNE", clustersNumber, "clusters" , geneName, commentName,
                                                                 "tSNEpicture", tSNEpicture, "_alpha", alpha,
                                                                 sep="_"), ".pdf")),
          width=width, height=height, onefile=onefile, # not changed by default
          family=family, title=title, fonts=fonts, version=version,
          paper=paper, encoding=encoding, bg=bg, fg=fg, pointsize=pointsize,
          pagecentre=pagecentre, colormodel=colormodel,
          useDingbats=useDingbats, useKerning=useKerning,
          fillOddEven=fillOddEven, compress=compress)
  }
  tmp <- ggplot2::ggplot(tSNECoords, aes(x=tSNECoords[,1],
                                y=tSNECoords[,2], color=expression)) +
    geom_point(size=I(pointSize), alpha = alpha) + theme_bw() +
    scale_colour_gradientn(colours=alpha(colorRampPalette(palette)(100), 0.8),
                           limits = limits) +
      ggtitle(geneName)
  #brewer.pal(9, "OrRd")[0:9]
  print(tmp)

  if(savePlot){
      dev.off()
  }

  if(returnPlot){
      return(tmp)
  }
}

#' plotGeneExpression
#'
#' The function saves a t-SNE plot colored by expression of a given gene. 
#' Warning: filename with t-SNE results is hardcoded, so please don't rename the output file.
#'
#' @param geneName name of the gene you want to plot.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which appears in filenames (supposed to be the same for one experiment during the workflow).
#' @param graphsDirectory name of the subdirectory where to put graphs. Default is "dataDirectory/pictures".
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param tSNEpicture number of the picture you want to use for plotting. 
#' Please check "dataDirectory/tsnes" or "dataDirectory/pictures/tSNE_pictures/clusters" to get the number, it is usually from 1 to 14.
#' @param commentName comment you want to specify in the filename.
#' @param palette color palette for the legend.
#' @param returnPlot boolean, should the function return a ggplot object or not.
#' @param savePlot boolean, should the function export the plot to pdf or not.
#' @param alpha opacity of the points of the plot.
#' @param limits range of the gene expression shown in the legend.
#' This option allows generating t-SNE plots with equal color
#' scale to compare the expression of different genes. By default, limits are the range
#' of expression of a selected gene.
#' @param pointSize size of the point.
#' @param width plot width.
#' @param height plot height.
#' @param ... other parameters of the pdf() function.
#'
#' @return A ggplot object of the plot if needed.
#' @export
plotGeneExpression <- function(geneName, experimentName, dataDirectory,
                               graphsDirectory = "pictures",
                               sceObject, tSNEpicture=1,
                               commentName = "", palette = c("grey","red",
                                                             "#7a0f09",
                                                             "black"),
                               returnPlot = FALSE,
                               savePlot = TRUE,
                               alpha = 1, limits = NA,
                               pointSize = 1,
                               width=6, height=5, ...){

  .plotGeneExpression(geneName, experimentName, dataDirectory,
                      graphsDirectory = graphsDirectory,
                      sceObject, tSNEpicture=tSNEpicture,
                      commentName = commentName, palette = palette,
                      returnPlot = returnPlot,
                      savePlot = savePlot,
                      alpha = alpha, limits = limits,
                      pointSize = pointSize,
                      width=width, height=height, ...)
}

#' exportClusteringResults
#'
#' The function saves clustering results into a table. Row names are cell names in the same order as in the sceObject.
#'
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which appears at the beginning of the file name 
#' (supposed to be the same for one experiment during the workflow).
#' @param fileName the rest of output file name.
#'
#' @export
exportClusteringResults <- function(sceObject, dataDirectory,
                                    experimentName, fileName){

  tableData <- S4Vectors::DataFrame(clusters = SummarizedExperiment::colData(sceObject)$clusters,
                         row.names = SummarizedExperiment::colData(sceObject)$cellName)
  write.table(tableData,
              file = file.path(dataDirectory, "output_tables",
                               paste0(experimentName,"_", fileName)),
              sep = "\t", quote = FALSE)
}

#' addClusteringManually
#'
#' The function replaces the content of the column "clusters" in the colData(sceObject) 
#' with the clustering provided in the user table.
#' The function will return the sceObject with cells which intersect with the cells from the input table.
#'
#' @param fileName a file with the clustering solution (for example, from previous CONCLUS runs).
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which appears in filenames (supposed to be the same for one experiment during the workflow).
#' @param columnName name of the column with the clusters.
#'
#' @return A SingleCellExperiment object with the created/renewed column "clusters" in the colData(sceObject).
#' @export
addClusteringManually <- function(fileName, sceObject, dataDirectory,
                                  experimentName, columnName = "clusters"){

  tableData <- read.table(file.path(dataDirectory, "output_tables",
                                paste0(experimentName,"_", fileName)), sep="\t")

  if(all(rownames(SummarizedExperiment::colData(sceObject)) %in% rownames(tableData))){
    if(ncol(tableData) == 1){
      SummarizedExperiment::colData(sceObject)$clusters <-
          factor(tableData[rownames(SummarizedExperiment::colData(sceObject)),])
    } else {
      SummarizedExperiment::colData(sceObject)$clusters <-
          factor(tableData[rownames(SummarizedExperiment::colData(sceObject)),][,columnName])
    }

    return(sceObject)
  } else {
    message("Rownames in colData are not equal to rownames in table.
Returning SCE object with cells intersecting with clusters_table.")
      sceObject <- sceObject[ ,colnames(sceObject) %in%
                                 intersect(colnames(sceObject),
                                           rownames(tableData))]
      tableData$randomColumn <- NA
      tableData <- tableData[rownames(tableData) %in%
                                 intersect(colnames(sceObject),
                                           rownames(tableData)), ]
      SummarizedExperiment::colData(sceObject)$clusters <-
          factor(tableData[rownames(SummarizedExperiment::colData(sceObject)),][,columnName])
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
                          genomeAnnot, ensemblPattern, rowData = NULL,
                          databaseDir = system.file("extdata", package = "conclus")){
    if(databaseDir == FALSE){
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
    }else{
        # only for mouse
        ensemblPattern <- "ENSMUSG"
        database <- read.delim(file.path(databaseDir,
                                         "Mmus_gene_database_secretedMol.tsv"),
                               stringsAsFactors = FALSE)
        database <- database[!duplicated(database$Symbol),]

        ensemblGenes <- rownames(countMatrix)[grep(ensemblPattern,
                                                   rownames(countMatrix))]
        ensemblGenesInternal <- gsub(paste0(".*_", ensemblPattern),
                             ensemblPattern, ensemblGenes)
        symbolGenes <- rownames(countMatrix)[!grepl(ensemblPattern,
                                                    rownames(countMatrix))]

        rowdataEnsembl <- data.frame(ensemblGenesInternal = ensemblGenesInternal,
                                     nameInCountMatrix = ensemblGenes)
        rowdataSymbol <- data.frame(nameInCountMatrix = symbolGenes)

        message(paste0("Annotating ",length(ensemblGenes), " genes containing ",
                       ensemblPattern, " pattern."))

        rowdataEnsembl <- merge(rowdataEnsembl, database,
                         by.x = "ensemblGenesInternal", by.y = "Ensembl",
                         all.x = TRUE, all.y = FALSE, sort = FALSE)
        rowdataEnsembl <- rowdataEnsembl[,-1]

        rowdataEnsembl$Ensembl <- rowdataEnsembl$nameInCountMatrix
        rowdataSymbol$Symbol <- rowdataSymbol$nameInCountMatrix

        message("Annotating rest ", length(symbolGenes), " genes
    considering them as SYMBOLs.")

        rowdataSymbol <- merge(rowdataSymbol, database,
                               by.x = "nameInCountMatrix", by.y = "Symbol",
                               all.x = TRUE, all.y = FALSE, sort = FALSE)

        rowdata <- base::rbind(rowdataSymbol, rowdataEnsembl)
        colnames(rowdata)[colnames(rowdata) == "Ensembl"] <- "ENSEMBL"
        colnames(rowdata)[colnames(rowdata) == "Symbol"] <- "SYMBOL"
        colnames(rowdata)[colnames(rowdata) == "Name"] <- "GENENAME"
        colnames(rowdata)[colnames(rowdata) == "Feature.Type"] <- "gene_biotype"
    }

    if(!is.null(rowData)){
        rowData$nameInCountMatrix <- rownames(rowData)
        rowdata <- merge(rowData, rowdata,
                         by.x = "nameInCountMatrix", by.y = "nameInCountMatrix",
                         all.x = TRUE, all.y = TRUE, sort = FALSE)
    }

    rownames(rowdata) <- rowdata$nameInCountMatrix
    rowdata <- rowdata[rownames(countMatrix),]
    rowdata$SYMBOL[(S4Vectors::isEmpty(rowdata$SYMBOL)) | (rowdata$SYMBOL == "")] <- NA
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

    # internal function, filters genes which are more than in 10 cells and less than (all-10) cells

    selRows <- ((rowSums(countMatrix[,] >= 1)) > 10)
    countMatrix <- countMatrix[selRows,]
    rowData <- rowData[rowData$nameInCountMatrix %in% rownames(countMatrix),]

    return(list(countMatrix, rowData))
}

# Do not export this function.
.normaliseCountMatrix <- function(countMatrix,
                                 species,
                                 method="default",
                                 sizes=c(20,40,60,80,100),
                                 rowData=NULL,
                                 colData=NULL,
                                 alreadyCellFiltered = FALSE,
                                 runQuickCluster = TRUE,
                                 databaseDir = system.file("extdata", package = "conclus")){
    # Does normalisation of count matrix with.
    # There are 2 possible methods: "default" or "census"
    # The function returns SCE object with normalised count matrix
    if(method == "default"){
        rowData <- annotateGenes(countMatrix, species = species,
                                 rowData = rowData, databaseDir = databaseDir)
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
            SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(countMatrix)),
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
        print(summary(SingleCellExperiment::sizeFactors(sceNorm)))
        if(length(SingleCellExperiment::sizeFactors(sceNorm)[SingleCellExperiment::sizeFactors(sceNorm) <= 0]) > 0){
            message("Cells with negative sizeFactors will be deleted before the
downstream analysis.")
        }
        sceNorm <- sceNorm[, SingleCellExperiment::sizeFactors(sceNorm) > 0]
        sceNorm <- scater::normalize(sceNorm)
        rm(sce)

        return(sceNorm)

    }else if(method == "census"){
         message("Method 'census' is currently unavailable. Please select 'default'.")
         message("Unmodified count matrix returned.")
         return(countMatrix)
    #    sceObject <- normalize_dataset(as.matrix(countMatrix))
    #    SummarizedExperiment::colData(sceObject)$cellName = rownames(SummarizedExperiment::colData(sceObject))
    #    return(sceObject)
    }else{
        message("Wrong method. Unmodified count matrix returned.")
        return(countMatrix)
    }
}

#' normaliseCountMatrix
#'
#' Create a SingleCellExperiment object and perform normalization. The same as conclus::normalizeCountMatrix.
#'
#' @param countMatrix a matrix with non-normalised gene expression.
#' @param species either 'mmu' or 'human'.
#' @param method a method of clustering: available option is "default" using scran and scater.
#' @param sizes a vector of size factors from scran::computeSumFactors() function.
#' @param rowData a data frame with information about genes
#' @param colData a data frame with information about cells
#' @param alreadyCellFiltered if TRUE, cells quality check and filtering will not be applied. 
#' However, the function may delete some cells if they have negative size factors after scran::computeSumFactors.
#' @param runQuickCluster if scran::quickCluster() function must be applied.
#' Usually, it allows to improve normalization for medium-size count matrices. 
#' However, it is not recommended for datasets with less than 200 cells and
#' may take too long for datasets with more than 10000 cells.
#' @param databaseDir a path to annotation database provided with CONCLUS called 
#' "Mmus_gene_database_secretedMol.tsv" (only for MusMusculus 'mmu').
#' The function will work also without the database but slower because it will retrieve genes info from biomaRt.
#'
#' @return A SingleCellExperiment object with normalized gene expression, colData, and rowData.
#' @export
normaliseCountMatrix <- function(countMatrix,
                                 species,
                                 method="default",
                                 sizes=c(20,40,60,80,100),
                                 rowData=NULL,
                                 colData=NULL,
                                 alreadyCellFiltered = FALSE,
                                 runQuickCluster = TRUE,
                                 databaseDir = system.file("extdata", package = "conclus") # FALSE for not using the database but download from biomaRt
                                 ){

    .normaliseCountMatrix(countMatrix,
                                 species,
                                 method=method,
                                 sizes=sizes,
                                 rowData=rowData,
                                 colData=colData,
                                 alreadyCellFiltered = alreadyCellFiltered,
                                 runQuickCluster = runQuickCluster,
                                 databaseDir = databaseDir)
}
# deleted from the description of the normaliseCountMatrix() function:
# #' @param method a method of clustering: "default" (using scran and scater) or "census" (using Census from Monocle).

#' @export
normalizeCountMatrix <- function(countMatrix,
                                 species,
                                 method="default",
                                 sizes=c(20,40,60,80,100),
                                 rowData=NULL,
                                 colData=NULL,
                                 alreadyCellFiltered = FALSE,
                                 runQuickCluster = TRUE,
                                 databaseDir = ""){

    .normaliseCountMatrix(countMatrix,
                                 species,
                                 method="default",
                                 sizes=c(20,40,60,80,100),
                                 rowData=NULL,
                                 colData=NULL,
                                 alreadyCellFiltered = FALSE,
                                 runQuickCluster = TRUE,
                                 databaseDir = "")
}

#' Collect genes information to one table.
#'
#' The function takes a data frame containing gene symbols and (or) ENSEMBL IDs and returns
#' a data frame with such information as gene name, feature type, chromosome,
#' gene IDs in different annotations, knockout information from MGI, a summary from NCBI 
#' and UniProt, and whether or not a gene belongs to GO terms containing proteins on the cell surface or 
#' involved in secretion.
#'
#' @param genes a data frame with the first column called "geneName" containing gene symbols and (or) ENSEMBL IDs.
#' Other columns are optional. For example, the second column could be "clusters" with the name of the cluster 
#' for which the gene is a marker.
#' @param databaseDir a path to the database provided with CONCLUS called "Mmus_gene_database_secretedMol.tsv".
#' @param groupBy a column in the input table used for grouping the genes in the output tables.
#' This option is useful if a table contains genes from different clusters.
#' @param orderGenes if "initial" then the order of genes will not be changed.
#' @param getUniprot boolean, whether to get information from UniProt or not. Default is TRUE.
#' Sometimes, the connection to the website is not reliable. 
#' If you tried a couple of times and it failed, select FALSE. 
#' @param silent whether to show messages from intermediate steps or not.
#' @param coresGenes maximum number of jobs that the function can run in parallel.
#'
#' @return Returns a data frame.
#' @export
getGenesInfo <- function(genes, databaseDir = system.file("extdata", package = "conclus"), 
                         groupBy = "clusters",
                         orderGenes = "initial",
                         getUniprot = TRUE,
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
    database <- database[!duplicated(database$Symbol),]

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

    myCluster <- parallel::makeCluster(coresGenes, # number of cores to use
                             type = "PSOCK") # type of cluster
    doParallel::registerDoParallel(myCluster)

    # MGI vector
    if(!silent){
        message("Collecting knockout phenotype information from MGI.")
    }

    MGI <- unname(unlist(foreach::foreach(MGIid=result$MGI.Gene.Marker.ID) %dopar% getMGIentry(MGIid)))

    #NCBI vector
    if(!silent){
        message("Retrieving info from NCBI.")
    }

    NCBI <- unname(unlist( foreach::foreach(NCBIid=result$Entrez.Gene.ID) %dopar% getNCBIentry(NCBIid) ))

    if(getUniprot){
        # Uniprot
        if(!silent){
            message("Getting summary from Uniprot.")
        }
        UniprotFunction <- unname(unlist( foreach::foreach(UniprotID=result$Uniprot.ID) %dopar% getUniprotEntry(UniprotID) ))
    }

    parallel::stopCluster(myCluster)

    #MGI <- data.frame(MGI, stringsAsFactors = FALSE)
    #NCBI <- data.frame(NCBI, stringsAsFactors = FALSE)
    #UniprotFunction <- data.frame(UniprotFunction, stringsAsFactors = FALSE)

    if(getUniprot){
        result <- cbind(result, MGI, NCBI, UniprotFunction)
    }else{
        result <- cbind(result, MGI, NCBI)
    }

    result$MGI <- as.character(result$MGI)
    result$NCBI <- as.character(result$NCBI)
    if(getUniprot){
        result$UniprotFunction <- as.character(result$UniprotFunction)
    }

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
    if(getUniprot){
        result <- result[,c("geneName", "clusters", "Name",
                        "Feature.Type",  "go_id",
                        "name_1006", "MGI", "NCBI",
                        "UniprotFunction", "chromosome_name",
                        "Symbol", "Ensembl", "MGI.Gene.Marker.ID",
                        "Entrez.Gene.ID", "Uniprot.ID")]
    }else{
        result <- result[,c("geneName", "clusters", "Name",
                        "Feature.Type",  "go_id",
                        "name_1006", "MGI", "NCBI",
                        "chromosome_name", # no "UniprotFunction"
                        "Symbol", "Ensembl", "MGI.Gene.Marker.ID",
                        "Entrez.Gene.ID", "Uniprot.ID")]
    }
    return(result)
}

#' Save gene information into a table or tables for multiple inputs.
#'
#' This function runs conclus::getGenesInfo() function for all tables into the inputDir 
#' and saves the result into the outputDir.
#'
#' @param dataDirectory a directory with CONCLUS output. You can specify either 
#' dataDirectory, then inputDir and outputDir will be hardcoded, or inputDir and outputDir only.
#' The first is recommended during running CONCLUS workflow when the second option
#' is comfortable when you created input tables with genes manually.
#' @param inputDir input directory containing text files. These files can be obtained by 
#' applying conclus::saveMarkersLists() function or created manually. Each file must be a 
#' data frame with the first column called "geneName" containing gene symbols and (or) ENSEMBL IDs.
#' @param pattern a pattern of file names to take.
#' @param outputDir output directory.
#' @param databaseDir a path to the database "Mmus_gene_database_secretedMol.tsv". It is provided with the conclus package.
#' @param sep a parameter of read.delim() function.
#' @param header whether or not your input files have a header.
#' @param startFromFile number of the input file to start with. The function approaches files one by one.
#' It uses web scraping method to collect publicly available info from MGI, NCBI and UniProt websites.
#' Sometimes, if the Internet connection is not reliable, the function can drop. 
#' In this case, it is comfortable to start from the failed file and not to redo the previous ones.
#' @param groupBy a column in the input table used for grouping the genes in the output tables.
#' @param orderGenes if "initial" then the order of genes will not be changed.
#' @param getUniprot boolean, whether to get information from UniProt or not. Default is TRUE.
#' Sometimes, the connection to the website is not reliable. 
#' If you tried a couple of times and it failed, select FALSE. 
#' @param silent whether to show messages from intermediate steps or not.
#' @param coresGenes maximum number of jobs that the function can run in parallel.
#'
#' @return It saves text files either in the 'dataDirectory/marker_genes/saveGenesInfo' or outputDir 
#' depending on whether you specify dataDirectory or (inpitDir and outputDir) explicitly.
#' @export
saveGenesInfo <- function(dataDirectory = "",
                          inputDir = "", 
                          outputDir = "", 
                          pattern = "", 
                          databaseDir = system.file("extdata", package = "conclus"),
                          sep = ";", header = TRUE,
                          startFromFile = 1, #outputFormat = c("csv", "xlsx"),
                          groupBy = "clusters", # getGenesInfo params
                          orderGenes = "initial",
                          getUniprot = TRUE,
                          silent = FALSE, coresGenes = 20){

    if(dataDirectory != ""){
        inputDir = file.path(dataDirectory, "/marker_genes/markers_lists")
        outputDir = file.path(dataDirectory, "/marker_genes/saveGenesInfo")
        pattern = "markers.csv"
        dir.create(outputDir, showWarnings=F)
    }

    filesList <- list.files(inputDir, pattern = pattern)
    filePrefix <- do.call(rbind, strsplit(filesList, "[.]"))[,1]

    for(file in startFromFile:length(filesList)) {
        cat(file, "\n")
        genes <- read.delim(file.path(inputDir, filesList[file]),
                            sep = sep, header = header,
                            stringsAsFactors = FALSE)

        result <- getGenesInfo(genes, databaseDir,
                               groupBy = groupBy,
                               orderGenes = orderGenes,
                               getUniprot = getUniprot,
                               silent = silent, coresGenes = coresGenes)

        message("Writing the output file number ", file, "\n")

        #if(any(outputFormat %in% "csv")){
            write.table(result, file = file.path(outputDir,
                                                 paste0(filePrefix[file],
                                                        "_genesInfo.csv")),
                        quote = FALSE, sep = ";", row.names = FALSE)
        #}else if(any(outputFormat %in% "xlsx")){
            #write.xlsx2(result, file = file.path(outputDir,
            #                                     paste0(filePrefix[file],
            #                                            "_genesInfo.xlsx")),
            #            row.names=FALSE)
        #}
    }
}

#' Save top N marker genes for each cluster into a format suitable for conclus::saveGenesInfo() function.
#' 
#' The function takes the output files of conclus::rankGenes(), extracts top N markers and saves
#' them into the first "geneName" column of the output table. The second column "clusters" contains the 
#' name of the corresponding cluster.
#'
#' @param experimentName name of the experiment which appears at the beginning of the file name 
#' (supposed to be the same for one experiment during the workflow).
#' @param dataDirectory experiment directory (supposed to be the same for one experiment during the workflow).
#' @param inputDir input directory, usually "marker_genes" created automatically after conclus::runCONCLUS().
#' @param outputDir output directory.
#' @param pattern a pattern of the input file names to take.
#' @param Ntop number of top markers to take from each cluster.
#'
#' @return It saves files into the outputDir. The number of files is equal to the number of clusters.
#' @export
saveMarkersLists <- function(experimentName, dataDirectory,
                             inputDir = file.path(dataDirectory, "marker_genes"),
                             outputDir = file.path(dataDirectory,
                                                   paste0("marker_genes/markers_lists")),
                             pattern = "genes.tsv", Ntop = 100){

    dir.create(outputDir, showWarnings=F)

    filesList <- list.files(outputDir, pattern = "_markers.csv")
    deletedFiles <- sapply(filesList, function(fileName) deleteFile(fileName,
                                                                    outputDir))
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
