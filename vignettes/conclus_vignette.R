## ------------------------------------------------------------------------
# required R >= 3.4. In the current vignette, we used R 3.5.0.
library(conclus)

## ------------------------------------------------------------------------
# setting necessary parameters
# dataDirectory is output directory
dataDirectory <- "YourOutputDirectory"
experimentName <- "Bergiers"

## ------------------------------------------------------------------------
# please do not change the path and filenames
countMatrix <- read.delim(file.path(system.file("extdata", package = "conclus"),
                                    "Bergiers_counts_matrix_filtered.tsv"), 
                          stringsAsFactors = FALSE)
colData <- read.delim(file.path(system.file("extdata", package = "conclus"), 
                                "Bergiers_colData_filtered.tsv"))

## ------------------------------------------------------------------------
# 1. Normalisation
sceObject <- conclus::normaliseCountMatrix(countMatrix, species = "mmu", 
                                  colData = colData)

## ------------------------------------------------------------------------
# checking what changed after the normalisation
dim(sceObject)

# show first columns and rows of the count matrix
SingleCellExperiment::counts(sceObject)[1:5,1:5]

# show first columns and rows of the normalized count matrix
Biobase::exprs(sceObject)[1:5,1:5]

# visualize first rows of metadata (coldata)
coldataSCE <- as.data.frame(SummarizedExperiment::colData(sceObject))
head(coldataSCE)

# visualize beginning of the rowdata containing gene information
rowdataSCE <- as.data.frame(SummarizedExperiment:::rowData(sceObject))
head(rowdataSCE)

## ----fig.height=5, fig.width=6-------------------------------------------
# 2. Test step (optional)
p <- conclus::testClustering(sceObject, dataDirectory, experimentName)

## ----fig.height=5, fig.width=6-------------------------------------------
# saved as "dataDirectory/test_clustering/test_tSNE.pdf"
p[[1]]
# saved as "dataDirectory/test_clustering/test_clustering.pdf"
p[[3]]

## ------------------------------------------------------------------------
initialisePath(dataDirectory)
# default parameters, can be selected by a user
PCs=c(4, 6, 8, 10, 20, 40, 50)
perplexities=c(30, 40)
randomSeed = 42
tSNEResults <- generateTSNECoordinates(sceObject, dataDirectory, 
                                          experimentName, PCs=PCs, 
                                          perplexities=perplexities,
                                          randomSeed = randomSeed)

## ------------------------------------------------------------------------
ncol(tSNEResults)
# the third matrix of t-SNE coordinates with PC = 8 and perplixities = 30
# it is saved as "tsnes/Bergiers_tsne_coordinates_3_8PCs_30perp.tsv"
head(tSNEResults[1,3][[1]])

## ------------------------------------------------------------------------
SummarizedExperiment::colData(sceObject)$clusters = factor(c(rep(1, 100), rep(2, 200), rep(3, (ncol(sceObject)-300) ) ))
table(SummarizedExperiment::colData(sceObject)$clusters)

## ------------------------------------------------------------------------
epsilon=c(1.3, 1.4, 1.5)
minPoints=c(3, 4)
cores=14
message("Running dbscan using ", cores, " cores.")
dbscanResults <- conclus::runDBSCAN(tSNEResults, sceObject, dataDirectory, 
                           experimentName, epsilon=epsilon, 
                           minPoints=minPoints,
                           cores=cores)

## ------------------------------------------------------------------------
dim(dbscanResults)
dbscanResults[1:7, 1:10]

## ------------------------------------------------------------------------
clusteringMethod="ward.D2"
k=10 # parameter for cutree
message("Calculating cells similarity matrix.")
cellsSimilarityMatrix <- conclus::clusterCellsInternal(dbscanResults, sceObject, clusterNumber=k, 
                                  deepSplit=deepSplit, cores=cores,
                                  clusteringMethod=clusteringMethod)[[2]]
sceObjectFiltered <- sceObject

print(table(SummarizedExperiment::colData(sceObjectFiltered)$clusters, 
              dnn=list("Cells distribuion by clusters")))

## ----fig.height=6, fig.width=7-------------------------------------------
colorPalette="default"
statePalette="default"
plotPDFcellSim = TRUE
orderClusters = FALSE
clustersNumber <- length(unique(SummarizedExperiment::colData(sceObjectFiltered)$clusters))
colorPalette <- conclus::choosePalette(colorPalette, clustersNumber)

# Plotting stability of clusters
conclus::plotCellSimilarity(sceObjectFiltered, cellsSimilarityMatrix, dataDirectory,
                 experimentName, colorPalette, 
                 orderClusters = orderClusters, 
                 statePalette = statePalette, 
                 clusteringMethod = clusteringMethod,
                 plotPDF = plotPDFcellSim,
                 returnPlot = TRUE)

## ------------------------------------------------------------------------
deepSplit = 0 # 0 to avoid cutreeDynamic, 1 to 4 to use it
deleteOutliers = FALSE
epsilon=c(1.3, 1.4, 1.5)
minPoints=c(3, 4)
cores=14
clusteringMethod="ward.D2"
k=10 # split the dendrogram with cutree function into 10 groups
clusteringResults <- conclus::runClustering(tSNEResults, sceObject, dataDirectory, 
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

print(table(SummarizedExperiment::colData(sceObjectFiltered)$clusters, 
              dnn=list("Cells distribuion by clusters")))

## ----fig.height=6, fig.width=7-------------------------------------------
colorPalette="default"
statePalette="default"
plotPDFcellSim = TRUE
orderClusters = FALSE
clustersNumber <- length(unique(SummarizedExperiment::colData(sceObjectFiltered)$clusters))
colorPalette <- conclus::choosePalette(colorPalette, clustersNumber)

# Plotting cluster stablility
conclus::plotCellSimilarity(sceObjectFiltered, cellsSimilarityMatrix, dataDirectory,
                 experimentName, colorPalette, 
                 orderClusters = orderClusters, 
                 statePalette = statePalette, 
                 clusteringMethod = clusteringMethod,
                 plotPDF = plotPDFcellSim,
                 returnPlot = TRUE)

## ------------------------------------------------------------------------
tSNEclusters <- conclus::plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
                    PCs=PCs, perplexities=perplexities, colorPalette = colorPalette,
                    columnName = "clusters", returnPlot = TRUE)
tSNEnoColor <- conclus::plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
                PCs=PCs, perplexities=perplexities, colorPalette = colorPalette,
                columnName = "noColor", returnPlot = TRUE)
if(any(colnames(SummarizedExperiment::colData(sceObjectFiltered)) %in% "state")){
  tSNEstate <- conclus::plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
                    PCs=PCs, perplexities=perplexities, colorPalette = colorPalette,
                    columnName = "state", returnPlot = TRUE)
}

## ------------------------------------------------------------------------
tSNEclusters[[5]]
tSNEnoColor[[5]]
tSNEstate[[5]]

## ----fig.height=5.3, fig.width=6.5---------------------------------------
clustersSimilarityMatrix <- 
      conclus::calculateClustersSimilarity(cellsSimilarityMatrix, 
          sceObject = sceObjectFiltered,
          clusteringMethod = "ward.D2")[[1]]
  
conclus::plotClustersSimilarity(clustersSimilarityMatrix, 
                       sceObjectFiltered,
                       dataDirectory = dataDirectory, 
                       experimentName = experimentName, 
                       colorPalette = colorPalette,
                       statePalette = statePalette,
                       clusteringMethod = clusteringMethod,
                       returnPlot = TRUE)

## ------------------------------------------------------------------------
conclus::rankGenes(sceObjectFiltered, clustersSimilarityMatrix, dataDirectory, 
            experimentName)
rankedGenesClus5 <- read.delim(file.path(dataDirectory, "marker_genes",
                               "Bergiers_cluster_5_genes.tsv"),
                               stringsAsFactors = FALSE)
head(rankedGenesClus5, n = 10)

## ------------------------------------------------------------------------
conclus::exportMatrix(cellsSimilarityMatrix, dataDirectory, experimentName, 
           "cellsSimilarityMatrix")
conclus::exportMatrix(clustersSimilarityMatrix, dataDirectory, experimentName, 
           "clustersSimilarityMatrix")

## ------------------------------------------------------------------------
# Reminder where we stopped.
# 2. Test step (optional)
#testClustering(sceObject, dataDirectory, experimentName)

# 3. Running the analysis and saving the consensus clustering solution
sceObjectCONCLUS <- conclus::runCONCLUS(sceObject, dataDirectory, experimentName, 
                               plotPDFcellSim = TRUE, # FALSE for > 2500 cells
                               k = 10,
                               cores = 14, # 14 for servers, 1 for PC
                               statePalette = c("bisque", "cadetblue2", 
                                                "coral1", "cornflowerblue"),
                               deleteOutliers = FALSE  # TRUE takes more time
                               )
conclus::exportClusteringResults(sceObjectCONCLUS, dataDirectory, experimentName, 
                        "clusters_table.tsv")

## ------------------------------------------------------------------------
# 4. Plotting heatmaps
genesNumber <- 10
markersClusters <- conclus::getMarkerGenes(dataDirectory, sceObjectCONCLUS, 
                                  experimentName = experimentName,
                                  genesNumber = genesNumber)
orderClusters <- T # F will apply hierarchical clustering to all cells
orderGenes <- T    # F will apply hierarchical clustering to all genes
meanCentered <- T  # F to show normalized counts
conclus::plotCellHeatmap(markersClusters, sceObjectCONCLUS, dataDirectory, 
                experimentName, 
                paste0("clusters",
                       length(levels(SummarizedExperiment::colData(sceObjectCONCLUS)$clusters)),
                       "_meanCentered",meanCentered,
                       "_orderClusters",orderClusters,
                       "_orderGenes",orderGenes,"_top",
                       genesNumber, "markersPerCluster"), 
                meanCentered = meanCentered, 
                colorPalette = RColorBrewer::brewer.pal(10, "Paired"),
                orderClusters = orderClusters,
                orderGenes = orderGenes,
                fontsize_row = 4,
                statePalette = c("bisque", "cadetblue2", 
                                 "coral1", "cornflowerblue"),
                color = colorRampPalette(c("#023b84","#4b97fc", 
                                           "#FEE395", 
                                           "#F4794E", "#D73027",
                                           "#a31008","#7a0f09"))(100),
                returnPlot = TRUE,
                width = 7.5, height = 6.5)

## ------------------------------------------------------------------------
orderClusters <- T # F will apply hierarchical clustering to all cells
orderGenes <- T    # F will apply hierarchical clustering to all genes
meanCentered <- F  # F to show normalized counts
conclus::plotCellHeatmap(markersClusters, sceObjectCONCLUS, dataDirectory, 
                experimentName, 
                paste0("clusters",
                       length(levels(SummarizedExperiment::colData(sceObjectCONCLUS)$clusters)),
                       "_meanCentered",meanCentered,
                       "_orderClusters",orderClusters,
                       "_orderGenes",orderGenes,"_top",
                       genesNumber, "markersPerCluster"), 
                meanCentered = meanCentered, 
                colorPalette = RColorBrewer::brewer.pal(10, "Paired"),
                orderClusters = orderClusters,
                orderGenes = orderGenes,
                fontsize_row = 4,
                statePalette = c("bisque", "cadetblue2", 
                                 "coral1", "cornflowerblue"),
                color = colorRampPalette(c("#023b84","#4b97fc", 
                                           "#FEE395", 
                                           "#F4794E", "#D73027",
                                           "#a31008","#7a0f09"))(100),
                returnPlot = TRUE)

## ------------------------------------------------------------------------
clustersTable <- read.delim(file.path(dataDirectory, "output_tables",                                                    paste0(experimentName, "_clusters_table.tsv")), 
                            stringsAsFactors = FALSE)
clustersTable$clusters[clustersTable$clusters == "9"] <- "newCluster"
clustersTable$clusters[clustersTable$clusters == "10"] <- "newCluster"
write.table(clustersTable, file.path(dataDirectory, "output_tables",                                                    paste0(experimentName, "_clusters_table.tsv")), 
            quote = FALSE, sep = "\t")

## ------------------------------------------------------------------------
# 5. Correcting clustering manually (optional)
sceObjectCONCLUS <- conclus::addClusteringManually(fileName = "clusters_table.tsv", 
    dataDirectory = dataDirectory, 
    experimentName = experimentName,
    sceObject = sceObjectCONCLUS, 
    columnName = "clusters")

# 5.1 Redo the analysis with manual clustering (optional)
sceObjectCONCLUS <- conclus::runCONCLUS(sceObjectCONCLUS, dataDirectory, experimentName, 
                        preClustered = TRUE,
                        tSNEalreadyGenerated = TRUE, # to use t-SNE coords from dataDirectory/tsnes
                        tSNEresExp = experimentName,
                        cores = 14, # 14 for servers, 1 for PC
                        statePalette = c("bisque", "cadetblue2", 
                                         "coral1", "cornflowerblue"))

## ------------------------------------------------------------------------
genesNumber <- 10
markersClusters <- getMarkerGenes(dataDirectory, sceObjectCONCLUS, 
                                  experimentName = experimentName,
                                  genesNumber = genesNumber)
orderClusters <- T # F will apply hierarchical clustering to all cells
orderGenes <- T    # F will apply hierarchical clustering to all genes
meanCentered <- T  # F to show normalized counts
conclus::plotCellHeatmap(markersClusters, sceObjectCONCLUS, dataDirectory, 
                experimentName, 
                paste0("clusters",
                       length(levels(SummarizedExperiment::colData(sceObjectCONCLUS)$clusters)),
                       "_meanCentered",meanCentered,
                       "_orderClusters",orderClusters,
                       "_orderGenes",orderGenes,"_top",
                       genesNumber, "markersPerCluster"), 
                meanCentered = meanCentered, 
                colorPalette = RColorBrewer::brewer.pal(10, "Paired"),
                orderClusters = orderClusters,
                orderGenes = orderGenes,
                fontsize_row = 4,
                statePalette = c("bisque", "cadetblue2", 
                                 "coral1", "cornflowerblue"),
                color = colorRampPalette(c("#023b84","#4b97fc", 
                                           "#FEE395", 
                                           "#F4794E", "#D73027",
                                           "#a31008", "#921912"))(100),
                returnPlot = TRUE,
                width = 7, height = 5.5)

## ------------------------------------------------------------------------
# 6. Plot gene expression in a selected tSNE plot
conclus::plotGeneExpression("Ccl3", experimentName, dataDirectory, 
                   sceObject = sceObjectCONCLUS,
                   tSNEpicture = 10, returnPlot = TRUE)
conclus::plotGeneExpression("Lama1", experimentName, dataDirectory, 
                   sceObject = sceObjectCONCLUS,
                   tSNEpicture = 10, returnPlot = TRUE)
tSNEstate[[10]]

## ------------------------------------------------------------------------
# 7. getGenesInfo example
result <- getGenesInfo(markersClusters, groupBy = "clusters",
                       getUniprot = FALSE) # please change to getUniprot = TRUE

# 7.1 save the result
outputDir <- file.path(dataDirectory, "/marker_genes/getGenesInfo")
dir.create(outputDir, showWarnings=F)
write.table(result, file = file.path(outputDir, 
                                     "Bergiers_markersClusters_top10_clusters9_genesInfo.csv"),
            quote = FALSE, sep = ";", row.names = FALSE)

## ------------------------------------------------------------------------
# 8. saveGenesInfo example
saveMarkersLists(experimentName, dataDirectory)

## ------------------------------------------------------------------------
saveGenesInfo(dataDirectory, sep = ";", header = TRUE, 
              startFromFile = 9, getUniprot = FALSE) # please change to getUniprot = TRUE

## ------------------------------------------------------------------------
saveRDS(sceObjectCONCLUS, file.path(dataDirectory, paste0("output_tables/", 
                                experimentName, "_sceObjectAfterCONCLUS.rds")))

## ------------------------------------------------------------------------
sessionInfo()

