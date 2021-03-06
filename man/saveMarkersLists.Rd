% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualisation_and_clustering_functions.R
\name{saveMarkersLists}
\alias{saveMarkersLists}
\title{Save top N marker genes for each cluster into a format suitable for conclus::saveGenesInfo() function.}
\usage{
saveMarkersLists(experimentName, dataDirectory,
  inputDir = file.path(dataDirectory, "marker_genes"),
  outputDir = file.path(dataDirectory,
  paste0("marker_genes/markers_lists")), pattern = "genes.tsv",
  Ntop = 100)
}
\arguments{
\item{experimentName}{name of the experiment which appears at the beginning of the file name 
(supposed to be the same for one experiment during the workflow).}

\item{dataDirectory}{experiment directory (supposed to be the same for one experiment during the workflow).}

\item{inputDir}{input directory, usually "marker_genes" created automatically after conclus::runCONCLUS().}

\item{outputDir}{output directory.}

\item{pattern}{a pattern of the input file names to take.}

\item{Ntop}{number of top markers to take from each cluster.}
}
\value{
It saves files into the outputDir. The number of files is equal to the number of clusters.
}
\description{
The function takes the output files of conclus::rankGenes(), extracts top N markers and saves
them into the first "geneName" column of the output table. The second column "clusters" contains the 
name of the corresponding cluster.
}
