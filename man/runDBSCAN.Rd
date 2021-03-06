% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualisation_and_clustering_functions.R
\name{runDBSCAN}
\alias{runDBSCAN}
\title{Run clustering iterations with selected parameters using DBSCAN.}
\usage{
runDBSCAN(tSNEResults, sceObject, dataDirectory, experimentName,
  cores = 14, epsilon = c(1.3, 1.4, 1.5), minPoints = c(3, 4))
}
\arguments{
\item{tSNEResults}{results of conclus::generateTSNECoordinates() function.}

\item{sceObject}{a SingleCellExperiment object with your experiment.}

\item{dataDirectory}{output directory for CONCLUS (supposed to be the same for one experiment during the workflow).}

\item{experimentName}{name of the experiment which will appear in filenames 
(supposed to be the same for one experiment during the workflow).}

\item{cores}{maximum number of jobs that CONCLUS can run in parallel.}

\item{epsilon}{a fpc::dbscan() parameter.}

\item{minPoints}{a fpc::dbscan() parameter.}
}
\value{
A matrix of DBSCAN results.
}
\description{
This function returns a matrix of clustering iterations of DBSCAN.
}
