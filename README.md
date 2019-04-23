# ScRNA-seq workflow CONCLUS: from CONsensus CLUSters to a meaningful CONCLUSion

## How to install

#Required R version: 3.4.0 <= R < 3.7.0  
#install.packages("devtools")  
library(devtools)  
install_github("PolinaPavlovich/CONCLUS")  
library(conclus)  

Note: if you install conclus on a local Linux machine you might need to install rJava. A simple instruction how to do it you can find here https://www.r-bloggers.com/installing-rjava-on-ubuntu/. On Windows PC, during the installation of conclus, in case if needed, R will return you a link to the website of Java Oracle to download the update.

## Introduction
CONCLUS is a tool for robust clustering and positive marker features selection.
The motivation was to create a user-friendly tool for single-cell RNA sequencing (scRNA-seq) data which can distinguish not only cell types from different tissues but also characterize heterogeneity within a tissue or capture in detail stages of a 
differentiation process. Most of the existing tools offer only one clustering solution when a slight change of
input parameters can lead to significant changes. CONCLUS shows how the solution change in dozens of clustering attempts which frees the user from multiple restarting the program and exhausting manual adjustment of settings. CONCLUS allows looking at the data at different 
levels from major groups (families of clusters) to the most detailed grouping. The consensus solution allows avoiding
picking technical artifacts of one concrete clustering or t-SNE. CONCLUS uses topological analysis for clustering and statistics for markers selection. Union of two methods allows to reveal and avoid artifacts of each of them. A user can judge the validity of the consensus clustering solution with a heatmap with top positive markers. CONCLUS provides user-friendly functions which help to manually merge and rename clusters that make it an indispensable tool for creating concise and beautiful figures for scientific publications.

CONCLUS includes functions getGenesInfo and saveGenesInfo which facilitate manual cell types annotation. They take a list of genes
and return a table with genes and chromosome names, GO term showing whether a gene codes a cell surface or secreted protein, etc. 
GetGenesInfo performs web scraping and collects information about knockout phenotype from MGI, summaries from NCBI and UniProt websites.
SaveGenesInfo is a wrapper of getGeneInfo which allows saving tables for multiple input files. These two functions can be applied not only to marker genes from scRNA-seq but also for DESeq2 results in bulk RNA-seq.


CONCLUS was developed by Polina Pavlovich and Konstantin Chukreev in the laboratory of Christophe Lancrin
in the European Molecular Biology Laboratory, EMBL Rome.

The pdf reference manual for CONCLUS is 'conclus_manual.pdf'.  
The vignette 'conclus_vignette.nb.html' is located in the 'vignettes' folder. We are sorry that currently we cannot render our vignette directly on GitHub because the file is too big but you can open the html file with your browser after cloning or downloading the CONCLUS folder to your computer.

If you do not want to download the folder with the entire package but quickly get only the documentation, you can check the project https://github.com/PolinaPavlovich/CONCLUS_documentation with the files 'conclus_manual.pdf' and 'conclus_vignette.nb.html.'

If you have any questions about the tool or suggestions, please contact Polina Pavlovich by *pavlovich@phystech.edu*. 

## Authors contribution:

Polina Pavlovich: development of a concept, writing the algorithm, refining the coding style, writing documentation, building the package.  
Konstantin Chukreev: creating a frame of a user-friendly tool, selecting the coding style, writing documentation.  
Dr. Christophe Lancrin: project leader, inspiration, ideas, testing the code, writing documentation.

Andreas Buness: critics, ideas.  
Dr. Maya Shvartsman: testing the code, comments.  
Kerstin Ganter: testing the code.  
Dr. Nicolas Descostes: comments.

## Acknowledgment:

We want to thank Dr. Georg Gasteiger and Mi Lian from the Institute of System Immunobiology in Würzburg for
providing freedom, time, checking the results of CONCLUS on their unpublished data during the rotation period of P.P. which allowed to improve the functionality of the tool. 

We thank Dr. Nina Cabezas-Wallscheid from the Max Planck Institute of Immunobiology and Epigenetics in Freiburg for provided freedom, time, and support in working on the tool. 

We are grateful to Dr. Emmy Tsang from EMBL Rome for useful discussion and idea for the algorithm.
