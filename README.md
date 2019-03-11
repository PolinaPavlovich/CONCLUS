# CONCLUS: from CONsensus CLUSters to a meaningful CONCLUSion

CONCLUS is a tool for robust clustering and positive marker features selection.
The motivation was to create a user-friendly tool for single-cell RNA sequencing (scRNA-seq) data which can distinguish not only cell types 
from different tissues but also characterize heterogeneity within a tissue or capture in detail stages of a 
differentiation process. Most of the existing tools offer only one clustering solution when a slight change of
input parameters can lead to significant changes. Programs with consensus approach lack of a proper marker selection.
CONCLUS shows the dynamic how the solution change in dozens of clustering attempts which frees the user from multiple 
restarting the program and exhausting manual adjustment of settings. CONCLUS allows looking at the data at different 
levels from major groups (families of clusters) to the most detailed grouping. The consensus solution allows avoiding
picking technical artifacts of one concrete clustering or t-SNE. CONCLUS uses topological analysis for clustering and statistics for markers selection.
Union of two methods allows to reveal and avoid artefacts of each of them. A user can judge the validity of the consensus clustering 
solution with a heatmap with top positive markers. CONCLUS provides user-friendly functions which help to manually merge 
and rename clusters that make it an indispensable tool for creating concise and beautiful figures for scientific publications.
CONCLUS includes functions getGenesInfo and saveGenesInfo which facilitate manual cell state annotation. They take a list of genes
and return a table with genes and chromosome names, GO term showing whether a gene codes a cell surface or secreted protein, etc. 
GetGenesInfo performs web scraping and collects information about knockout phenotype from MGI, summaries from NCBI and UniProt.
SaveGenesInfo is a wrapper of getGeneInfo which allows saving tables for multiple input files.
These two functions can be applied not only to marker genes from scRNA-seq but also for DESeq2 results in bulk RNA-seq.


CONCLUS was developed by Polina Pavlovich and Konstantin Chukreev in the laboratory of Christophe Lancrin
in the European Molecular Biology Laboratory, EMBL Rome.

The vignette, an [example R script](https://github.com/PolinaPavlovich/CONCLUS/blob/master/Example_full_workflow.R) with all input 
files are available in the GitHub folder. We are sorry that currently we cannot render our vignette directly on GitHub because the file is too big 
but you can open the '20181221_CONCLUS_vignette.html' file with your browser after cloning or downloading the CONCLUS folder to your computer.

If you have any questions about the tool or suggestions, please contact Polina Pavlovich by *pavlovich@phystech.edu*. 

## Authors contribution:

Polina Pavlovich: development of a concept, writing the algorithm, refining the coding style, writing documentation.
Konstantin Chukreev: creating a frame of a user-friendly tool, selecting the coding style.
Christophe Lancrin: project leader, inspiration, ideas, testing the code.
Andreas Buness: critics, ideas.

Maya Shvartsman: testing the code, comments.
Kerstin Ganter: testing the code.

## Acknowledgment:

We want to thank Dr. Georg Gasteiger and Mi Lian from the Institute of System Immunobiology in Würzburg for
providing freedom, time, checking the results of CONCLUS on their unpublished data during the rotation period of P.P which allowed
 to improve the functionality of the tool. 

We thank Dr. Nina Cabezas-Wallscheid from the Max Planck Institute of Immunobiology and Epigenetics in Freiburg for provided freedom, time, 
and support in working on the tool.

We are grateful to Alexey Samosyuk from the SkolTech Institute in Moscow for useful discussions and 
kindly provided a script for an alternative method of normalization. 
