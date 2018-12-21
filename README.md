# CONCLUS

CONCLUS (CONsensus CLUStering) is a tool for robust clustering and positive marker features selection.
The motivation was to create a user-friendly tool for scRNA-seq data which can distinguish not only cell types 
from different tissues but also characterize heterogeneity within a tissue or capture in detail stages of a 
differentiation process. Most of the existing tools offer only one clustering solution when a slight change of
input parameters can lead to significant changes. Tools with consensus approach lack of a proper marker selection.
CONCLUS shows the dynamic how the solution change in dozens of clustering attempts which frees the user from multiple 
restarting the program and exhausting manual adjustment of settings. CONCLUS allows to look at the data at different 
levels from major groups (families of clusters) to the most detailed grouping. The consensus solution allows to avoid
picking technical artefacts of one concrete clustering or t-SNE. CONCLUS uses topological analysis for clustering and statistics for markers selection.
Union of two methods allows to reveal and avoid artefacts of each of them. A user can judge the validity of the consensus clustering 
solution with a heatmap with top positive markers. CONCLUS provides a user-friendly functions which help to manually merge 
and rename clusters that makes it an indispensable tool for creating concise and beautiful figures for scientific publications.


CONCLUS was developed by Polina Pavlovich and Konstantin Chuckreev in the laboratory of Christophe Lancrin
in the European Molecular Biology Laboratory, EMBL Rome.

The vignette (link), an example R script (link) with all input and output files are available in the GitHub folder.

If you have any questions, please contact Polina Pavlovich by pavlovich@phystech.edu. 

Authors contribution:

Polina Pavlovich: development of a concept, writing the algorithm, refining the coding style, writing documentation
Konstantin Chuckreev: creating a frame of a user-friendly tool, selecting the coding style
Christophe Lancrin: project leader, inspiration, ideas, testing the code
Andreas Buness: critics, ideas

Maya Shvartsman: testing the code, comments
Kerstin Ganter: testing the code

Acknowledgment:

We want to thank Dr. Georg Gasteiger and Mi Lian from the Institute of System Immunobiology in Würzburg for
providing freedom, time, checking the results of CONCLUS on their unpublished data during the rotation period of P.P which allowed
 to improve the functionality of the tool. 

We thank Dr. Nina Cabezas-Wallscheid from the Max Planck Institute of Immunogiology and Epigenetics in Freiburg for provided freedom, time, 
and support in working on the tool.

We are grateful to Alexey Samosyuk from the SkolTech Institute from Moscow for useful discussions and 
kindly provided a script for an alternative method of normalization. 
