# Phenotypic and genetic correlation between sleep, behavior, and macroscale cortical grey matter - analysis

This repository is a companion to our article "Phenotypic and genetic correlation between sleep, behavior, and macroscale cortical grey matter". Here, we provide the basic scripts of the analysis as well as outcomes of the genetic correlation analysis. 
Phenotypic analysis can only be reran after downloading the respective data from the Human Connectome Project (https://www.humanconnectome.org/study/hcp-young-adult) and enhanced NKI dataset (http://fcon_1000.projects.nitrc.org/indi/enhanced/studies.html).

# How to run

The analysis are performed with help of: colorbrewer (https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab) - for colormaps, Spearman correlations and confidence intervals were computed using 'Spearman.m', which is downloaded from https://github.com/CPernet/robustcorrtool, SurfStat:(http://www.math.mcgill.ca/keith/surfstat/), solar-eclipse: https://www.nitrc.org/projects/se_linux/, as well as pedigree information, PLS was run based on the code of https://miplab.epfl.ch/index.php/software/PLS (included), and nifti's for the brain map analysis printed using https://de.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image, brainmap toolbox can be found here: http://www.brainmap.org. 


