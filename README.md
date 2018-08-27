# Supplementary-files-for-SCC
The repository includes the codes and data files to reproduce the figures and tables in the manuscript.

Some important notes are as follows:

(1) To run the codes, please first install the following packages: 
     glmnet package for matlab to perform lasso (https://web.stanford.edu/~hastie/glmnet_matlab/download.html),
     jplv7 package for matlab to perform GWR (http://www.spatial-econometrics.com),
     RNetCDF package for R to handle data files in netcdf format, 
     spgwr package for R to perform GWR,
     glmnet package for R to perform lasso,
     devtools package for R to download penreg package that solves the generalized lasso package，
     R.matlab package for R to read data stored in matlab data format (.mat).
(2)  The master matlab script run_SCC.m can be used to reproduce the figures and tables in the manuscript. Detailed 
     instructions are included in run_SCC.m.
(3)  Reproducing Figure 3 requires solving the generalized lasso, which may take days to run.
(4)  Please use R version 3.3.3. The compilation of package penreg may not be successful on other versions.  
