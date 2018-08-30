library(Matrix)
library(devtools)
devtools::install_github("jaredhuling/penreg")
library(penreg)
library(RNetCDF)
# Description: This R script is used to perform the SCC method based on the radius 
#              based nearest neighbor (denoted as SCC-RNN in the manuscript) with 
#              a radius of 0.05 and on 4-nearest neighbour (denoted as SCC-KNN)
# Please modify the paths on line 31, 41, 56, and 66 as necessary

# INPUT ARGUMENTS
# X:       array storing the p explanatory variables on n locations
# y:       vector storing the response variables on n locations
# lon:     vector storing the longitudes of n locations
# lat:     vector storing the latitudes of n locations
# D:       penalty matrix
# The input arguments for SCC_KNN and SCC_RNN are stored in the file genlasso.nc and
# genlasso_dc.nc.They will be automitically loaded by running the script.

# OUTPUT ARGUMENTS
# beta:    array storing the estimated regression coefficients at different tuning parameters
# The output arguments for SCC_KNN and SCC_RNN will be output to the file genlasso_out.nc and
# genlasso_dc_out.nc.


# Perform 100 simulations
for (i in 1:100) {
  print(i)
  
# Estimate the regression coefficient using SCC-RNN
  filename=paste0("./output_R/genlasso_dc",as.character(i),".nc")
  nc<-open.nc(filename)
  X<-var.get.nc(nc,"X")
  y<-var.get.nc(nc,"y")
  D<-var.get.nc(nc,"L")
  
  
  close.nc(nc)
  res <- admm.genlasso(X, y, D = D,standardize = FALSE, intercept = FALSE,nlambda=200,lambda.min.ratio=1e-8)
  
  filesave=paste0("./output_R/genlasso_dc_out",as.character(i),".nc")
  nc <- create.nc(filesave)
  varnum<-dim(D)
  varnum<-varnum[2]
  dim.def.nc(nc, "vdim", varnum)
  dim.def.nc(nc, "tdim", 200)
  ##  Create two variables, one as coordinate variable
  var.def.nc(nc, "lambda", "NC_DOUBLE", 1)
  var.def.nc(nc, "beta", "NC_DOUBLE", c(0,1))
  var.put.nc(nc,"lambda", res$lambda)
  var.put.nc(nc,"beta", res$beta[2:(varnum+1),1:200])
  close.nc(nc)
#-------------------------------------------------------

  # Estimate the regression coefficient using SCC-KNN
  filename=paste0("./output_R/genlasso",as.character(i),".nc")
  nc<-open.nc(filename)
  X<-var.get.nc(nc,"X")
  y<-var.get.nc(nc,"y")
  D<-var.get.nc(nc,"L")
  
  
  close.nc(nc)
  res <- admm.genlasso(X, y, D = D,standardize = FALSE, intercept = FALSE,nlambda=200,lambda.min.ratio=1e-8)
  
  filesave=paste0("./output_R/genlasso_out",as.character(i),".nc")
  nc <- create.nc(filesave)
  varnum<-dim(D)
  varnum<-varnum[2]
  dim.def.nc(nc, "vdim", varnum)
  dim.def.nc(nc, "tdim", 200)
  ##  Create two variables, one as coordinate variable
  var.def.nc(nc, "lambda", "NC_DOUBLE", 1)
  var.def.nc(nc, "beta", "NC_DOUBLE", c(0,1))
  var.put.nc(nc,"lambda", res$lambda)
  var.put.nc(nc,"beta", res$beta[2:(varnum+1),1:200])
  close.nc(nc)
  #-------------------------------------------------------  
  
}

