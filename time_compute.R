library(Matrix)
library(devtools)
#devtools::install_github("jaredhuling/penreg")
library(penreg)
library(RNetCDF)
library(spgwr) 
library(glmnet)
library(R.matlab)



# Description
# This R script is used to evalute the computational time for SCC-MST, SCC-KNN, SCC-RNN and GWR. Note that the 
# computational time is machine dependent.  
# Please modify the paths on line 43, 66, and 86 as necessary

data=readMat("./data/cluster_w.mat");
X=data$x;y=data$y;lon=data$lon;lat=data$lat

s=10 # number of simulation

# Compute the time for GWR
dt_GWR=vector(,s)
for (i in 1:s) {
  
  print(i)
  Xi=X[i,,]
  yi=y[i,]
  coord=cbind(lon,lat)
  ptm=proc.time()
  bwG <- gwr.sel(yi ~ -1+Xi,  coords=coord, gweight = gwr.Gauss, verbose = FALSE)
  ccc=proc.time()-ptm
  dt_GWR[i]=ccc[1]
}
#----------------------------

# Compute the time for MST
dt_MST=vector(,s)
for (i in 1:s) {
  
  print(i)
  
  filename=paste0("./output_R/lasso",as.character(i),".nc")
  nc<-open.nc(filename)
  G<-var.get.nc(nc,"G")
  y<-var.get.nc(nc,"y")
  close.nc(nc)
  varnum=dim(G)
  p=varnum[2]
  q=p/varnum[1]
  ptm=proc.time()
  glmnet(G, y, family=c("gaussian"), alpha = 1, nlambda = 200, lambda.min.ratio = 1e-8, 
         standardize = FALSE, intercept=FALSE, penalty.factor = c(rep(1, p-q),rep(0,q)))
  ccc=proc.time()-ptm
  dt_MST[i]=ccc[1]
}
#-----------------------------


# Compute the time for KNN
dt_KNN=vector(,s)
for (i in 1:s) {
  
  print(i)
  
  filename=paste0("./output_R/genlasso",as.character(i),".nc")
  nc<-open.nc(filename)
  X<-var.get.nc(nc,"X")
  y<-var.get.nc(nc,"y")
  D<-var.get.nc(nc,"L")
  close.nc(nc)
  
  ptm=proc.time()
  res <- admm.genlasso(X, y, D = D,standardize = FALSE, intercept = FALSE,nlambda=200,lambda.min.ratio=1e-8)
  ccc=proc.time()-ptm
  dt_KNN[i]=ccc[1]
}
#----------------------------

# Compute the time for RNN
dt_RNN=vector(,s)
for (i in 1:s) {
  
  print(i)
  
  filename=paste0("./output_R/genlasso_dc",as.character(i),".nc")
  nc<-open.nc(filename)
  X<-var.get.nc(nc,"X")
  y<-var.get.nc(nc,"y")
  D<-var.get.nc(nc,"L")
  close.nc(nc)
  
  ptm=proc.time()
  res <- admm.genlasso(X, y, D = D,standardize = FALSE, intercept = FALSE,nlambda=200,lambda.min.ratio=1e-8)
  ccc=proc.time()-ptm
  dt_RNN[i]=ccc[1]
}
#---------------------------
