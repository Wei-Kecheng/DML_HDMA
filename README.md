# Debiased machine learning for ultra-high dimensional mediation analysis
------------------------------------------------
# Description
"ML_SS_ranger" is an R function that utilizes random forests to estimate nonparametric functions. Users can replace it with other machine learning approaches as needed. "DML_HDMA" is an R function designed for ultra-high dimensional mediation analysis, incorporating screening, regularization, estimation, and inference.

# Usage
- ML_SS(Z,Y,z,tree)
- DML_HDMA(x,m,y,z,K)
  
# Arguments
ML_SS(Z,Y,z,tree)
- Z: Predictors in training data, with matrix $n_1\times q$ (sample size $n_1$ and dimension $q$)
- Y: Outcome in training data, with vector $n_1$
- Z: Predictors in testing data, with matrix $n_2\times q$ (sample size $n_2$)
- Tree: Number of trees
  
DML_HDMA(x,m,y,z,K)
- x: Exposure, with vector $n$ (sample size $n$)
- m: Mediators, with matrix $n\times p$ (dimension $p$)
- y: Outcome, with vector $n$
- z: Confounders, with matrix $n\times q$ (dimension $q$)
- k: Number of folds

# Value
ML_SS(Z,Y,z,tree)
- Output predicted outcome in testing data, with vector $n_2$

DML_HDMA(x,m,y,z,K)
- Output a matrix where each row represents the selected mediators after screening and regularization, and the columns contain the estimation and inference results for the direct and indirect effects

# Examples
```r
library(ranger) #random forests
library(glmnet) #lasso
n=200 #sample size
p=100 #dimension of mediators
q=10 #dimension of confounders
#nonparametric confounding effects###########################################
z=matrix(runif(n*q,-1,1),ncol=q)
h=z[,1]+z[,2]^2+4*(z[,3]>0)
g=h
f=h
#true coeffcients############################################################
alpha=1
gamma=c(1,1,1,1,1,1,0,0,rep(0,p-8))*(-1)
beta=c(1,1,1,1,0,0,1,1,rep(0,p-8))*(-1)
#exposure+mediators+outcome###################################################
x=h+rnorm(n)
m=as.matrix(x)%*%gamma+matrix(rep(g,p),ncol=p)+matrix(rnorm(n*p),ncol=p)
y=x*alpha+m%*%beta+f+rnorm(n)
##############################################################################
source("ML_SS_ranger.R")
source("DML_HDMA.R")
result=DML_HDMA(x,m,y,z,K=5)
