###This demonsrate the method for sparse functional data based 
### on test statistics in eq (4) of the manuscript. It is applied on the
### ebay dataset.





library(MASS)
library(fda)
library(svd)
library(gskat)
library(far)
library(cvTools)
library(fdapace)
source("supplementary_functions.R")


#### preparing the data for the FPCA function available in fdapace package###

dat=read.csv("live_bids_open.csv")
ID=unique(dat$ID)
n=length(unique(ID))
Lt=Lx=list()
max(dat$Time)
Y=NULL



for(i in 1: n)
{
  tmp=which(dat$ID==ID[i])
  Lt[[i]]=dat$Time[tmp]
  Lx[[i]]=(dat$Bids[tmp])
  Y[i]=(dat$Opening[tmp[1]])
}


#### using FPCA to obtain funcitions from the predicted scores  (see page 6) of the manuscript)

res= FPCA(Lx, Lt,list(dataType='Sparse',  verbose=TRUE))
pos=res$workGrid
phi=res$phi
mu=res$mu
fpc=res$xiEst

xtilde=matrix(0,nrow=n,ncol=length(pos))
for( i in 1:n)
{
  
  xtilde[i,]=mu+apply(phi,1,function(v) crossprod(v,fpc[i,]))
  
}



#### implementing the proposed test using test statistic (4) #####

x=xtilde
y=Y
int=pos[2]-pos[1]
fu_p=FU(x,y,n,int)


###This is the p-value###
fu_p





