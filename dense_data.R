###This demonsrate the method for dense functional data based 
### on test statistics in eq (2) and (3) of the manuscript. Preprocessed Diffusion Tensor Imaging
### dataset is provided.

library(MASS)
library(far)
library("refund")
library("gee")
library("cvTools")
library("Matrix")

source("supplementary_functions.R")

func=as.matrix(read.csv("DTI_x.csv",header = F))
scalar=as.matrix(read.csv("DTI_y.csv",header = F))
pos=as.matrix(read.csv("DTI_pos.csv",header = F))

int=pos[2,1]-pos[1,1]
n=nrow(func)



##output the p-value for## 
FU=FU(func,scalar,n,int)
FU



