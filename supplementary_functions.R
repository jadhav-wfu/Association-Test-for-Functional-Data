
count_st=function(colno,y,x,n)
{
  yt=y
  xt=x[,colno]
  tmp1=matrix(0,ncol=n,nrow=n)
  for( i in 1: (n-1) )
  {
    #print(i)
    tmp2=(i+1):n
    tmp1[i,tmp2]=sapply(tmp2,function(v) as.numeric((xt[i]-xt[v])*(yt[i]-yt[v])>=0))
  }
  
  return(sum(tmp1))
  
}



CDF=function(a,ai)
{
  tmp=which(a<=ai)
  return(length(tmp)/length(a))
}


FU=function(x,y,n,int)
{
  stat=sapply(1:ncol(x),function(v) count_st(v,y,x,n))
  stat=stat/dim(combn(n,2))[2]
  
  
  stat=(stat-0.5)*sqrt(n)/2
  stat_norm=crossprod(stat,stat)*int
  
  rvy=sapply(y,function(v) CDF(y,v))
  z=matrix(0,ncol=ncol(x),nrow=nrow(x))
  for( t in 1: ncol(x))
  {
    tmp1=sapply(x[,t],function(v) CDF(x[,t],v))
    z[,t]=sapply(1:n, function(v) tmp1[v]*rvy[v]+ (1-tmp1[v])*(1-rvy[v]) )
    
  }
  
  z=z-0.5
  covz=cov(z)
  svd=eigen(covz)
  base=svd$vectors   ### cols are eigenvectors###
  eigval=svd$values/ncol(z)
  base=sqrt(ncol(z))*orthonormalization(base, basis=T, norm=TRUE)
  tot=sum(eigval)
  cumsum=sapply(1:length(eigval), function(v) sum(eigval[1:v])/tot)
  pn=which(cumsum>0.99)[1]
  
  tmp=rchisq(pn*10000,1,0)
  nor=matrix(tmp,nrow=10000,ncol=pn)
  dist=apply(nor,1,function(v)crossprod(eigval[1:pn],v))
  fu_p=length(which(dist>as.numeric(stat_norm)))/10000
  
  return(fu_p)
  
}