# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

mvclb<-function(x,y,lambda1,lambda2,lambda3,A,H,normalize=TRUE,err=0.0001)
{
  ##coordinate wise multivariate regression of  cumulated laplacian matrix
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  mx = colMeans(x)
  my = colMeans(y)
  ###求取初始beta0在lambda1下的lasso值##
  beta0=matrix(1,p,q)
  for(i in 1:q)
  {
    bb=cwlasso(x,y[,i],lambda = lambda1)
    beta0[,i]<-as.numeric(bb)
  }
  for(i in 1:p)
    for(j in 1:q)
    {if(beta0[i,j]=="NaN"){beta0[i,j]=0}}
  for(i in 1:p)
    for(j in 1:q)
    {
      if(beta0[i,j]<0.001){beta0[i,j]=0}
    }
  ###数据标准化
  if (normalize) {
    mmx=matrix(NA,n,p)
    mmy=matrix(NA,n,q)
    for(i in 1:p)
    { for(j in 1:n){mmx[j,i]=mx[i]}}
    for(i in 1:q)
    {for(j in 1:n){mmy[j,i]=my[i]}}
    x1=x-mmx
    y1 =y-mmy
  }
  x2=matrix(NA,n,p)
  for(i in 1:p)
  {x2[,i]=x1[,i]/(sqrt(sum(x1[,i]^2)))}
  for(i in 1:n)
    for(j in 1:p)
    {if(x2[i,j]=="NaN"){x2[i,j]=0}}
  y=y1;x=x2

  ######cycle######
  beta=beta0
  k=1;j=1
  a=A[j,]
  h=H[k,]
  b=t(x[,k])%*%(y[,j]-x[,-k]%*%beta[-k,j])-lambda2*(a[-j]%*%beta[k,-j])-lambda3*(h[-k]%*%beta[-k,j])
  beta[k,j]=si(b,lambda1)/(1+lambda2*A[j,j]+lambda3*H[k,k])
  t=1;k=k+1
  repeat{
    if(k>p&j==q)
    {k=1;j=1
    if(sum((beta-beta0)^2)>err)
    {
      beta0[k,j]=beta[k,j]
      a=A[j,]
      h=H[k,]
      b=t(x[,k])%*%(y[,j]-x[,-k]%*%beta[-k,j])-lambda2*(a[-j]%*%beta[k,-j])-lambda3*(h[-k]%*%beta[-k,j])
      beta[k,j]=si(b,lambda1)/(1+lambda2*A[j,j]+lambda3*H[k,k])
      k=k+1
    }
    else break
    }else if(k>p&j<q)
    {  k=1;j=j+1
    if(sum((beta-beta0)^2)>err)
    {
      beta0[k,j]=beta[k,j]
      a=A[j,]
      h=H[k,]
      b=t(x[,k])%*%(y[,j]-x[,-k]%*%beta[-k,j])-lambda2*(a[-j]%*%beta[k,-j])-lambda3*(h[-k]%*%beta[-k,j])
      beta[k,j]=si(b,lambda1)/(1+lambda2*A[j,j]+lambda3*H[k,k])
      k=k+1
    }
    else break
    }else{
      if(sum((beta-beta0)^2)>err)
      {
        beta0[k,j]=beta[k,j]
        a=A[j,]
        h=H[k,]
        b=t(x[,k])%*%(y[,j]-x[,-k]%*%beta[-k,j])-lambda2*(a[-j]%*%beta[k,-j])-lambda3*(h[-k]%*%beta[-k,j])
        beta[k,j]=si(b,lambda1)/(1+lambda2*A[j,j]+lambda3*H[k,k])
        k=k+1
      }
      else break
    }
    t=t+1
    if(t>50*p*q){break}
  }
  dd=matrix(0,p,p)
  for(i in 1:p)
  {dd[i,i]=1/sqrt(sum(x1[,i]^2)) }
  beta=dd%*%beta
  for(i in 1:p)
    for(j in 1:q)
    {if(beta[i,j]=="NaN"){beta[i,j]=0}}
  for(i in 1:p)
    for(j in 1:q)
    {
      if(abs(beta[i,j])<0.1){beta[i,j]=0}
    }
  return(list(beta=beta,t=t))
}
