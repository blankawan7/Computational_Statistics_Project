item = 30 # items
m = 100000 # examinees
set.seed(1)
theta=rnorm(m,0,1) # ability of examinees
delta=runif(item,-1,1) # difficulty of item 
# Simulate data
x=matrix(0,m,item)
for (i in 1:item) x[,i]=1*(rlogis(m,0,1)<=(theta-delta[i]))
sum_item=colSums(x)
sum_m=rowSums(x)
# Elementary symmetric function
Gamma<-function(b)
{
  N=length(b)
  gamma=matrix(0,N+1)
  gamma[1]=1
  gamma[2]=b[1]
  for (j in 2:N){
    for (s in (j+1):2)
    { 
      gamma[s]=gamma[s]+gamma[s-1]*b[j]
    }
  }
  return(gamma)
}
# constant for parameter bi 
cst_item<-function(b,lambda)
{
  N=length(b)
  g=Gamma(b)
  num=0
  denom=g[1]*lambda[1]
  for (s in 2:(N+1))
  {
    num=num+g[s-1]*lambda[s]
    denom=denom+g[s]*lambda[s]
  }
  return(num/denom)
}
# constant for parameter lambdat
cst_score<-function(b,lambda,t)
{ 
  N=length(b)
  g=Gamma(b)
  num=g[t]
  denom=0
  
  for (s in 1:(N+1))
  {
    if (s!=t) denom=denom+g[s]*lambda[s]
  }
  return(num/denom)
}

# Gibbs sampler
set.seed(1)
# prior
b=runif(item)
lambda=runif(item+1)
# initial values
B=b
LAMBDA=lambda
#score=sapply(0:item,function(p)sum(sum_m==p))
T=1000
# iteration
for (iter in 1:T){
  # items
  for (i in 1:item)
  {
    y=rbeta(1,sum_item[i]+1,m-sum_item[i]+1)
    c1=cst_item(b[-i],lambda)
    b[i]=y/(c1*(1-y))
  }
  # scores
  for (s in 1:(item+1))
  { 
    s=3
    score=length(which(rowSums(x)==s))
    z=rbeta(1,score+1,m-score+1)
    c2=cst_score(b,lambda,s)
    lambda[s]=z/(c2*(1-z))
  }
  lambda=lambda*b[1]^(1:(item+1))
  b=b/b[1]
  lambda=lambda/lambda[1]
  B=cbind(B,b)
  LAMBDA=cbind(LAMBDA,lambda)
}


## ACF function and trace plot for lambda2 lambda4 lambda5
par(mfrow=c(2,3))
 acf(LAMBDA[2,],lag=50)
 acf(LAMBDA[4,],lag=50)
 acf(LAMBDA[5,],lag=50)
 plot(c(1:length(LAMBDA[2,])),LAMBDA[2,], type = "l")
 plot(c(1:length(LAMBDA[4,])),LAMBDA[4,], type = "l")
 plot(c(1:length(LAMBDA[5,])),LAMBDA[5,], type = "l")
 
## ACF function and trace plot for b6 b7 b8
par(mfrow=c(2,3))
  acf(B[6,],lag=50)
  acf(B[16,],lag=50)
  acf(B[26,],lag=50)
  plot(c(1:length(B[6,])),B[6,], type = "l")
  plot(c(1:length(B[16,])),B[16,], type = "l")
  plot(c(1:length(B[26,])),B[26,], type = "l")
