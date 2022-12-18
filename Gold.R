rm(list=ls(all=TRUE)) 
#setwd("C:/Users/lcgd4/iCloudDrive/Documents/Mizzou/Dissertation 1 - plm/Programs/cpp")
setwd('/Users/lichen/Documents/Mizzou/Dissertation 1 - plm/Programs/cpp')

library(MASS)
library(Matrix)
library(stargazer)
library(Rcpp)

sourceCpp("cpart_2.cpp")

Nsample=1000
n=900
g=-2
b=2
a=1

t=list()
x=list()
y=list()
z=list()
b.hat=matrix(0, nrow=Nsample, ncol=3)
se.b=b.hat
sample=1

K1=function(x,h){
  y=0.75*(1-(x/h)^2)/h*(abs(x)<h)
  #y=dnorm(x/h)/h
  return(y)
}

S.matrix=function(ty,h){
  ty=as.matrix(unlist(ty))
  n=length(ty)
  K=matrix(NA,ncol=1,nrow=n)
  S=matrix(NA, nrow=n,ncol=n)
  for(i in 1:n){
    K=K1(ty[i]-ty,h)
    #w1=as.numeric(t(ty-ty[i])%*%K); w2=as.numeric(t((ty-ty[i])^2)%*%K); w=as.numeric(t(ty-ty[i]-w2/w1)%*%K); S[i,]=t((ty-ty[i]-w2/w1)*K)/w
    s1=as.numeric(t(ty[i]-ty)%*%K); s2=as.numeric(t((ty[i]-ty)^2)%*%K); w=as.numeric(K*(s2-(ty[i]-ty)*s1)); S[i,]=t(w)/sum(w)
  }
  return(S)
}

start=proc.time()

repeat{
  nobs=rpois(n,5)+1
  for(i in 1:n){
    t[[i]]=as.matrix(runif(nobs[i]))   # x obs time
    corr=exp(-abs(rep(1,nobs[i])%*%t(t[[i]])-t[[i]]%*%t(rep(1,nobs[i]))))
    corr.e=2^(-abs(rep(1,nobs[i])%*%t(t[[i]])-t[[i]]%*%t(rep(1,nobs[i]))))
    
    MX=sqrt(t[[i]])
    #MX=rep(1,nobs[i])
    
    #MZ=.5+t[[i]]^.5
    #MZ=rep(2, nobs[i])
    MZ=2*sin(2*pi*t[[i]])
    
    x[[i]]=mvrnorm(1, MX, corr)
    z[[i]]=mvrnorm(1, MZ, corr)
    eps=as.matrix(mvrnorm(1,rep(0,nobs[i]),corr.e))
    y[[i]]=a+g*z[[i]]+x[[i]]*b+eps
  }
  
  x.uls=as.matrix(unlist(x))
  z.uls=as.matrix(unlist(z))
  X=cbind(rep(1, length(x.uls)), z.uls, x.uls)
  Y=unlist(y)
  
  b.hat[sample,]=solve(t(X)%*%X)%*%(t(X)%*%Y)
  
  U.sq.g=0;U.dot.g=0
  U.sq.b=0;U.g=0;U.b=0
  U1=0;U2=0
  for (i in 1:n) {
    x.temp=cbind(1, z[[i]], x[[i]])
    temp=t(x.temp)%*%(y[[i]]-x.temp%*%b.hat[sample,])
    U.sq.g=U.sq.g+temp%*%t(temp)
    U.dot.g=U.dot.g+t(x.temp)%*%x.temp
  }
  #b.hat[sample]=U1/U2
  se.b[sample,]=diag(solve(U.dot.g)%*%U.sq.g%*%solve(U.dot.g))^.5
  
  if(sample %% 50 == 0) print(sample)
  sample=sample+1
  if(sample>Nsample) break
}

proc.time()-start

B=matrix(c(a,g,b), ncol= 3, nrow=Nsample, byrow = T)
bias=apply(b.hat-B, 2, mean)
SD=apply(b.hat-B, 2, sd)
SE=apply(se.b, 2, mean)
CP=round(apply(B>b.hat-1.96*se.b & B<b.hat+1.96*se.b, 2, mean)*100)


stargazer(c(n,t(matrix(t(matrix(c(bias, SD, SE, CP), nrow=3)), byrow=T))))

