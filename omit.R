rm(list=ls(all=TRUE)) 
#setwd("C:/Users/lcgd4/iCloudDrive/Documents/Mizzou/Dissertation 1 - plm/Programs/cpp")
setwd('/Users/lichen/Documents/Mizzou/Dissertation 1 - plm/Programs/cpp')

library(MASS)
library(Matrix)
library(stargazer)
library(Rcpp)

sourceCpp("cpart_2.cpp")

Nsample=100
n=200
h=n^-.6
g=-1
b=2
a=1
xz.corr=T
xz.indep=F

t=list()
x=list()
y=list()
z=list()
b.hat=matrix(0, nrow=Nsample, ncol=2)
b.hat.plm=rep(0, Nsample)
b.hat.2sg=b.hat.plm
a.hat.plm=list()
se.b=b.hat
se.b.plm=b.hat.plm
se.b.2sg=se.b.plm

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
    
    #MZ=.5+t[[i]]
    #MZ=rep(2, nobs[i])
    MZ=2*sin(2*pi*t[[i]])
    #MZ=(2*t[[i]]^2+1)^-.5
    
    if(xz.corr==F){
      if(xz.indep==T){
        x[[i]]=mvrnorm(1, MX, corr)
        z[[i]]=mvrnorm(1, MZ, corr)
        eps=as.matrix(mvrnorm(1,rep(0,nobs[i]),corr.e))
      } else {
        eta=mvrnorm(1, rep(0, nobs[i]), corr)
        z[[i]]=MZ+eta
        x[[i]]=MX+rnorm(1)*eta
        eps=rnorm(1)*eta
      }
    } else {
      Corr=matrix(0, nrow=2*nobs[i], ncol=2*nobs[i])
      Corr[1:nobs[i], 1:nobs[i]]=corr
      Corr[-(1:nobs[i]), -(1:nobs[i])]=corr
      Corr[1:nobs[i], -(1:nobs[i])]=corr/2
      Corr[-(1:nobs[i]), 1:nobs[i]]=corr/2
      xz.temp=mvrnorm(1, c(MX, MZ), Corr)
      x[[i]]=xz.temp[1:nobs[i]]
      z[[i]]=xz.temp[-(1:nobs[i])]
      eps=as.matrix(mvrnorm(1,rep(0,nobs[i]),corr.e))
    }
    y[[i]]=a+g*z[[i]]+x[[i]]*b+eps
  }
  
  ####################### naive   ####################### 
  
  x.uls=as.matrix(unlist(x))
  X=cbind(rep(1, length(x.uls)), x.uls)
  Y=unlist(y)
  
  b.hat[sample,]=solve(t(X)%*%X)%*%(t(X)%*%Y)
  
  U.sq.g=0;U.dot.g=0
  U.sq.b=0;U.g=0;U.b=0
  U1=0;U2=0
  for (i in 1:n) {
    x.temp=cbind(1, x[[i]])
    temp=t(x.temp)%*%(y[[i]]-x.temp%*%b.hat[sample,])
    U.sq.g=U.sq.g+temp%*%t(temp)
    U.dot.g=U.dot.g+t(x.temp)%*%x.temp
    if(0){
      for (j in 1:ny[i]) {
        for (k in 1:nx[i]) {
          K=K1(tx[[i]][k]-ty[[i]][j], h)
          U1=U1+K*x[[i]][k]*(y[[i]][j]-g.hat[sample]*z[[i]][j])
          U2=U2+K*x[[i]][k]^2
        }
      }
    }
  }
  #b.hat[sample]=U1/U2
  se.b[sample,]=diag(solve(U.dot.g)%*%U.sq.g%*%solve(U.dot.g))^.5
  
  
  ####################### PLM  ####################### 
  
  x.uls=as.matrix(unlist(x))
  y.uls=as.matrix(unlist(y))
  t.uls=as.matrix(unlist(t))
  N=length(x.uls)
  
  S=S.matrix(t.uls,h)
  beta=solve(t(x.uls)%*%t(diag(1,N)-S)%*%(diag(1,N)-S)%*%x.uls)%*%
    t(x.uls)%*%t(diag(1,N)-S)%*%(diag(1,N)-S)%*%y.uls
  b.hat.plm[sample]=beta
  a.hat.plm[[sample]]=S%*%(y.uls-x.uls%*%beta)
  #sum(is.na(S))
  
  if(1){
    e=list(NA)
    j=0
    for(i in 1:n){
      x[[i]]=na.omit(x[[i]])
      y[[i]]=na.omit(y[[i]])
      t[[i]]=na.omit(t[[i]])
      ind=c(j+1,j+length(x[[i]]))
      e[[i]]=y[[i]]-(a.hat.plm[[sample]][ind[1]:ind[2]]+x[[i]]%*%beta)
      e[[i]]=e[[i]]%*%t(e[[i]])
      j=j+length(x[[i]])
    }
    
    D=t(x.uls)%*%t(diag(1,N)-S)%*%(diag(1,N)-S)%*%x.uls
    C=bdiag(e)
    V=t(x.uls)%*%t(diag(1,N)-S)%*%C%*%(diag(1,N)-S)%*%x.uls
    se.b.plm[sample]=sqrt(as.numeric(solve(D)%*%V%*%solve(D)))
  }
  
  
  ####################### two stage  ####################### 
  
  if(1){
    x.center=x
    y.center=y
    if(1){
      x.uls=unlist(x)
      t.uls=unlist(t)
      y.uls=unlist(y)
      Ex=list()
      Ey=list()
      for(i in 1:n){
        Ex[[i]]=x[[i]]
        Ey[[i]]=y[[i]]
        for (k in 1:nobs[i]) {
          temp=K1(t[[i]][k]-t.uls, h)/sum(K1(t[[i]][k]-t.uls, h))
          Ex[[i]][k]=t(x.uls)%*%temp
          Ey[[i]][k]=t(y.uls)%*%temp
        }
        x.center[[i]]=x[[i]]-Ex[[i]]
        y.center[[i]]=y[[i]]-Ey[[i]]
      }
    }
    
    X=unlist(x.center)
    Y=unlist(y.center)
    
    b.hat.2sg[sample]=t(X)%*%Y/(t(X)%*%X)
    
    U.sq.g=0;U.dot.g=0
    U.sq.b=0;U.g=0;U.b=0
    U1=0;U2=0
    for (i in 1:n) {
      U.sq.g=U.sq.g+(sum(x.center[[i]]*(y.center[[i]]-x.center[[i]]*b.hat.2sg[sample])))^2
      U.dot.g=U.dot.g+t(x.center[[i]])%*%x.center[[i]]
      if(0){
        for (j in 1:ny[i]) {
          for (k in 1:nx[i]) {
            K=K1(tx[[i]][k]-ty[[i]][j], h)
            U1=U1+K*x[[i]][k]*(y[[i]][j]-g.hat[sample]*z[[i]][j])
            U2=U2+K*x[[i]][k]^2
          }
        }
      }
    }
    #b.hat[sample]=U1/U2
    se.b.2sg[sample]=U.sq.g^.5/U.dot.g
    
  }
  if(sample %% 50 == 0) print(sample)
  
  sample=sample+1
  if(sample>Nsample) break
}

proc.time()-start

if(0){
  #naive
  apply(b.hat, 2, mean)-b
  apply(b.hat, 2, sd)
  apply(se.b, 2, mean)
  sum(b>b.hat[,2]-1.96*se.b[,2] & b<b.hat[,2]+1.96*se.b[,2])/Nsample
  
  #PLM
  mean(b.hat.plm)-b
  sd(b.hat.plm)
  mean(se.b.plm)
  plot(t.uls,a.hat.plm[[sample-1]])
  sum(b>b.hat.plm-1.96*se.b.plm & b<b.hat.plm+1.96*se.b.plm)/Nsample
  
  #Two stage
  mean(b.hat.2sg)-b
  sd(b.hat.2sg)
  mean(se.b.2sg)
  sum(b>b.hat.2sg-1.96*se.b.2sg & b<b.hat.2sg+1.96*se.b.2sg)/Nsample
}

b.sum=rep(0, 13)
b.sum[1]=n
b.sum[2]=apply(b.hat, 2, mean)[2]-b
b.sum[3]=apply(b.hat, 2, sd)[2]
b.sum[4]=apply(se.b, 2, mean)[2]
b.sum[5]=round(sum(b>b.hat[,2]-1.96*se.b[,2] & b<b.hat[,2]+1.96*se.b[,2])/Nsample*100)
b.sum[6]=mean(b.hat.plm)-b
b.sum[7]=sd(b.hat.plm)
b.sum[8]=mean(se.b.plm)
b.sum[9]=round(mean(b>b.hat.plm-1.96*se.b.plm & b<b.hat.plm+1.96*se.b.plm)*100)
b.sum[10]=mean(b.hat.2sg)-b
b.sum[11]=sd(b.hat.2sg)
b.sum[12]=mean(se.b.2sg)
b.sum[13]=round(mean(b>b.hat.2sg-1.96*se.b.2sg & b<b.hat.2sg+1.96*se.b.2sg)*100)

stargazer(b.sum, type='latex')
