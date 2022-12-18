##############################################################################
#### automatic bandwidth selection only applies to KS in PLM+KS and JRSSB ####
##############################################################################

rm(list=ls(all=TRUE)) 
#setwd("/Users/lichen/Documents/Mizzou/Dissertation 2 - confounder/Program/")
setwd('C:/Users/lcgd4/iCloudDrive/Documents/Mizzou/Dissertation 2 - confounder/Program')

library(MASS)
library(mvtnorm)
library(Matrix)
library(stargazer)
#library(Rcpp)

#sourceCpp("cpart.cpp")

####SIMULATE THE DATA

start=proc.time()

N=100
Nsample=1000
b=2
g=-1
a=1
EX=3
h.can=seq(N^(-.8), N^(-.6), length=10)
#.can=seq(.001, .1, length=20)
MSE=matrix(0, nrow=Nsample, ncol=10)
h=rep(0, Nsample)
Iter=20
Btsp=ceiling(N/2)

ty=list()
tz=list()
corr.x=list()
x=list()
y=list()
z=list()
corr.e=list()

K1=function(x,h){
  0.75*(1-(x/h)^2)/h*(abs(x)<h)
}

LOCF=function(x,y,z,ty,tz){
  n=length(ty)
  ind=rep(T, n)
  z.LOCF=y
  for(i in 1:n){
    for(j in 1:length(ty[[i]])){
      if(sum(tz[[i]]<ty[[i]][j])!=0){
        z.LOCF[[i]][j]=z[[i]][tz[[i]]==max(tz[[i]][tz[[i]]<ty[[i]][j]])]
      }else{
        z.LOCF[[i]][j]=NA
        y[[i]][j]=NA
        x[[i]][j]=NA
        ty[[i]][j]=NA
      }
    }
    z.LOCF[[i]]=z.LOCF[[i]][!is.na(z.LOCF[[i]])]
    y[[i]]=y[[i]][!is.na(y[[i]])]
    x[[i]]=x[[i]][!is.na(x[[i]])]
    ty[[i]]=ty[[i]][!is.na(ty[[i]])]
    if(length(x[[i]])==0) ind[i]=F
  }
  return(list(x,y,z.LOCF,ty, ind))
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

u=function(tx, ty, x, y, h, n){
  U.dot=0
  U=0
  for(i in 1:n){
    for(j in 1:length(ty[[i]])){
      for(k in 1:length(tx[[i]])){
        k1=K1(tx[[i]][k]-ty[[i]][j], h)
        U=U+k1*x[[i]][k,]*y[[i]][j]
        U.dot=U.dot+k1*x[[i]][k,]%*%t(x[[i]][k,])
      }
    }
  }
  return(list(U,U.dot))
}

b.hat.lvcf=rep(NA,Nsample)
g.hat.lvcf=b.hat.lvcf
a.hat.lvcf=b.hat.lvcf
b.hat.2sg=rep(NA,Nsample)
g.hat.2sg.lvcf=b.hat.2sg
a.hat.2sg.lvcf=b.hat.2sg
g.hat.2sg.ks=b.hat.2sg
a.hat.2sg.ks=b.hat.2sg
b.hat.jrssb=rep(NA, Nsample)
g.hat.jrssb=b.hat.jrssb
a.hat.jrssb=b.hat.jrssb
b.se.lvcf=rep(NA,Nsample)
g.se.lvcf=b.se.lvcf
a.se.lvcf=b.se.lvcf
b.se.2sg=rep(NA,Nsample)
g.se.2sg.lvcf=b.se.2sg
a.se.2sg.lvcf=b.se.2sg
g.se.2sg.ks=b.se.2sg
a.se.2sg.ks=b.se.2sg
b.se.jrssb=rep(NA, Nsample)
g.se.jrssb=b.se.jrssb
a.se.jrssb=b.se.jrssb

sample=1
na.count=0

start=proc.time()

repeat{
  ny=rpois(N,5)+1  # # of y obs
  nz=rpois(N,5)+1  # # of x obs
  #Z=rbinom(N,1,.5)#*2-1
  for(i in 1:N){
    ty[[i]]=as.matrix(runif(ny[i]))    # y obs time
    tz[[i]]=as.matrix(runif(nz[i]))   # x obs time
    t.temp=rbind(tz[[i]],ty[[i]])
    n.temp=nz[i]+ny[i]
    corr=exp(-abs(rep(1,n.temp)%*%t(t.temp)-t.temp%*%t(rep(1,n.temp))))
    #corr.z=exp(-abs(rep(1,ny[i])%*%t(ty[[i]])-ty[[i]]%*%t(rep(1,ny[i]))))
    corr.e=2^(-abs(rep(1,n.temp)%*%t(t.temp)-t.temp%*%t(rep(1,n.temp))))
    #corr.x[[i]]=exp(-abs(rep(1,nx[i])%*%t(tx[[i]])-tx[[i]]%*%t(rep(1,nx[i]))))
    
    MX=t.temp^.5
    #MX=rep(1,nobs[i])
    
    #MZ=.5+t.temp
    MZ=rep(2, n.temp)
    #MZ=2*sin(2*pi*t.temp)
    #MZ=1-2*(t.temp-.5)^2
    #MZ=2*abs(2*t.temp-1)-1
    #MZ=(2*t.temp^2+1)^-.5
    
    x.temp=mvrnorm(1,MX,corr)
    z.temp=mvrnorm(1,MZ, corr)
    #z.temp=rep(Z[i],n.temp)
    z[[i]]=as.matrix(z.temp[1:nz[i]])
    x[[i]]=as.matrix(x.temp[-(1:nz[i])])
    y.temp=a+g*z.temp+x.temp*b+as.matrix(mvrnorm(1,rep(0,n.temp),corr.e))
    y[[i]]=as.matrix(y.temp[-(1:nz[i])])
  }
  
  
  ############## LVCF ##############
  
  if(1){
    data=LOCF(x,y,z,ty,tz)
    X=data[[1]][data[[5]]]
    Y=data[[2]][data[[5]]]
    Z=data[[3]][data[[5]]]
    x.LOCF=unlist(X)
    y.LOCF=unlist(Y)
    z.LOCF=unlist(Z)
    ty.LOCF=unlist(data[[4]][data[[5]]])
    #N=sum(x.LOCF!=0)
    n=length(X)
    
    Cov=cbind(1, z.LOCF, x.LOCF)
    #Cov=Cov[x.LOCF!=0,]
    
    est=solve(t(Cov)%*%Cov)%*%t(Cov)%*%y.LOCF
    b.hat.lvcf[sample]=est[3]
    g.hat.lvcf[sample]=est[2]
    a.hat.lvcf[sample]=est[1]
    
    U.sq.g=0;U.dot.g=0
    U.sq.b=0;U.g=0;U.b=0
    U1=0;U2=0
    for (i in 1:n) {
      Cov.temp=cbind(1, Z[[i]], X[[i]])
      temp=t(Cov.temp)%*%(Y[[i]]-Cov.temp%*%est)
      U.sq.g=U.sq.g+temp%*%t(temp)
      U.dot.g=U.dot.g+t(Cov.temp)%*%Cov.temp
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
    se=diag(solve(U.dot.g)%*%U.sq.g%*%solve(U.dot.g))^.5
    a.se.lvcf[sample]=se[1]
    b.se.lvcf[sample]=se[3]
    g.se.lvcf[sample]=se[2]
  }


############## Centering + LVCF ##############
  
  gr=sample(1:5, N, replace=T)
  mse=matrix(0,nrow=5, ncol=10)
  for(hh in 1:10){
    h.temp=h.can[hh]
    for(cv in 1:5){
      x0=x[gr!=cv]; y0=y[gr!=cv]; ty0=ty[gr!=cv]; z0=z[gr!=cv]; tz0=tz[gr!=cv]; ny0=ny[gr!=cv]
      x.center=x0
      y.center=y0
      if(1){
        x.uls=unlist(x0)
        t.uls=unlist(ty0)
        y.uls=unlist(y0)
        Ex=list()
        Ey=list()
        for(i in 1:sum(gr!=cv)){
          Ex[[i]]=x0[[i]]
          Ey[[i]]=y0[[i]]
          for (k in 1:ny0[i]) {
            temp=K1(ty0[[i]][k]-t.uls, h.temp)/sum(K1(ty0[[i]][k]-t.uls, h.temp))
            Ex[[i]][k]=t(x.uls)%*%temp
            Ey[[i]][k]=t(y.uls)%*%temp
          }
          x.center[[i]]=x0[[i]]-Ex[[i]]
          y.center[[i]]=y0[[i]]-Ey[[i]]
        }
      }
      
      X=unlist(x.center)
      Y=unlist(y.center)
      
      b0=as.numeric(t(X)%*%Y/(t(X)%*%X))
      
      rsd=list()
      Z=list()
      for(i in 1:sum(gr!=cv)){
        rsd[[i]]=y0[[i]]-x0[[i]]*b0
        Z[[i]]=cbind(1, z0[[i]])
      }
      U=u(tz0, ty0, Z, rsd, h.temp, sum(gr!=cv))
      est0=solve(U[[2]])%*%U[[1]]
      
      x1=x[gr==cv]; y1=y[gr==cv]; ty1=ty[gr==cv]; z1=z[gr==cv]; tz1=tz[gr==cv]
      mse.temp=0
      KK=0
      temp=0
      for(i in 1:sum(gr==cv)){
        for(j in 1:length(ty1[[i]])){
          for(k in 1:length(tz1[[i]])){
            KK=KK+K1(tz1[[i]][k]-ty1[[i]][j], h.temp)
            temp=temp+K1(tz1[[i]][k]-ty1[[i]][j], h.temp)*(y1[[i]][j]-est0[1]-est0[2]*z1[[i]][k]-b0*x1[[i]][j])^2
          }
        }
      }
      if(KK!=0) {mse.temp=mse.temp+temp/KK}
      mse[cv, hh]=mse.temp
    }
    
    if(0){
      x.uls=as.matrix(unlist(x))
      y.uls=as.matrix(unlist(y))
      ty.uls=as.matrix(unlist(ty))
      n=length(x.uls)
      S=S.matrix(ty.uls,h)
      beta=as.numeric(solve(t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%x.uls)%*%
                        t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%y.uls)
      rsd=list()
      Z=list()
      for(i in 1:N){
        rsd[[i]]=y[[i]]-x[[i]]*beta
        Z[[i]]=cbind(1, z[[i]])
      }
      U=u(tz, ty, Z, rsd, h, N)
      est=solve(U[[2]])%*%U[[1]]
      
      x1=x[gr==1]; y1=y[gr==1]; ty1=ty[gr==1]; z1=z[gr==1]; tz1=tz[gr==1]
      x.uls=as.matrix(unlist(x1))
      y.uls=as.matrix(unlist(y1))
      ty.uls=as.matrix(unlist(ty1))
      n=length(x.uls)
      S=S.matrix(ty.uls,h)
      b1=as.numeric(solve(t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%x.uls)%*%
                      t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%y.uls)
      
      x0=x[gr==0]; y0=y[gr==0]; ty0=ty[gr==0]; z0=z[gr==0]; tz0=tz[gr==0]
      x.uls=as.matrix(unlist(x0))
      y.uls=as.matrix(unlist(y0))
      ty.uls=as.matrix(unlist(ty0))
      n=length(x.uls)
      S=S.matrix(ty.uls,h)
      b0=as.numeric(solve(t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%x.uls)%*%
                      t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%y.uls)
      
      #if(is.na(b0)==1 | is.na(b1)==1) next
      rsd=list()
      Z=list()
      for(i in 1:sum(gr==1)){
        rsd[[i]]=y1[[i]]-x1[[i]]*b1
        Z[[i]]=cbind(1, z1[[i]])
      }
      U=u(tz1, ty1, Z, rsd, h, sum(gr==1))
      est1=solve(U[[2]])%*%U[[1]]
      
      rsd=list()
      Z=list()
      for(i in 1:sum(gr==0)){
        rsd[[i]]=y0[[i]]-x0[[i]]*b0
        Z[[i]]=cbind(1, z0[[i]])
      }
      U=u(tz0, ty0, Z, rsd, h, sum(gr==0))
      est0=solve(U[[2]])%*%U[[1]]
      
      b.can[hh]=beta
      se.can[hh]=(b1-b0)^2/4
    }
    
  }
  #h.reg=h.can^2
  #C=lm(b.can~h.can^2)$coef[2]
  #temp=C^2*h.can^4+se.can
  #h[sample]=h.can[which.min(temp)]
  MSE[sample,]=apply(mse, 2, mean)
  h[sample]=h.can[which.min(MSE[sample,])]
  
  if(1){
    x.center=x
    y.center=y
    if(1){
      x.uls=unlist(x)
      t.uls=unlist(ty)
      y.uls=unlist(y)
      Ex=list()
      Ey=list()
      for(i in 1:n){
        Ex[[i]]=x[[i]]
        Ey[[i]]=y[[i]]
        for (k in 1:ny[i]) {
          temp=K1(ty[[i]][k]-t.uls, h[sample])/sum(K1(ty[[i]][k]-t.uls, h[sample]))
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
    for (i in 1:N) {
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
    b.se.2sg[sample]=U.sq.g^.5/U.dot.g
    
    if(1){
      rsd=list()
      for(i in 1:N){
        rsd[[i]]=y[[i]]-x[[i]]*b.hat.2sg[sample]
      }
      data=LOCF(x,rsd,z,ty,tz)
      Rsd=data[[2]][data[[5]]]
      Z=data[[3]][data[[5]]]
      Rsd.LOCF=unlist(Rsd)
      z.LOCF=unlist(Z)
      ty.LOCF=unlist(data[[4]][data[[5]]])
      #N=sum(x.LOCF!=0)
      n=length(Z)
      
      Cov=cbind(1, z.LOCF)
      #Cov=Cov[x.LOCF!=0,]
      
      est=solve(t(Cov)%*%Cov)%*%t(Cov)%*%Rsd.LOCF
      #b.hat[sample]=est[3]
      g.hat.2sg.lvcf[sample]=est[2]
      a.hat.2sg.lvcf[sample]=est[1]
      
      U.sq.g=0;U.dot.g=0
      U.sq.b=0;U.g=0;U.b=0
      U1=0;U2=0
      for (i in 1:n) {
        Cov.temp=cbind(1, Z[[i]])
        temp=t(Cov.temp)%*%(Rsd[[i]]-Cov.temp%*%est)
        U.sq.g=U.sq.g+temp%*%t(temp)
        U.dot.g=U.dot.g+t(Cov.temp)%*%Cov.temp
      }
      #b.hat[sample]=U1/U2
      se=diag(solve(U.dot.g)%*%U.sq.g%*%solve(U.dot.g))^.5
      a.se.2sg.lvcf[sample]=se[1]
      #b.se[sample]=se[3]
      g.se.2sg.lvcf[sample]=se[2]
    }
    
 ############## Centering + Ksmooth ##############
    
    if(1){
      rsd=list()
      Z=list()
      for(i in 1:N){
        rsd[[i]]=y[[i]]-x[[i]]*b.hat.2sg[sample]
        Z[[i]]=cbind(1, z[[i]])
      }
      U=u(tz, ty, Z, rsd, h[sample], N)
      est=solve(U[[2]])%*%U[[1]]
      g.hat.2sg.ks[sample]=est[2]
      a.hat.2sg.ks[sample]=est[1]
      g.temp=rep(NA, Iter)
      a.temp=rep(NA, Iter)
      next.ind=0
      for(iter in 1:Iter){
        ind=sample(N, size=Btsp, replace=T)
        X=x[ind]
        Y=y[ind]
        Z=z[ind]
        Ty=ty[ind]
        X.center=X
        Y.center=Y
        X.uls=unlist(X)
        T.uls=unlist(Ty)
        Y.uls=unlist(Y)
        EX=list()
        EY=list()
        for(i in 1:Btsp){
          EX[[i]]=X[[i]]
          EY[[i]]=Y[[i]]
          for (k in 1:length(Y[[i]])) {
            temp=K1(Ty[[i]][k]-T.uls, h[sample])/sum(K1(Ty[[i]][k]-T.uls, h[sample]))
            EX[[i]][k,]=t(X.uls)%*%temp
            EY[[i]][k]=t(Y.uls)%*%temp
          }
          X.center[[i]]=X[[i]]-EX[[i]]
          Y.center[[i]]=Y[[i]]-EY[[i]]
        }
        
        XX=unlist(X.center)
        YY=unlist(Y.center)
        
        beta=as.numeric(solve(t(XX)%*%XX)%*%(t(XX)%*%YY))
        
        rsd.Btsp=list()
        Z.Btsp=list()
        for(i in 1:Btsp){
          rsd.Btsp[[i]]=Y[[i]]-X[[i]]*beta
          Z.Btsp[[i]]=cbind(1, Z[[i]])
        }
        tz.Btsp=tz[ind]
        ty.Btsp=ty[ind]
        
        U=u(tz.Btsp, ty.Btsp, Z.Btsp, rsd.Btsp, h[sample], Btsp)
        est.temp=solve(U[[2]])%*%U[[1]]
        if(sum(is.na(est.temp)!=0)) {
          next.ind=1
          break
        }
        g.temp[iter]=est.temp[2]
        a.temp[iter]=est.temp[1]
      }
      
      if(next.ind) next
      g.se.2sg.ks[sample]=sd(g.temp)*(Btsp/N)^.5
      a.se.2sg.ks[sample]=sd(a.temp)*(Btsp/N)^.5
    }
  }
  
  
############## Cao JRSSB ##############
  
  if(1){
    U.dot=0
    U=0
    for(i in 1:N){
      for(j in 1:length(ty[[i]])){
        for(k in 1:length(tz[[i]])){
          temp=c(1,x[[i]][j],z[[i]][k])
          k1=K1(tz[[i]][k]-ty[[i]][j], h[sample])
          U=U+k1*temp*y[[i]][j]
          U.dot=U.dot+k1*temp%*%t(temp)
        }
      }
    }
    est=solve(U.dot)%*%U
    a.hat.jrssb[sample]=est[1]
    b.hat.jrssb[sample]=est[2]
    g.hat.jrssb[sample]=est[3]
    U.sq=0
    for(i in 1:N){
      temp.U=0
      for(j in 1:length(ty[[i]])){
        for(k in 1:length(tz[[i]])){
          temp=c(1,x[[i]][j],z[[i]][k])
          k1=K1(tz[[i]][k]-ty[[i]][j], h[sample])
          temp.U=temp.U+k1*temp*as.numeric(y[[i]][j]-temp%*%est)
        }
      }
      U.sq=U.sq+temp.U%*%t(temp.U)
    }
    se=diag(solve(U.dot)%*%U.sq%*%solve(U.dot))^.5
    a.se.jrssb[sample]=se[1]
    b.se.jrssb[sample]=se[2]
    g.se.jrssb[sample]=se[3]
  }
  
  
  if(sample %% 50 == 0) print(sample)
  
  sample=sample+1
  if(sample>Nsample) break
  
}

proc.time()-start

b.sum=rep(0, 17)
b.sum[1]=N
b.sum[2]=mean(b.hat.lvcf)-b
b.sum[3]=sd(b.hat.lvcf)
b.sum[4]=mean(b.se.lvcf)
b.sum[5]=round(mean(b>b.hat.lvcf-1.96*b.se.lvcf & b<b.hat.lvcf+1.96*b.se.lvcf)*100)
b.sum[6]=mean(b.hat.2sg)-b
b.sum[7]=sd(b.hat.2sg)
b.sum[8]=mean(b.se.2sg)
b.sum[9]=round(mean(b>b.hat.2sg-1.96*b.se.2sg & b<b.hat.2sg+1.96*b.se.2sg)*100)
b.sum[14]=mean(b.hat.jrssb)-b
b.sum[15]=sd(b.hat.jrssb)
b.sum[16]=mean(b.se.jrssb)
b.sum[17]=round(mean(b>b.hat.jrssb-1.96*b.se.jrssb & b<b.hat.jrssb+1.96*b.se.jrssb)*100)

g.sum=rep(0, 17)
g.sum[1]=N
g.sum[2]=mean(g.hat.lvcf)-g
g.sum[3]=sd(g.hat.lvcf)
g.sum[4]=mean(g.se.lvcf)
g.sum[5]=round(mean(g>g.hat.lvcf-1.96*g.se.lvcf & g<g.hat.lvcf+1.96*g.se.lvcf)*100)
g.sum[6]=mean(g.hat.2sg.lvcf)-g
g.sum[7]=sd(g.hat.2sg.lvcf)
g.sum[8]=mean(g.se.2sg.lvcf)
g.sum[9]=round(mean(g>g.hat.2sg.lvcf-1.96*g.se.2sg.lvcf & g<g.hat.2sg.lvcf+1.96*g.se.2sg.lvcf)*100)
g.sum[10]=mean(g.hat.2sg.ks)-g
g.sum[11]=sd(g.hat.2sg.ks)
g.sum[12]=mean(g.se.2sg.ks)
g.sum[13]=round(mean(g>g.hat.2sg.ks-1.96*g.se.2sg.ks & g<g.hat.2sg.ks+1.96*g.se.2sg.ks)*100)
g.sum[14]=mean(g.hat.jrssb)-g
g.sum[15]=sd(g.hat.jrssb)
g.sum[16]=mean(g.se.jrssb)
g.sum[17]=round(mean(g>g.hat.jrssb-1.96*g.se.jrssb & g<g.hat.jrssb+1.96*g.se.jrssb)*100)

a.sum=rep(0, 17)
a.sum[1]=N
a.sum[2]=mean(a.hat.lvcf)-a
a.sum[3]=sd(a.hat.lvcf)
a.sum[4]=mean(a.se.lvcf)
a.sum[5]=round(mean(a>a.hat.lvcf-1.96*a.se.lvcf & a<a.hat.lvcf+1.96*a.se.lvcf)*100)
a.sum[6]=mean(a.hat.2sg.lvcf)-a
a.sum[7]=sd(a.hat.2sg.lvcf)
a.sum[8]=mean(a.se.2sg.lvcf)
a.sum[9]=round(mean(a>a.hat.2sg.lvcf-1.96*a.se.2sg.lvcf & a<a.hat.2sg.lvcf+1.96*a.se.2sg.lvcf)*100)
a.sum[10]=mean(a.hat.2sg.ks)-a
a.sum[11]=sd(a.hat.2sg.ks)
a.sum[12]=mean(a.se.2sg.ks)
a.sum[13]=round(mean(a>a.hat.2sg.ks-1.96*a.se.2sg.ks & a<a.hat.2sg.ks+1.96*a.se.2sg.ks)*100)
a.sum[14]=mean(a.hat.jrssb)-a
a.sum[15]=sd(a.hat.jrssb)
a.sum[16]=mean(a.se.jrssb)
a.sum[17]=round(mean(a>a.hat.jrssb-1.96*a.se.jrssb & a<a.hat.jrssb+1.96*a.se.jrssb)*100)

stargazer(rbind(b.sum, g.sum, a.sum), type='latex')
plot(h.can, apply(MSE,2,mean)/N, type='o', xlab='Bandwidth', ylab='Squared prediction error', main='Squared prediction error against bandwidth')

proc.time()-start