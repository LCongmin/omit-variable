rm(list=ls(all=TRUE)) 
library(stargazer)
library(Matrix)
options(warn=2)
DATA=read.csv('/Users/lichen/Documents/Mizzou/Dissertation 2 - confounder/Program/Cmv.csv')

diag.dateX=as.Date(as.character(DATA$HIV), '%m/%d/%y')
base.cd4date=as.Date(as.character(DATA$baseCD4date), '%m/%d/%y')
base.hivdate=as.Date(as.character(DATA$baseHIVdate), '%m/%d/%y')
cd4date=as.Date(as.character(DATA$lastCD4date), '%m/%d/%y')
hivdate=as.Date(as.character(DATA$lastHIVdate), '%m/%d/%y')
cd4count=log(DATA$lastCD4)
hivcount=log(DATA$lastHIV)
patient=as.numeric(DATA$CMVPATIENT)
race=(as.character(DATA$ETHNIC)=='BLACK')+0
gender=(as.character(DATA$GENDER)=='M')+0
age=difftime(as.Date(as.character(DATA$HIV), '%m/%d/%y'), as.Date(paste(DATA$monthbirth, DATA$daybirth, DATA$yearbirth, sep='/'), '%m/%d/%Y'), units='days')/365.25
visitage=DATA$visitage/365.2422

diag.ind=c(1,diff(patient))-!is.na(diag.dateX)
age=cbind(patient[!is.na(diag.dateX)], age[!is.na(diag.dateX)])
diag.date=cbind(patient[!is.na(diag.dateX)], diag.dateX[!is.na(diag.dateX)])
base.cd4date=cbind(patient[!is.na(diag.dateX)], base.cd4date[!is.na(diag.dateX)])
base.hivdate=cbind(patient[!is.na(diag.dateX)], base.hivdate[!is.na(diag.dateX)])

na.cd4=!(is.na(cd4date) | is.na(cd4count)) & (as.character(DATA$ETHNIC)=='BLACK' | as.character(DATA$ETHNIC)=='WHITE') & !diag.ind #& !is.na(DATA$age)
cd4date=cd4date[na.cd4]
cd4count=cd4count[na.cd4]
cd4patient=patient[na.cd4]
date.ind=(cd4date>='1997-1-1')
cd4date=cd4date[date.ind]
cd4count=cd4count[date.ind]
cd4patient=cd4patient[date.ind]
na.hiv=!(is.na(hivdate) | is.na(hivcount)) & (as.character(DATA$ETHNIC)=='BLACK' | as.character(DATA$ETHNIC)=='WHITE') & !diag.ind #& !is.na(DATA$age)
hivdate=hivdate[na.hiv]
visitage=visitage[na.hiv]
hivcount=hivcount[na.hiv]
hivpatient=patient[na.hiv]
gender=gender[na.hiv]
race=race[na.hiv]
#age=age[na.hiv]
date.ind=(hivdate>='1997-1-1')
hivdate=hivdate[date.ind]
visitage=visitage[date.ind]
hivcount=hivcount[date.ind]
hivpatient=hivpatient[date.ind]
gender=gender[date.ind]
race=race[date.ind]
#age=age[date.ind]



ty=list()
tz=list()
x=list()
y=list()
z=list()

K1=function(x,h){
  0.75*(1-(x/h)^2)/h*(abs(x)<h)
}
K2=function(x,y,h1,h2){
  0.5625*(1-(x/h1)^2)*(1-(y/h2)^2)/h1/h2*(abs(x)<h1)*(abs(y)<h2)
}

LOCF=function(x,y,z,ty,tz){
  n=length(ty)
  ind=rep(T, n)
  z.LOCF=y
  for(i in 1:n){
    for(j in 1:length(ty[[i]])){
      if(sum(tz[[i]]<ty[[i]][j])!=0){
        z.LOCF[[i]][j]=z[[i]][tz[[i]]==max(tz[[i]][tz[[i]]<=ty[[i]][j]])]
      }else{
        z.LOCF[[i]][j]=NA
        y[[i]][j]=NA
        x[[i]][j]=NA
        ty[[i]][j]=NA
      }
    }
    z.LOCF[[i]]=z.LOCF[[i]][!is.na(z.LOCF[[i]])]
    y[[i]]=y[[i]][!is.na(y[[i]])]
    x[[i]]=x[[i]][!is.na(x[[i]][,1]),]
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

b.hat.lvcf=NA
g.hat.lvcf=b.hat.lvcf
a.hat.lvcf=b.hat.lvcf
b.hat.2sg=NA
g.hat.2sg.lvcf=b.hat.2sg
a.hat.2sg.lvcf=b.hat.2sg
g.hat.2sg.ks=b.hat.2sg
a.hat.2sg.ks=b.hat.2sg
b.hat.jrssb=NA
g.hat.jrssb=b.hat.jrssb
a.hat.jrssb=b.hat.jrssb
b.se.lvcf=NA
g.se.lvcf=b.se.lvcf
a.se.lvcf=b.se.lvcf
b.se.2sg=NA
g.se.2sg.lvcf=b.se.2sg
a.se.2sg.lvcf=b.se.2sg
g.se.2sg.ks=b.se.2sg
a.se.2sg.ks=b.se.2sg
b.se.jrssb=NA
g.se.jrssb=b.se.jrssb
a.se.jrssb=b.se.jrssb

sample=1
na.count=0

ind=1
cd4patient.ls=list()
hivpatient.ls=list()
tz0=list()
ty0=list()
for(i in 1:max(c(cd4patient, hivpatient))){
  if(sum(age[,1]==i)==0) next
  if(sum(cd4count[cd4patient==i])!=0 & sum(hivcount[hivpatient==i])!=0 & !is.na(age[age[,1]==i,2])
     ) {
    z[[ind]]=cd4count[cd4patient==i]
    tz[[ind]]=as.numeric(cd4date[cd4patient==i])
    z[[ind]]=z[[ind]][!c(0, diff(tz[[ind]])<=0)]
    tz[[ind]]=tz[[ind]][!c(0, diff(tz[[ind]])<=0)]-diag.date[diag.date[,1]==i, 2]
    y[[ind]]=hivcount[hivpatient==i]
    ty[[ind]]=as.numeric(hivdate[hivpatient==i])
    y[[ind]]=y[[ind]][!c(0, diff(ty[[ind]])<=0)]
    x[[ind]]=matrix(cbind(gender[hivpatient==i], race[hivpatient==i]), ncol=2)
    x[[ind]]=matrix(x[[ind]][!c(0, diff(ty[[ind]])<=0),], ncol=2)
    ty[[ind]]=ty[[ind]][!c(0, diff(ty[[ind]])<=0)]-diag.date[diag.date[,1]==i, 2]
    cd4patient.ls[[ind]]=rep(ind, length(z[[ind]]))
    hivpatient.ls[[ind]]=rep(ind, length(y[[ind]]))
    #if(length(tz[[ind]])!=length(z[[ind]])) break
    #if(sum(tz[[ind]]<0)+sum(ty[[ind]]<0)>0) break
    ind=ind+1
  }
}

############ outlier
z=z[-166]
tz=tz[-166]
y=y[-166]
ty=ty[-166]
x=x[-166]

plot(unlist(cd4patient.ls), unlist(tz), pch=19, ylab='visit time (day)', xlab='patient', col='#FF0000')
points(unlist(hivpatient.ls), unlist(ty), pch=17, col='#0000FF')
legend('bottomright', legend = c('CD4', 'HIV'), pch=c(19,17), col=c('#FF0000', '#0000FF'))

N=ind-2
Q3=quantile(c(unlist(tz),unlist(ty)))[4]
Q1=quantile(c(unlist(tz),unlist(ty)))[2]
HH=10
h.can=seq(N^(-.8), N^(-.6), length=HH)*2*(Q3-Q1)
fold=5
mse=matrix(0,nrow=fold, ncol=HH)
MSE=rep(0, HH)
Iter=20
Btsp=ceiling(N/2)

############## Naive ##############

X=cbind(1,matrix(unlist(sapply(x, t)), ncol=2, byrow=T))
Y=unlist(y)

b.hat.naive=solve(t(X)%*%X)%*%(t(X)%*%Y)

U.sq.g=0;U.dot.g=0
U.sq.b=0;U.g=0;U.b=0
U1=0;U2=0
for (i in 1:N) {
  x.temp=cbind(1, x[[i]])
  temp=t(x.temp)%*%(y[[i]]-x.temp%*%b.hat.naive)
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
b.se.naive=diag(solve(U.dot.g)%*%U.sq.g%*%solve(U.dot.g))^.5
2*pnorm(-abs(b.hat.naive/b.se.naive))

############## LVCF ##############

if(1){
  data=LOCF(x,y,z,ty,tz)
  X=data[[1]][data[[5]]]
  Y=data[[2]][data[[5]]]
  Z=data[[3]][data[[5]]]
  x.LOCF=matrix(unlist(sapply(X, t)), ncol=2, byrow=T)
  y.LOCF=unlist(Y)
  z.LOCF=unlist(Z)
  ty.LOCF=unlist(data[[4]][data[[5]]])
  #N=sum(x.LOCF!=0)
  n=length(X)

  Cov=cbind(1, z.LOCF, x.LOCF)
  #Cov=cbind(1,z.LOCF)

  est=solve(t(Cov)%*%Cov)%*%t(Cov)%*%y.LOCF
  b.hat.lvcf=est[3:4]
  g.hat.lvcf=est[2]
  a.hat.lvcf=est[1]
  
  U.sq.g=0;U.dot.g=0
  U.sq.b=0;U.g=0;U.b=0
  U1=0;U2=0
  for (i in 1:n) {
    Cov.temp=matrix(c(rep(1, length(Z[[i]])), Z[[i]], X[[i]]), ncol=4, byrow=F)
    #Cov.temp=cbind(1, Z[[i]])
    temp=t(Cov.temp)%*%(Y[[i]]-Cov.temp%*%est)
    U.sq.g=U.sq.g+temp%*%t(temp)
    U.dot.g=U.dot.g+t(Cov.temp)%*%Cov.temp
    if(0){
      for (j in 1:ny[i]) {
        for (k in 1:nx[i]) {
          K=K1(tx[[i]][k]-ty[[i]][j], h)
          U1=U1+K*x[[i]][k]*(y[[i]][j]-g.hat*z[[i]][j])
          U2=U2+K*x[[i]][k]^2
        }
      }
    }
  }
  #b.hat[sample]=U1/U2
  se=diag(solve(U.dot.g)%*%U.sq.g%*%solve(U.dot.g))^.5
  a.se.lvcf=se[1]
  b.se.lvcf=se[3:4]
  g.se.lvcf=se[2]
}
############## Centering + LVCF ##############

ny=sapply(y, length)

if(0){ # JRSSB method
  gr=sample(1:2, N, rep=T)
  b.can=rep(0,10);se.can=rep(0,10)
  for(hh in 1:10){
    h.temp=h.can[hh]
    x.uls=as.matrix(unlist(x))
    y.uls=as.matrix(unlist(y))
    ty.uls=as.matrix(unlist(ty))
    n=length(x.uls)
    S=S.matrix(ty.uls,h.temp)
    beta=as.numeric(solve(t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%x.uls)%*%
                      t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%y.uls)
    rsd=list()
    Z=list()
    for(i in 1:N){
      rsd[[i]]=y[[i]]-x[[i]]*beta
      Z[[i]]=cbind(1, z[[i]])
    }
    U=u(tz, ty, Z, rsd, h.temp, N)
    est=solve(U[[2]])%*%U[[1]]
    
    x1=x[gr==1]; y1=y[gr==1]; ty1=ty[gr==1]; z1=z[gr==1]; tz1=tz[gr==1]
    x.uls=as.matrix(unlist(x1))
    y.uls=as.matrix(unlist(y1))
    ty.uls=as.matrix(unlist(ty1))
    n=length(x.uls)
    S=S.matrix(ty.uls,h.temp)
    b1=as.numeric(solve(t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%x.uls)%*%
                    t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%y.uls)
    #if(is.na(b0)==1 | is.na(b1)==1) next
    rsd=list()
    Z=list()
    for(i in 1:sum(gr==1)){
      rsd[[i]]=y1[[i]]-x1[[i]]*b1
      Z[[i]]=cbind(1, z1[[i]])
    }
    U=u(tz1, ty1, Z, rsd, h.temp, sum(gr==1))
    est1=solve(U[[2]])%*%U[[1]]
    
    x0=x[gr==2]; y0=y[gr==2]; ty0=ty[gr==2]; z0=z[gr==2]; tz0=tz[gr==2]
    x.uls=as.matrix(unlist(x0))
    y.uls=as.matrix(unlist(y0))
    ty.uls=as.matrix(unlist(ty0))
    n=length(x.uls)
    S=S.matrix(ty.uls,h.temp)
    b0=as.numeric(solve(t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%x.uls)%*%
                    t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%y.uls)
    rsd=list()
    Z=list()
    for(i in 1:sum(gr==2)){
      rsd[[i]]=y0[[i]]-x0[[i]]*b0
      Z[[i]]=cbind(1, z0[[i]])
    }
    U=u(tz0, ty0, Z, rsd, h.temp, sum(gr==2))
    est0=solve(U[[2]])%*%U[[1]]
    
    b.can[hh]=est[1]
    se.can[hh]=(est1[1]-est0[1])^2/4
  }
  h.reg=h.can^2
  C=lm(b.can~h.can^2)$coef[2]
  MSE=C^2*h.can^4+se.can
  h=h.can[which.min(MSE)]
}

if(1){ # CV
  gr=sample(1:fold, N, replace=T)

  for(hh in 1:HH){
    h.temp=h.can[hh]
    for(cv in 1:fold){
        x0=x[gr!=cv]; y0=y[gr!=cv]; ty0=ty[gr!=cv]; z0=z[gr!=cv]; tz0=tz[gr!=cv]; ny0=ny[gr!=cv]
        x.center=x0
        y.center=y0
        x.uls=matrix(unlist(sapply(x0, t)), ncol=2, byrow=T)
        t.uls=unlist(ty0)
        y.uls=unlist(y0)
        Ex=list()
        Ey=list()
        for(i in 1:sum(gr!=cv)){
          Ex[[i]]=x0[[i]]
          Ey[[i]]=y0[[i]]
          for (k in 1:ny0[i]) {
            temp=K1(ty0[[i]][k]-t.uls, h.temp)/sum(K1(ty0[[i]][k]-t.uls, h.temp))
            Ex[[i]][k,]=t(x.uls)%*%temp
            Ey[[i]][k]=t(y.uls)%*%temp
          }
          x.center[[i]]=x0[[i]]-Ex[[i]]
          y.center[[i]]=y0[[i]]-Ey[[i]]
        }
        
        X=matrix(unlist(sapply(x.center, t)), ncol=2, byrow=T)
        Y=unlist(y.center)
        
        b0=solve(t(X)%*%X)%*%(t(X)%*%Y)
        
        rsd=list()
        Z=list()
        for(i in 1:sum(gr!=cv)){
          rsd[[i]]=y0[[i]]-x0[[i]]%*%b0
          Z[[i]]=cbind(1, z0[[i]])
        }
        U=u(tz0, ty0, Z, rsd, h.temp, sum(gr!=cv))
        est0=solve(U[[2]])%*%U[[1]]
        
        x1=x[gr==cv]; y1=y[gr==cv]; ty1=ty[gr==cv]; z1=z[gr==cv]; tz1=tz[gr==cv]
        mse.temp=0
        temp=0
        KK=0
        for(i in 1:sum(gr==cv)){
          for(j in 1:length(ty1[[i]])){
            for(k in 1:length(tz1[[i]])){
              KK=KK+K1(tz1[[i]][k]-ty1[[i]][j], h.temp)
              temp=temp+K1(tz1[[i]][k]-ty1[[i]][j], h.temp)*(y1[[i]][j]-est0[1]-est0[2]*z1[[i]][k]-x1[[i]][j,]%*%b0)^2
            }
          }
        }
        if(KK!=0) {mse.temp=mse.temp+temp/KK}
        mse[cv, hh]=mse.temp
      }
  }
  MSE=apply(mse, 2, mean)
  h=h.can[which.min(MSE)]
}


if(1){
  x.center=x
  y.center=y
  if(1){
    x.uls=matrix(unlist(sapply(x, t)), ncol=2, byrow=T)
    t.uls=unlist(ty)
    y.uls=unlist(y)
    Ex=list()
    Ey=list()
    for(i in 1:N){
      Ex[[i]]=x[[i]]
      Ey[[i]]=y[[i]]
      for (k in 1:ny[i]) {
        temp=K1(ty[[i]][k]-t.uls, h)/sum(K1(ty[[i]][k]-t.uls, h))
        Ex[[i]][k,]=t(x.uls)%*%temp
        Ey[[i]][k]=t(y.uls)%*%temp
      }
      x.center[[i]]=x[[i]]-Ex[[i]]
      y.center[[i]]=y[[i]]-Ey[[i]]
    }
  }
  
  X=matrix(unlist(sapply(x.center, t)), ncol=2, byrow=T)
  Y=unlist(y.center)
  
  b.hat.2sg=solve(t(X)%*%X)%*%(t(X)%*%Y)
  
  U.sq.g=0;U.dot.g=0
  U.sq.b=0;U.g=0;U.b=0
  U1=0;U2=0
  for (i in 1:N) {
    U.sq.g=U.sq.g+t(x.center[[i]])%*%(y.center[[i]]-x.center[[i]]%*%b.hat.2sg)%*%t(t(x.center[[i]])%*%(y.center[[i]]-x.center[[i]]%*%b.hat.2sg))
    U.dot.g=U.dot.g+t(x.center[[i]])%*%x.center[[i]]
    if(0){
      for (j in 1:ny[i]) {
        for (k in 1:nx[i]) {
          K=K1(tx[[i]][k]-ty[[i]][j], h)
          U1=U1+K*x[[i]][k]*(y[[i]][j]-g.hat*z[[i]][j])
          U2=U2+K*x[[i]][k]^2
        }
      }
    }
  }
  #b.hat[sample]=U1/U2
  b.se.2sg=diag(solve(U.dot.g)%*%U.sq.g%*%solve(U.dot.g))^.5

  if(1){
    rsd=list()
    for(i in 1:N){
      rsd[[i]]=y[[i]]-x[[i]]%*%b.hat.2sg
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
    g.hat.2sg.lvcf=est[2]
    a.hat.2sg.lvcf=est[1]
    
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
    a.se.2sg.lvcf=se[1]
    #b.se[sample]=se[3]
    g.se.2sg.lvcf=se[2]
  }


############## Centering + Ksmooth ##############
  
if(1){
    rsd=list()
    Z=list()
    for(i in 1:N){
      rsd[[i]]=y[[i]]-x[[i]]%*%b.hat.2sg
      Z[[i]]=cbind(1, z[[i]])
    }
    U=u(tz, ty, Z, rsd, h, N)
    est=solve(U[[2]])%*%U[[1]]
    g.hat.2sg.ks=est[2]
    a.hat.2sg.ks=est[1]
    g.temp=rep(NA, Iter)
    a.temp=rep(NA, Iter)
    next.ind=0
    iter=1
    repeat{
      ind=sample(N, size=Btsp, replace=T)
      X=x[ind]
      Y=y[ind]
      Z=z[ind]
      Ty=ty[ind]
      X.center=X
      Y.center=Y
        X.uls=matrix(unlist(sapply(X, t)), ncol=2, byrow=T)
        T.uls=unlist(Ty)
        Y.uls=unlist(Y)
        EX=list()
        EY=list()
        for(i in 1:Btsp){
          EX[[i]]=X[[i]]
          EY[[i]]=Y[[i]]
          for (k in 1:length(Y[[i]])) {
            temp=K1(Ty[[i]][k]-T.uls, h)/sum(K1(Ty[[i]][k]-T.uls, h))
            EX[[i]][k,]=t(X.uls)%*%temp
            EY[[i]][k]=t(Y.uls)%*%temp
          }
          X.center[[i]]=X[[i]]-EX[[i]]
          Y.center[[i]]=Y[[i]]-EY[[i]]
        }
      
      XX=matrix(unlist(sapply(X.center, t)), ncol=2, byrow=T)
      YY=unlist(Y.center)
      
      beta=solve(t(XX)%*%XX)%*%(t(XX)%*%YY)
      
      rsd.Btsp=list()
      Z.Btsp=list()
      for(i in 1:Btsp){
        rsd.Btsp[[i]]=Y[[i]]-X[[i]]%*%beta
        Z.Btsp[[i]]=cbind(1, Z[[i]])
      }
      tz.Btsp=tz[ind]
      ty.Btsp=ty[ind]
      
      U=u(tz.Btsp, ty.Btsp, Z.Btsp, rsd.Btsp, h, Btsp)
      est.temp=solve(U[[2]])%*%U[[1]]
      if(sum(is.na(est.temp)!=0)) {
        next.ind=1
      } else {
        g.temp[iter]=est.temp[2]
        a.temp[iter]=est.temp[1]
        iter=iter+1
      }
      if(iter>Iter) break
    }
    
    g.se.2sg.ks=sd(g.temp)*(Btsp/N)^.5
    a.se.2sg.ks=sd(a.temp)*(Btsp/N)^.5
  }
}


############## Cao JRSSB ##############

if(1){
  U.dot=0
  U=0
  for(i in 1:N){
    for(j in 1:length(ty[[i]])){
      for(k in 1:length(tz[[i]])){
        temp=c(1,z[[i]][k],x[[i]][j,])
        #temp=c(1,z[[i]][k])
        k1=K1(tz[[i]][k]-ty[[i]][j], h)
        U=U+k1*temp*y[[i]][j]
        U.dot=U.dot+k1*temp%*%t(temp)
      }
    }
  }
  est=solve(U.dot)%*%U
  a.hat.jrssb=est[1]
  b.hat.jrssb=est[3:4]
  g.hat.jrssb=est[2]
  U.sq=0
  for(i in 1:N){
    temp.U=0
    for(j in 1:length(ty[[i]])){
      for(k in 1:length(tz[[i]])){
        temp=c(1,z[[i]][k],x[[i]][j,])
        #temp=c(1,z[[i]][k])
        k1=K1(tz[[i]][k]-ty[[i]][j], h)
        temp.U=temp.U+k1*temp*as.numeric(y[[i]][j]-temp%*%est)
      }
    }
    U.sq=U.sq+temp.U%*%t(temp.U)
  }
  se=diag(solve(U.dot)%*%U.sq%*%solve(U.dot))^.5
  a.se.jrssb=se[1]
  b.se.jrssb=se[3:4]
  g.se.jrssb=se[2]
}

b.sum=matrix(0, nrow=2, ncol=15)
b.sum[,1]=b.hat.naive
b.sum[,2]=b.se.naive
b.sum[,3]=2*pnorm(-abs(b.hat.naive/b.se.naive))
b.sum[,4]=b.hat.lvcf
b.sum[,5]=b.se.lvcf
b.sum[,6]=2*pnorm(-abs(b.hat.lvcf/b.se.lvcf))
b.sum[,7]=b.hat.2sg
b.sum[,8]=b.se.2sg
b.sum[,9]=2*pnorm(-abs(b.hat.2sg/b.se.2sg))
b.sum[,10]=b.hat.2sg
b.sum[,11]=b.se.2sg
b.sum[,12]=2*pnorm(-abs(b.hat.2sg/b.se.2sg))
b.sum[,13]=b.hat.jrssb
b.sum[,14]=b.se.jrssb
b.sum[,15]=2*pnorm(-abs(b.hat.jrssb/b.se.jrssb))

g.sum=rep(0, 15)
g.sum[4]=g.hat.lvcf
g.sum[5]=g.se.lvcf
g.sum[6]=2*pnorm(-abs(g.hat.lvcf/g.se.lvcf))
g.sum[7]=g.hat.2sg.lvcf
g.sum[8]=g.se.2sg.lvcf
g.sum[9]=2*pnorm(-abs(g.hat.2sg.lvcf/g.se.2sg.lvcf))
g.sum[10]=g.hat.2sg.ks
g.sum[11]=g.se.2sg.ks
g.sum[12]=2*pnorm(-abs(g.hat.2sg.ks/g.se.2sg.ks))
g.sum[13]=g.hat.jrssb
g.sum[14]=g.se.jrssb
g.sum[15]=2*pnorm(-abs(g.hat.jrssb/g.se.jrssb))

a.sum=rep(0, 15)
a.sum[4]=a.hat.lvcf
a.sum[5]=a.se.lvcf
a.sum[6]=2*pnorm(-abs(a.hat.lvcf/a.se.lvcf))
a.sum[7]=a.hat.2sg.lvcf
a.sum[8]=a.se.2sg.lvcf
a.sum[9]=2*pnorm(-abs(a.hat.2sg.lvcf/a.se.2sg.lvcf))
a.sum[10]=a.hat.2sg.ks
a.sum[11]=a.se.2sg.ks
a.sum[12]=2*pnorm(-abs(a.hat.2sg.ks/a.se.2sg.ks))
a.sum[13]=a.hat.jrssb
a.sum[14]=a.se.jrssb
a.sum[15]=2*pnorm(-abs(a.hat.2sg.ks/a.se.2sg.ks))

stargazer(rbind(b.sum, g.sum, a.sum), type='latex', digits=3)
#MSE=MSE*h.can^.5
plot(h.can, MSE/N, type='o', xlab='Bandwidth', ylab='Squared prediction error', main='Squared prediction error against bandwidth')

#plot(h.can, mse[1,]/N,type='l', ylim=range(mse/N), main='Squred prediction error of each fold', ylab='Squared prediction error', xlab='Bandwidth'); lines(h.can, mse[2, ]/N); lines(h.can, mse[3, ]/N); lines(h.can, mse[4, ]/N); lines(h.can, mse[5, ]/N)
#hist(c(unlist(tz), unlist(ty)), xlab='Days after diagnosis', main='')
#hist(unlist(ty))

#plot(x[[1]][,1]-min(x[[1]][,1]), y[[1]]-y[[1]][which.min(x[[1]][,1])], type='l', xlim=c(0,4), ylim=c(-15,10))
#for(i in 2:N){
#  lines(x[[i]][,1]-min(x[[i]][,1]), y[[i]]-y[[i]][which.min(x[[i]][,1])])
#}




Tseq=seq(Q1, Q3, length=500)
Tseq=seq(1500, 3000, length=500)
yerr=list()
y.uls=unlist(y)
ty.uls=unlist(ty)
SST=rep(0, length(Tseq))
SSE=SST
ybar=SST
ysq=SST
BD=200

#for(ii in 1:length(Tseq)){
#  #for(i in 1:length(ty.uls)){
#  ybar[ii]=K1(Tseq[ii]-ty.uls,h)%*%y.uls/sum(K1(Tseq[ii]-ty.uls,h))
#  ysq[ii]=K1(Tseq[ii]-ty.uls,h)%*%y.uls^2/sum(K1(Tseq[ii]-ty.uls,h))
#  SST[ii]=ysq[ii]-ybar[ii]^2
#  #}
#}

B=b.hat.2sg; G=c(a.hat.2sg.ks, g.hat.2sg.ks)

for(ii in 1:length(Tseq)){
  KK=0
  temp=0
  for(i in 1:N){
    for(j in 1:length(y[[i]])){
      for(k in 1:length(z[[i]])){
        KK=KK+K2(ty[[i]][j]-Tseq[ii], tz[[i]][k]-Tseq[ii], h ,h)
        temp=temp+K2(ty[[i]][j]-Tseq[ii], tz[[i]][k]-Tseq[ii], h ,h)*(y[[i]][j]-G[1]-G[2]*z[[i]][k]-x[[i]][j,]%*%B)^2
      }
    }
  }
  if(KK!=0) {SSE[ii]=temp/KK}
}
#plot(Tseq, 1-SSE/SST, type='l')
#lines(Tseq, 1-SSE/SST, lty=2)

TT=seq(min(Tseq), max(Tseq), length=50)
SSEsm=TT
SSTsm=TT
for(ii in 1:length(TT)){
  SSEsm[ii]=K1(TT[ii]-Tseq, BD)%*%SSE/sum(K1(TT[ii]-Tseq, h))
#  SSTsm[ii]=K1(TT[ii]-Tseq, BD)%*%SST/sum(K1(TT[ii]-Tseq, h))
}
#SST=rep(mean(y.uls), 500)
#SSTsm=rep(mean(y.uls), 50)
SST=rep(sum((y.uls-mean(y.uls))^2)/length(y.uls),500)
SSTsm=rep(sum((y.uls-mean(y.uls))^2)/length(y.uls),50)

plot(TT, 1-SSEsm/SSTsm, type='l', col='#FF0000', ylab=expression(R^2), xlab='Days', ylim=c(-.1,.7))
print(mean(1-SSE/SST))

B=b.hat.lvcf; G=c(a.hat.lvcf, g.hat.lvcf)
#B=b.hat.lvcf; G=c(a.hat.2sg.lvcf, g.hat.2sg.lvcf)
for(ii in 1:length(Tseq)){
  KK=0
  temp=0
  for(i in 1:N){
    for(j in 1:length(y[[i]])){
      for(k in 1:length(z[[i]])){
        KK=KK+K2(ty[[i]][j]-Tseq[ii], tz[[i]][k]-Tseq[ii], h ,h)
        temp=temp+K2(ty[[i]][j]-Tseq[ii], tz[[i]][k]-Tseq[ii], h ,h)*(y[[i]][j]-G[1]-G[2]*z[[i]][k]-x[[i]][j,]%*%B)^2
      }
    }
  }
  if(KK!=0) {SSE[ii]=temp/KK}
}
#plot(Tseq, 1-SSE/SST, type='l')
#lines(Tseq, 1-SSE/SST, lty=2)

#TT=seq(min(Tseq), max(Tseq), length=50)
#SSEsm=TT
#SSTsm=TT
for(ii in 1:length(TT)){
  SSEsm[ii]=K1(TT[ii]-Tseq, BD)%*%SSE/sum(K1(TT[ii]-Tseq, h))
  #SSTsm[ii]=K1(TT[ii]-Tseq, BD)%*%SST/sum(K1(TT[ii]-Tseq, h))
}
lines(TT, 1-SSEsm/SSTsm, lty=1, col='#BBBB00')
print(mean(1-SSE/SST))

if(0){
  B=b.hat.jrssb; G=c(a.hat.jrssb, g.hat.jrssb)
  for(ii in 1:length(Tseq)){
    KK=0
    temp=0
    for(i in 1:N){
      for(j in 1:length(y[[i]])){
        for(k in 1:length(z[[i]])){
          KK=KK+K2(ty[[i]][j]-Tseq[ii], tz[[i]][k]-Tseq[ii], h ,h)
          temp=temp+K2(ty[[i]][j]-Tseq[ii], tz[[i]][k]-Tseq[ii], h ,h)*(y[[i]][j]-G[1]-G[2]*z[[i]][k]-x[[i]][j,]%*%B)^2
        }
      }
    }
    if(KK!=0) {SSE[ii]=temp/KK}
  }
  #plot(Tseq, 1-SSE/SST, type='l')
  #lines(Tseq, 1-SSE/SST, lty=2)
  
  #TT=seq(min(Tseq), max(Tseq), length=50)
  #SSEsm=TT
  #SSTsm=TT
  for(ii in 1:length(TT)){
    SSEsm[ii]=K1(TT[ii]-Tseq, BD)%*%SSE/sum(K1(TT[ii]-Tseq, h))
    #SSTsm[ii]=K1(TT[ii]-Tseq, BD)%*%SST/sum(K1(TT[ii]-Tseq, h))
  }
  lines(TT, 1-SSEsm/SSTsm, lty=3)
  print(mean(1-SSE/SST))
}


B=b.hat.naive
for(ii in 1:length(Tseq)){
  KK=0
  temp=0
  for(i in 1:N){
    for(j in 1:length(y[[i]])){
      #for(k in 1:length(z[[i]])){
      KK=KK+K1(ty[[i]][j]-Tseq[ii], h)
      temp=temp+K1(ty[[i]][j]-Tseq[ii] ,h)*(y[[i]][j]-B[1]-x[[i]][j,]%*%B[-1])^2
      #}
    }
  }
  if(KK!=0) {SSE[ii]=temp/KK}
}
#plot(Tseq, 1-SSE/SST, type='l')
#lines(Tseq, 1-SSE/SST, lty=2)

#TT=seq(min(Tseq), max(Tseq), length=50)
#SSEsm=TT
#SSTsm=TT
for(ii in 1:length(TT)){
  SSEsm[ii]=K1(TT[ii]-Tseq, BD)%*%SSE/sum(K1(TT[ii]-Tseq, h))
  #SSTsm[ii]=K1(TT[ii]-Tseq, BD)%*%SST/sum(K1(TT[ii]-Tseq, h))
}
lines(TT, 1-SSEsm/SSTsm, lty=1, col='#559955')
print(mean(1-SSE/SST))



legend('topright', legend = c('Two-stage', 'LVCF', 'PLM (CD4 missing)', 'Naive'), lty=1, col=c('#FF0000', '#BBBB00', '#0000FF', '#559955'))



x.uls=matrix(unlist(sapply(x, t)), ncol=2, byrow=T)
y.uls=as.matrix(unlist(y))
ty.uls=as.matrix(unlist(ty))
n=length(y.uls)

S=S.matrix(ty.uls,h)
S[is.na(S)]=0
beta=solve(t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%x.uls)%*%
  t(x.uls)%*%t(diag(1,n)-S)%*%(diag(1,n)-S)%*%y.uls
b.plm=beta
#b.plm=b.hat.2sg
a.plm=S%*%(y.uls-x.uls%*%b.plm)
for(ii in 1:length(Tseq)){
  KK=0
  temp=0
  ind=0
  for(i in 1:N){
    a.temp=a.plm[(ind+1):(ind+length(ty[[i]]))]
    ind=ind+length(ty[[i]])
    for(j in 1:length(y[[i]])){
      #for(k in 1:length(z[[i]])){
      KK=KK+K1(ty[[i]][j]-Tseq[ii], h)
      temp=temp+K1(ty[[i]][j]-Tseq[ii] ,h)*(y[[i]][j]-a.temp[j]-x[[i]][j,]%*%b.plm)^2
      #}
    }
  }
  if(KK!=0) {SSE[ii]=temp/KK}
}
#plot(Tseq, 1-SSE/SST, type='l')
#lines(Tseq, 1-SSE/SST, lty=2)

#TT=seq(min(Tseq), max(Tseq), length=50)
#SSEsm=TT
#SSTsm=TT
for(ii in 1:length(TT)){
  SSEsm[ii]=K1(TT[ii]-Tseq, BD)%*%SSE/sum(K1(TT[ii]-Tseq, h))
  #SSTsm[ii]=K1(TT[ii]-Tseq, BD)%*%SST/sum(K1(TT[ii]-Tseq, h))
}
#plot(TT, 1-SSEsm/SSTsm, type='l', ylab=expression(R^2), xlab='Days')
lines(TT, 1-SSEsm/SSTsm, lty=1, col='#0000FF')
print(mean(1-SSE/SST))

