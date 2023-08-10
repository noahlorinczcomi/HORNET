datagen=function(n,p,propMiss=0.25,quant=0.4) {
  # propMiss is the proportion of missing on the gene at the center of the group
  # quant controls how close genes are to each other (larger quant=closer genes)
  rr=1:n; qq=quantile(rr,probs=c(quant,1-quant)); 
  rr1=qq[1]; rr2=qq[2]; rrs=round(seq(rr1,rr2,length.out=p))
  nmiss=floor(propMiss*n); nnonmiss=n-nmiss
  DataTrue=DataMiss=matrix(nr=n,nc=p); 
  for(ii in 1:p) {
    dat=dnorm(1:n,rrs[ii],sd=p/1.5); # plot(dat) # see for example
    dat=dat/sd(dat)
    DataTrue[,ii]=dat
    keepinds=(rrs[ii]-nnonmiss/2+1):(rrs[ii]+nnonmiss/2)
    keepinds=floor(keepinds)
    nonkeep=c(1:n)[c(1:n)%in%keepinds==FALSE]
    dat[nonkeep]=NA
    DataMiss[,ii]=dat
  }
  Omega=is.na(DataMiss)
  list(DataMiss=DataMiss,DataTrue=DataTrue,Omega=Omega)
}

addME=function(dgOut,S,R,amt=1/sqrt(nrow(dgOut$DataMiss))) {
  n=nrow(dgOut$DataMiss); p=ncol(dgOut$DataMiss)
  if(sum(R==diag(n))==nrow(R)^2) {
    newDat=t(apply(dgOut$DataMiss,1,function(h) h+rmvn(1,(1:p)*0,amt*S)))# rmvnorm(1,sigma=amt*S)))
  } else {
    K=kronecker(S*amt,R)
    # newDat=mvtnorm::rmvnorm(1,sigma=K) # slow(er)
    newDat=rmvn(1,mu=rep(0,ncol(K)),sigma=K) # fast(er)
  }
  newDat=matrix(c(newDat),nr=n,nc=p)
  dgOut$DataMiss=dgOut$DataMiss+newDat
  dgOut
}

dataPlot=function(dgOut,...) {
  n=nrow(dgOut$DataMiss);p=ncol(dgOut$DataMiss)
  pal=paletteer::paletteer_c("grDevices::Zissou 1", p)
  addMat=matrix(c(1:p),n,p,byrow=TRUE)
  matplot(dgOut$DataMiss+addMat,yaxt='n',ylab="data",xlab="SNP base pair position",col=pal,pch=16,...)
  for(ii in 1:p) {
    inds=which(dgOut$Omega[,ii])
    points(inds,dgOut$DataTrue[inds,ii]+ii,pch=16,col="gray70")
  }
  legend("topleft",legend="missing",pch=16,col="gray70")
}

soft=function(x,lambda) sign(x)*ifelse((abs(x)-lambda)<0,0,abs(x)-lambda)

inner=function(A,lambda,Suu=nrow(A)*diag(ncol(A))*0,R=diag(nrow(A)),
               eps=1e-1,max.iter=100,A0=NA,max.rank=ncol(A)-1) {
  n=nrow(A);
  ed=eigen(R)
  Thsq=ed$vectors%*%diag(1/sqrt(ed$values))%*%t(ed$vectors)
  Rsq=ed$vectors%*%diag(sqrt(ed$values))%*%t(ed$vectors)
  Omega=is.na(A)
  deco=svd(Suu)
  if(!is.matrix(A0)) {A0=A; A0[Omega]=0}
  A0=Thsq%*%A0
  adj=sqrt(1-sum(diag(Suu))/sum(diag(t(A0)%*%A0)))
  # A0sv=svd(A0);A0sv$d=A0sv$d*adj; A0=A0sv$u%*%diag(A0sv$d)%*%t(A0sv$v)
  es=c()
  k=0;error=eps+1;rank=ncol(A)+1
  while(k<max.iter & error>eps) {
    k=k+1
    ### no adj for measurement error
    # ss=svd(A0);dk=soft(ss$d,lambda);rank=sum(dk>0)
    ### adj for measurement error (not sure if its best way to adjust)
    ss=svd(A0);dk=soft(ss$d*adj,lambda); rank=sum(dk>0)
    Ak=ss$u%*%diag(dk)%*%t(ss$v)
    Ak[!Omega]=A[!Omega]
    error=norm(A0-Ak,'f'); es[k]=error
    A0=Ak
  }
  Ak=Rsq%*%Ak
  Ak[!Omega]=A[!Omega] # after back-transformation
  list(Ak=Ak,errors=es,kiter=k,rank=rank)
}

likefun=function(x,Theta) -t(x)%*%Theta%*%x

softImp=function(A,lambda.max,R=diag(dim(A)[1]),Theta=diag(prod(dim(A))),
                 Suu=diag(ncol(A)),nLambda=20,max.rank=ncol(A)-1,...) {
  p=ncol(A); pinds=1:p
  lambdas=seq(0,lambda.max,length.out=nLambda)
  era=ranks=like=c();RES=list()
  for(ll in 1:nLambda) {
    res=inner(A=A,lambda=lambdas[ll],Suu=Suu,R=R)
    era[ll]=tail(res$errors,1)
    like[ll]=likefun(c(res$Ak),Theta)
    RES[[ll]]=res$Ak
    ranks[ll]=res$rank
  }
  candidates=which(ranks<=max.rank); 
  #op=candidates[which.min(era[candidates])]
  op=candidates[which.max(like[candidates])]
  if(length(candidates)==0) stop("no solution: Must increase max.rank argument")
  out=list(Ak=RES[[op]],lambda=lambdas[op],errors=era,ranks=ranks,likes=like)
  out
}

###
library(mvnfast);library(visdat)
ar1=function(n,rho){h=rho^(toeplitz(1:n)-1);diag(h)=1;h}
n=100;p=10;pmiss=0.75;nmiss=floor(pmiss*n)
S=1/toeplitz(1:p)
amt=sqrt(n)*1/n^(1/3)
Suu=S*amt
R=ar1(n,0.5)
ps=seq(0,0.75,length.out=30)
niter=100
K=kronecker(S,R)
Theta=solve(K)
beta=rep(0,p);beta[1:floor(p/3)]=0.2 # can evaluate both Type I and II with this
qsig=qchisq(1-0.05,1)
for(i in 1:length(ps)) {
 k=0
 res1=res2=c()
 while(k<niter) {
  data=datagen(n,p,propMiss=ps[i])
  data=addME(data,S,R,amt=amt)
  Xmiss=data$DataMiss
  if(nrow(Xmiss)<p) res1[k]=NA
  X0=data$DataMiss; X0[is.na(X0)]=0; lambda.max=median(svd(X0)$d)*3/2
  Ximp=softImp(data$DataMiss,lambda.max=lambda.max,R=R,Theta=Theta,Suu=Suu/2+0.5*diag(p))$Ak
  y=data$DataTrue%*%beta+rnorm(n)
  fit1=lm(y~Xmiss-1)
  fit2=lm(y~Ximp-1)
  sig1=(coef(fit1)^2/diag(vcov(fit1)))>qsig
  sig2=(coef(fit2)^2/diag(vcov(fit2)))>qsig
 }
 print(i)
}

dataPlot(data,main='Estimated Associations')



















