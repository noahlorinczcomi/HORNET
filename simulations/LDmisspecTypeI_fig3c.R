pruning=function(bx,by,bxse,byse,R,S=diag(ncol(bx)),rho=0.25,robust=TRUE,intercept=TRUE) {
  n=nrow(as.matrix(bx)); p=ncol(bx)
  ps=apply(bx,1,function(h)1/(t(h)%*%h)) # not exactly P-values, but proportinoal to P-values
  ps=cbind(1:n,ps); ps=ps[order(ps[,2]),]
  ord=ps[,1]; ps=ps[,2]
  R0=R[ord,ord]; R0[lower.tri(R0,diag=TRUE)]=0
  drop=c()
  for(i in 1:n) {
    if(i %in% drop) next
    vR=R0[i,]
    w=which(abs(vR)>rho)
    if(length(w)==0) next else drop=c(drop,w)
  }
  if(sum(drop)>0) keep=c(1:n)[-ord[unique(drop)]] else keep=1:n
  bx=bx[keep,];by=by[keep];R0=R[keep,keep]
  if(intercept) bx=cbind(1,bx)
  sR=solve(R0)
  BTBi=solve(t(bx)%*%sR%*%bx)
  BTa=t(bx)%*%sR%*%by
  est=c(BTBi%*%BTa)
  if(robust) {
    D=diag(c(by-bx%*%est))
    se=sqrt(diag(BTBi%*%t(bx)%*%sR%*%D%*%R0%*%D%*%sR%*%bx%*%BTBi))
  } else {
    se=sqrt(diag(BTBi))
  }
  out=list(est=est,se=se,mstar=length(keep),keepinds=keep)
  return(out)
}
n=500
R=ar1(n,0.5)
ns=floor(seq(n+1,5000,length.out=100))
niter=1000
qsig=qchisq(1-0.05,1)
bx=mvnfast::rmvn(niter,rep(0,n),R)
res=c()
for(nval in 1:length(ns)) {
  dat=rWishart(niter,ns[nval],R)
  k=0
  issig=c()
  while(k!=niter) {
    k=k+1
    Rk=cov2cor(dat[,,k]); # W=solve(Rk)
    bxk=bx[k,]
    byk=rnorm(n)
    bxsek=rchisq(n,10000-1)*1/10/(10000-1)
    bysek=rchisq(n,10000-1)*1/10/(10000-1)
    fit=pruning(matrix(bxk),byk,matrix(bxsek),bysek,Rk,rho=0.3,robust=FALSE,intercept=FALSE)
    issig[k]=((fit$est/fit$se)^2)>qsig
  }
  plot(res,type='b')
  res[nval]=sum(issig)/niter
}
sm=lm(res~I(1:length(res))+I((1:length(res))^(2)))
xl='reference panel sample size'
yl='Type I error rate in eQTL-MR'
plot(ns,sm$fitted,type='b',lwd=2,col='royalblue')
library(ggplot2);library(dplyr)
library(ggplot2)
data.frame(x=ns,y=res) %>%
  ggplot(aes(x,y)) +
  stat_smooth(method='gam',level=0.999,fill='gray80',alpha=1,color='gray80') +
  geom_point(pch=4) +
  theme_bw() +
  labs(x='Size of LD reference panel',y='eQTL-MR false positive rate',
       title='Inflation in eQTL-MR') +
  theme(text=element_text(size=17)) +
  geom_hline(yintercept=0.05,lwd=1.2) +
  lims(y=c(0,0.3))