nonconsecutive=function(x,minsize=5) {
  res=1;
  for(i in 1:length(x)) {
    k=i;doskip=TRUE
    while(k<length(x) & doskip) {
      k=k+1
      if(x[k]<(x[k-1]+minsize)) next else doskip=FALSE
    }
    res=c(res,x[k])
  }
  unique(res)
}

blockFit=function(M,cs,penNBlocks=TRUE) {
  n=nrow(M);pen=c();
  for(i in 2:length(cs)) {
    ix1=cs[i-1];ix2=cs[i]
    b=M[ix1:ix2,ix1:ix2]^2
    bo=M[ix1:ix2,-c(ix1:ix2)]^2
    p1=prod(dim(b)); p2=prod(dim(bo))
    fin=norm(b,'f'); fout=norm(bo%*%t(bo),'f')
    if(sum(fout)==0) fout=1
    if(penNBlocks) lenpen=log(length(cs)) else lenpen=0
    pen[i-1]=log(fin)*log(p1)-log(fout)*log(p2)-lenpen
  }
  pen
}

blockify=function(M,maxBlocks=10,minsize=10,...) {
  ds=c(); for(i in 2:ncol(M)) ds[i-1]=log(det(M[1:i,1:i]))
  n=length(ds); D=cbind(1:n,abs(ds-lag(ds)))
  D=D[order(D[,2]),]
  Dtry=D[1:floor(maxBlocks*1.5),] # x1.5 bc some may be ~consecutive positions
  pens=c()
  for(i in 1:nrow(Dtry)) {
    csi=c(1,Dtry[1:i,1],n+1)
    csi=nonconsecutive(sort(csi),minsize)
    pens[i]=sum(blockFit(M,csi,...))
  }
  cs=Dtry[1:which.max(pens),1]+1
  cs=c(1,cs,n+1)
  nonconsecutive(sort(cs),minsize) # outputs block cutpoints
}

zeroifyNonBlocks=function(M,cs) {
  Mc=M*0
  for(i in 2:length(cs)) {
    ix=cs[i-1]:cs[i]
    Mc[ix,ix]=M[ix,ix]
  }
  Mc
}
