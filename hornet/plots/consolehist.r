consolehist=function(x,bins=10,size=10,range=c(min(x),max(x))) {
    x=x[x>=range[1] & x<=range[2]]
    x=na.omit(x)
    cs=seq(min(x),max(x),length.out=bins+1)
    n=c()
    for(i in 2:length(cs)) {
        xi=x[x>cs[i-1] & x<=cs[i]]
        n[i-1]=length(xi)
    }
    actualn=n # when rescaling, some small counts will be rounded to 0, but I don't want to forget that they originally were not 0
    n=(n-min(n))/(max(n)-min(n))*size
    n=ceiling(n)
    prefix=c()
    for(i in 1:bins) prefix[i]=paste0(round(cs[i],2), ' (', actualn[i], ')')
    nleftchars=nchar(prefix)
    startat=max(nleftchars)+1
    res=c()
    for(i in 1:bins) {
        space=paste(rep('@',n[i]),collapse='')
        nonspace=paste(rep('-',size-n[i]),collapse='')
        toadd=paste(c(space,nonspace),collapse='')
        res[i]=toadd
    }
    # all toadd's will be of length `size` and will start at the `startat` position. need to fill from the end of `prefix`
    ntofill=startat-nleftchars
    for(i in 1:bins) {
        addedspace=paste(rep(' ',startat-nchar(prefix[i])),collapse='')
        prefix[i]=paste(prefix[i],addedspace,collapse='')
        res[i]=paste(prefix[i],res[i],collapse='')
    }
    # print
    cat('\n')
    for(i in 1:bins) cat(res[i],'\n')
}