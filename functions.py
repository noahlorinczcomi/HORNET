
from scipy.stats import norm
import warnings
import scipy
import argparse
import os
import time
import platform # platform.system() -> 'Linux': Linux, 'Darwin': Mac, 'Windows': Windows
import re
import subprocess
import collections
from numpy import matrix
from scipy import stats
from statsmodels.regression.quantile_regression import QuantReg
import random
import pandas
import numpy
import math
from MNet import * # this is a source file in the HORNET directory

### optional turning off of warnings
# numpy.seterr(divide='ignore')

#def my_formatwarning(message, category, filename, lineno, line=None):
#  print(message, category)
#  # lineno is the line number you are looking for
#  print('file:', filename, 'line number:', lineno)

#warnings.formatwarning=my_formatwarning

##################################################################################################################################################
# so-called lower-level functions
def makeI(k): # function to define a kxk identity matrix
    I=numpy.zeros((k,k)); numpy.fill_diagonal(I,1)
    return I

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def unpackDictOfLists(ddict):
    els=[]
    for i in range(0,len(ddict)):
        els.append(list(ddict[list(ddict)[i]]))
    return flatten_list(els)

def quickAppend(nestedlist):
    out=[nestedlist[x].append(float('nan')) for x in range(0,len(nestedlist))]
    return out

def c(A): # function to vectorize (going down columns) a matrix A
    return t(numpy.asarray(A)).ravel() 

def whichMax(ll):
    return list(numpy.where(numpy.array(ll)==max(ll)))[0][0]

def subSamCalcR(Zmat,k,niter=1000,repl=True):
    # Rubins imputation standard errors for rs
    p=Zmat.shape[1];n=Zmat.shape[0]
    rs,vBetween,vWithin=numpy.zeros((p,p)),numpy.zeros((p,p)),numpy.zeros((p,p))
    muHat=numpy.array(Zmat.corr())
    for jj in range(0,niter):
        subSam=Zmat.iloc[:,:].sample(n=k,replace=repl,axis=0)
        cc=numpy.array(subSam.corr())
        rs=rs+cc
        vWithin=vWithin+((1-cc**2)**2)/n # see The Standard Deviation of the Correlation Coefficient, JASA
        vBetween=vBetween+(cc-muHat)**2
    rs=rs/niter; vWithin=vWithin/niter; vBetween=vBetween/niter;
    vTotal=vWithin+vBetween+vBetween/niter
    return rs,vTotal # correlation matrix estimate, estimate of variance of correlation matrix

# this one subsamples from within evenly-spaced blocks, assuming the data in Zmat are sorted according to BP
def blockSubSamCalcR(ZmatP,k,niter,verbose):
    if verbose==True:
        print("data must be sorted by base pair position")
    p=Zmat.shape[1];n=Zmat.shape[0]
    rs,vBetween,vWithin=numpy.zeros((p,p)),numpy.zeros((p,p)),numpy.zeros((p,p))
    blockInds=list(numpy.linspace(0,Zmat.shape[0]-1,k)); blockInds=[int(blockInds[x]) for x in range(0,len(blockInds))]
    muHat=Zmat.loc[blockInds,:].corr()
    for jj in range(0,niter):
        bbInds=[random.randint(blockInds[bb],blockInds[bb+1]) for bb in range(0, len(blockInds)-1)]
        subSam=Zmat.iloc[bbInds,:]
        cc=numpy.array(subSam.corr())
        rs=rs+cc
        vWithin=vWithin+((1-cc**2)**2)/n # see The Standard Deviation of the Correlation Coefficient, JASA
        vBetween=vBetween+(cc-muHat)**2
    rs=rs/niter; vWithin=vWithin/niter; vBetween=vBetween/niter;
    vTotal=vWithin+vBetween+vBetween/niter
    return rs,vTotal # correlation matrix estimate, estimate of variance of correlation matrix

def reshapeCorrelDf(merged,theseGenes,geneLabel,joinType):
    dd={};newColNames=[]
    for jj in range(0,len(theseGenes)):
        dfGene=merged[merged[geneLabel]==theseGenes[jj]]; dfGene=dfGene[['geneSNP','geneZ']];
        # rename geneZ to be gene-specific
        cn=theseGenes[jj]+'_Z'; newColNames.append(cn)
        dfGene=dfGene.rename(columns={'geneZ': cn})
        dd[theseGenes[jj]]=dfGene
        if jj>0:
            prevkey=list(dd)[jj-1]; prevdat=dd[prevkey]
            if joinType=='inner':
                binding=pandas.merge(prevdat,dfGene,how='inner',left_on='geneSNP',right_on='geneSNP')
            else:
                binding=pandas.merge(prevdat,dfGene,how='outer',left_on='geneSNP',right_on='geneSNP')
            dd[theseGenes[jj]]=binding
    # now add phenotype to this, first column
    toBindAgain=merged[['phenoSNP','phenoZ']].drop_duplicates(['phenoSNP','phenoZ'])
    binding0=pandas.merge(toBindAgain,binding,left_on='phenoSNP',right_on='geneSNP')
    mergedcut=merged.drop_duplicates(subset=['geneSNP'],keep='first')
    binding0=pandas.merge(mergedcut.drop(columns=['geneZ','phenoZ']),binding0,left_on='geneSNP',right_on='geneSNP')       
    return binding0

def flatten_list(_2d_list): # flatten list of lists into a single list
    flat_list=[]
    # Iterate through the outer list
    for element in _2d_list:
        if type(element) is list:
            # If the element is of type list, iterate through the sublist
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list

def subR(R,removeInds):
    a1=numpy.delete(R,removeInds,1)
    a2=numpy.delete(a1,removeInds,0)
    return a2

def SpleioP(Bhat,Ahat,UU,UV,VV,Thetahat,VThetahat):
    m=len(Bhat[:,0]);q=len(Ahat[0,:])
    I=makeI(q)
    Spleio=[]
    for i in range(0,m):
        K=numpy.kron(I,Bhat[i,:])
        num=Ahat[i,:]-Thetahat.T@Bhat[i,:]
        den=VV+Thetahat.T@UU@Thetahat+K@VThetahat@K.T-2*Thetahat.T@UV
        S=num@numpy.linalg.inv(den)@num.T; S=float([S][0])
        Spleio.append(S)
    P=[1-stats.chi2.cdf(Spleio[x],q) for x in range(0,m)]
    return P

def estMRBEE(Bhat,Ahat,UU,UV,VV,Rinv,subff=True):
    # Bhat,UU: exposure
    # Ahat,VV: outcome
    # Rinv: inverse of LD matrix
    m=len(Bhat[:,0]);q=len(Ahat[0]);p=Bhat.shape[1]
    est=numpy.linalg.inv(Bhat.T@Rinv@Bhat-m*UU)@(Bhat.T@Rinv@Ahat-m*UV)
    F=(Bhat.T@Rinv@Bhat-m*UU)
    res0=Ahat.squeeze()-(Bhat@est).squeeze()
    ResMat=numpy.zeros((len(res0),len(res0))); numpy.fill_diagonal(ResMat,res0.squeeze())
    if subff==True:
        ff=numpy.zeros((m,p))
        dd=flatten_list([l.tolist() for l in UV-UU@est])
        for i in range(0,ff.shape[0]):
            ff[i,:]=dd
        E=ResMat@Bhat-ff
    else:
        E=ResMat@Bhat
    V=E.T@Rinv@E
    Varest=numpy.linalg.inv(F)@V@numpy.linalg.inv(F)
    return est, Varest

def bootMRBEE(Bhat,Ahat,UU,UV,VV,R,niter=500):
    # R: LD matrix
    m=len(Bhat[:,0]); p=Bhat.shape[1]; ests=numpy.zeros((niter,p))
    for j in range(0,niter):
        ind=random.choices(range(0,m),k=m)
        kBhat=Bhat[ind,:];kAhat=Ahat[ind,:];
        kR1=R[ind,:];kR=kR1[:,ind]; kR=numpy.linalg.pinv(kR) # generalized inverse because not guaranteed to be posdef
        est,V=estMRBEE(kBhat,kAhat,UU,UV,VV,kR)
        estl=[l.tolist() for l in est]
        ests[j,:]=flatten_list(estl)
    cc=pandas.DataFrame(ests).cov()
    ests=numpy.array([numpy.mean(ests[:,x]) for x in range(0,bx.shape[1])]); cc=numpy.array(cc)
    return ests,cc

def imrbee(Bhat,Ahat,UU,UV,VV,Rinv,R,PleioPThreshold,boot=False,niter=1000,initial="bee",max_iter=15):
    warnings.filterwarnings('ignore')
    # R: LD matrix
    # Rinv: inverse of LD matrix
    m=len(Bhat[:,0]); p=len(Bhat[0,:])+1; q=len(Ahat[0,:]) # +1 on p bc I will add an intercept without asking
    Bhat_=numpy.column_stack(([1]*Bhat.shape[0],Bhat))
    UU_=numpy.column_stack(([0]*(p-1),UU)); UU_=numpy.row_stack(([0]*p, UU_))
    UV_=numpy.row_stack((0,UV))
    if initial=="robust":
        mod=QuantReg(Ahat,Bhat_) # includes an intercept
        res=mod.fit(q=0.5)
        est0=res.params.reshape((p,1))
        V0=res.cov_params()
        pleioPs0=SpleioP(Bhat_,Ahat,UU_,UV_,VV,est0,V0)
    else:
        est0,V0=estMRBEE(Bhat_,Ahat,UU_,UV_,VV,Rinv)
        pleioPs0=SpleioP(Bhat_,Ahat,UU_,UV_,VV,est0,V0)
    k=0; thetadiff=1; diffs=[]; kR=Rinv.copy()
    mask=numpy.ones((m,),dtype='bool')
    Outliers0=[i for i in range(0,len(pleioPs0)) if pleioPs0[i]<0.05] # the initial threshold can be stricter than the one used in iteration
    mask[Outliers0]=False
    Outliers0=numpy.array(Outliers0); Outliers=Outliers0.copy().tolist()
    while (k<max_iter) & (thetadiff>(0.0001*(p*q+p))) & (len(Outliers)<(m-10)):   
        k=k+1
        kBhat=Bhat_[mask]
        kAhat=Ahat[mask]
        Rsub=R[mask,:][:,mask]
        kR=numpy.linalg.pinv(Rsub)
        fitkEst,fitkVar=estMRBEE(kBhat,kAhat,UU_,UV_,VV,kR)
        pleioPs0=SpleioP(Bhat_,Ahat,UU_,UV_,VV,fitkEst,fitkVar)
        Outliers=[i for i in range(0,len(pleioPs0)) if pleioPs0[i]<PleioPThreshold]
        mask[Outliers]=0
        if len(Outliers)==kBhat.shape[0]:
            thetadiff=0
            est0=numpy.array([1]*p)*float('nan')
            V0=numpy.zeros((p,p))*float('nan')
            Outliers=list(range(0,m))
            kR=float('nan')
        elif len(Outliers)==0:
            thetadiff=0
            est0,V0,Outliers,kR=fitkEst,fitkVar,Outliers0,kR
        else:
            thetadiff=numpy.sum(abs(fitkEst-est0))
            est0=fitkEst;V0=fitkVar
            diffs.append(thetadiff)
    if (boot is True) & (len(Outliers)<m):
        if len(Outliers)==0:
            kBhat=Bhat.copy(); kAhat=Ahat.copy(); Rsub=R.copy(); kR=numpy.linalg.inv(Rsub)
        else:
            kBhat=numpy.delete(Bhat, (Outliers), axis=0)
            kAhat=numpy.delete(Ahat, (Outliers), axis=0)
            Rsub=subR(R,Outliers)
            kR=numpy.linalg.pinv(Rsub)
        est0,V0=bootMRBEE(kBhat,kAhat,UU,UV,VV,kR,niter=niter)
    warnings.filterwarnings('default')
    return est0, V0, Outliers, kR, k

def makeMask(n,inds,starting=True):
    mask=numpy.ones((n,),dtype='bool') if starting==True else numpy.zeros((n,),dtype='bool')
    mask[inds]=(starting==False)
    return mask

def IVW(bx,by,Theta,ogZs,uniIVPThreshold=5e-8,has_intercept=False):
    # ogZs: original Z-stats (in the case of a transformation of bx by LD structure)
    # Theta: inverse of LD matrix
    # this will perform univariable and multivariable MR
    # uniIVPThreshold: if the user wants to restrict univariable MR IVs to only those significant at some threshold (recall bx and by are Z-scores)
    m=bx.shape[0];p=bx.shape[1]
    if has_intercept==False:
        bx_=numpy.column_stack(([1]*m,bx))
        # MVMR
    BTB=bx_.T@Theta@bx_
    hasInverse=numpy.linalg.matrix_rank(BTB)==BTB.shape[0]
    if hasInverse:
        BTBinv=numpy.linalg.inv(BTB)
    else:
        BTBinv=numpy.linalg.pinv(BTB)
    mvmrest=BTBinv@bx_.T@Theta@by
    mvmrest=mvmrest.reshape((p+1,))
    res=by.reshape((m,))-bx_@mvmrest
    sig2=sum(res**2)/(m-p-1)
    mvmrse=numpy.diag(sig2*BTBinv)**0.5
    mvmrest=mvmrest[1:] # take off intercept-relevant terms
    mvmrse=mvmrse[1:] # take off intercept-relevant terms
    # univariable
    UniEst=[];UniSE=[];q0=norm.ppf(1-uniIVPThreshold)
    for _ in range(0,p):
        bxj_=bx[:,_]; bxj_=numpy.column_stack(([1]*m,bxj_))
        mask=numpy.abs(ogZs[:,_])>q0
        newm=sum(mask)
        if newm<5: # if <5 SNPs meet P-value threshold (some genes may have SNPs like this since we keep all SNPs in MRJones with P<tau in a chi-square joint test)
            UniEst.append(float('nan'))
            UniSE.append(float('nan'))
        else:
            bxj_=bxj_[mask]; by_=by[mask]; Theta_=Theta[mask,:][:,mask]
            btbj=numpy.linalg.inv(bxj_.T@Theta_@bxj_)
            uniest=btbj@bxj_.T@Theta_@by_
            uniest=uniest.reshape((2,))
            res=by_.reshape((newm,))-bxj_@uniest
            sig2=sum(res**2)/(newm-2)
            unise=numpy.diag(sig2*btbj)**0.5
            UniEst.append(uniest[1])
            UniSE.append(unise[1])
    dd={'IVW_MVMR_Est': mvmrest, 'IVW_MVMR_SE': mvmrse, 'IVW_UNIMR_Est': UniEst, 'IVW_UNIMR_SE': UniSE}
    pp=pandas.DataFrame(dd)
    return pp

def LDShrinkUpper(ld,tau,snps,missMethod='drop'):
    # if missMethod=='drop':
    #     ld=numpy.nan_to_num(ld, nan=1)
    lld=ld.copy(); numpy.fill_diagonal(lld,0)
    mmaxs=numpy.max(lld,axis=0)
    drop=[-1]
    if numpy.max(mmaxs)<tau:
        return ld, snps
    else:
        if max(mmaxs)>tau:
            drop=[whichMax(mmaxs)] # finds first occurrence of largest value>tau
        lsnps=list(snps); k=0;
        while (-1 in drop)==False:
            k=k+1
            lld=numpy.delete(lld,drop[0],axis=0)
            lld=numpy.delete(lld,drop[0],axis=1)
            # lld=subR(lld,drop);
            del lsnps[drop[0]]
            mmaxs=numpy.max(lld,axis=0)
            if max(mmaxs)>tau:
                drop=[whichMax(mmaxs)]
            else:
                drop=[-1]
        numpy.fill_diagonal(lld,1)
        return lld, lsnps
# my method is 5x faster than the lower triangular method

def fPos(x):
    cc=x<0
    x[cc]=0
    return x

def scad(z,a,alambda):
    cond1=(abs(z)<=(2*alambda))
    cond2=(((2*alambda)<abs(z)) & (abs(z)<=(a*alambda)))
    cond3=abs(z)>(a*alambda)
    z[cond1]=numpy.sign(z[cond1])*fPos(abs(z[cond1])-alambda)
    z[cond2]=((a-1)*z[cond2]-numpy.sign(z[cond2])*a*alambda)/(a-2)
    z[cond3]=z[cond3]
    return z

# def LDShrinkLower(ld,a=2,lambdaMax=0.1,nLambda=30,fun='scad'):
#     lambdas=list(numpy.linspace(0,lambdaMax,nLambda))
#     ddet=numpy.linalg.det(ld);ddets=[];k=0
#     nlld=nlld0=ld
#     while ddet>0 and lambdas[k]!=max(lambdas):
#         k=k+1; nlld0=nlld
#         if fun=='scad':
#             nlld=scad(ld,a,lambdas[k])
#         else:
#             nlld=ld; nlld[abs(nlld)<lambdas[k]]=0
#         ddet=numpy.linalg.det(nlld);ddets.append(ddet)
#     return nlld0, lambdas[k], k-1


# function to detect columns indices corresponding to different values
def indexFinder(toFind, allValues):
    gg=[allValues.index(toFind[x]) for x in range(0,len(toFind)) if toFind[x] in allValues]
    return gg

def columnFinder(columnNames):
    columnNames=[x.lower() for x in columnNames]
    # columnNames must be a list of column names in either the gene expression or outcome data set
    # goal is to auto-detect what columns in the data mean
    # columns to detect: SNP, CHR, BP, A1, A2, Z-statistics
    # potential names
    snpNames=['rsid','snp','markername'] # need to make sure not in CHR:BP:A1:A2 format
    bpNames=['pos','bp','bp_hg19','bphg19','pos_hg19','poshg19','snppos','snpbp']
    chrNames=['chr','chrom','chr_h19','chrhg19','snpchr','snpchrom']
    a1Names=['a1','allele1','effectallele','effect_allele','assessed_allele','assessedallele','ref']
    a2Names=['a2','allele2','noneffectallele','noneffect_allele','nonassessed_allele','assessedallele','other_allele','otherallele','ref']
    zNames=['z','zstat','zstatistic','z_stat','z_statistic','t','tstat','t_stat','t_statistic','teststatistic','test_statistic','zscore','z_score','tscore','t_score']
    # find their indices
    snpInd=indexFinder(snpNames, columnNames)
    bpInd=indexFinder(bpNames, columnNames)
    chrInd=indexFinder(chrNames, columnNames)
    a1Ind=indexFinder(a1Names, columnNames)
    a2Ind=indexFinder(a2Names, columnNames)
    zInd=indexFinder(zNames, columnNames)
    # double-check all were found
    keys=['SNP','BP','CHR','A1','A2','Z']; vals=[snpInd,bpInd,chrInd,a1Ind,a2Ind,zInd]
    if len(flatten_list(vals))<len(keys):
        raise ValueError('Could not find column name for ' + keys[vals.index([])] + ' column')
    # make dictionary
    dict={'snpInd': snpInd[0], 'bpInd': bpInd[0], 'chrInd': chrInd[0], 'a1Ind': a1Ind[0], 'a2Ind': a2Ind[0], 'zInd': zInd[0]}
    return dict

# bias-corrected ridge: no removing pleiotropy
def MRBEERidge(Bhat,Ahat,UU,UV,VV,Rinv,nLambda=100):
    # estimator is simple; lambdaMax determined automatically
    m=Bhat.shape[0];p=Bhat.shape[1];
    lambdaMax=numpy.max(abs(Bhat.T@Rinv@Ahat)/m)**(1/2); lls=numpy.linspace(0,lambdaMax,nLambda); logmse=[]
    for i in range(0,len(lls)):
        ee=numpy.linalg.inv(Bhat.T@Rinv@Bhat-m*UU+lls[i]*numpy.eye(p))@(Bhat.T@Rinv@Ahat-m*UV)
        res=Ahat-Bhat@ee
        toAdd=numpy.log(sum(res**2)[0])
        logmse.append(toAdd)
        if i>1:
            if logmse[-1]>logmse[-2]:
                break
    # final estimation
    lambda0=lls[i-1]; Ip=makeI(p)
    ee0=numpy.linalg.inv(Bhat.T@Rinv@Bhat-m*UU+lambda0*Ip)@(Bhat.T@Rinv@Ahat-m*UV)
    # variance-covariance
    Fhat=(1/m)*(Bhat.T@Rinv@Bhat+lambda0*Ip-m*UU)
    ResMat=numpy.zeros((m,m))
    numpy.fill_diagonal(ResMat, flatten_list(Ahat-Bhat@ee0))
    ff=numpy.zeros((m,p))
    dd=flatten_list([l.tolist() for l in UV-UU@ee0+lambda0*ee0])
    for i in range(0,ff.shape[0]):
        ff[i,:]=dd
    E=ResMat@Bhat-ff
    Vhat=E.T@Rinv@E/m
    Varhat=1/m*numpy.linalg.inv(Fhat)@Vhat@numpy.linalg.inv(Fhat)
    return ee0, Varhat, lambda0

def chiTestInternal(v,R):
    v=v.reshape((v.shape[0],1))
    cc_=v.T@numpy.linalg.inv(R)@v; cc_=cc_[0]
    cp_=1-stats.chi2.cdf(cc_,v.shape[0])
    return cp_

def IVChiTest(Zmat, R):
    Z0=numpy.nan_to_num(Zmat,nan=0)
    chiPs=[]
    for _ in range(0,Z0.shape[0]):
        ta=chiTestInternal(Z0[_,:], R)
        chiPs.append(ta)
    chiPs=numpy.array(chiPs)
    return chiPs

def Sop(theta,lambda_):
    toReturn=[]
    for jj in range(0,len(theta)):
        toS=numpy.abs(theta[jj])-lambda_
        z=toS if toS>0 else 0
        toReturn.append(numpy.sign(theta[jj])*z)
    return toReturn

def CD(x,y,UU,UV,VV,Rinv,lambda_, alpha_,eps=0.00001,max_iter=100,norm_stoppage=100,pen='scad'):
    p=x.shape[1]; m=x.shape[0]; I=makeI(p)
    # no adjustment for horizontal pleiotropy here
    est0=numpy.linalg.inv(x.T@Rinv@x+lambda_*I-m*UU)@(x.T@Rinv@y-m*UV)
    nns=[]; nn=100
    for cycle_ in range(0,max_iter):
        normPrev=nn
        bjhat=[]
        for p_ in range(0,p):
            xp_=x[:,p_]
            r_j=y-numpy.delete(x,obj=p_,axis=1)@numpy.delete(est0,obj=p_,axis=0)
            rx=(1/m)*(r_j.T@xp_); rx=rx[0] # turn to int
            dx=(1/m)*(xp_.T@xp_)+(1-alpha_)*lambda_
            if pen.lower()=='lasso':
                toAdd=Sop(rx,lambda_)/dx
            else:
                toAdd=scad(numpy.array(rx),3.7,lambda_)/dx 
            bjhat.append(toAdd)
        bjhat=numpy.array(bjhat); bjhat=bjhat.reshape((bjhat.shape[0],1))    
        nn=((est0-bjhat).T@(est0-bjhat))**2; nn=nn[0][0]; nns.append(nn)
        est0=bjhat;
        # if norm(est0-bjhat)>norm_stoppage and this norm>previous norm, just stop - not going to converge
        if nn<(eps*p) or ((nn>norm_stoppage and nn>normPrev)):
            break
    # end loops
    dd={'est': est0, 'k': cycle_}
    return dd

def mad(a):
    b=1.483*numpy.median(abs(a-numpy.median(a)))
    return b

def dSCAD(a,lam,gamma=3.7):
    a=abs(a); a=a.reshape((a.shape[0],))
    z=a.copy()
    mask1=a<lam
    mask2=a>lam
    mask3=a>(gamma*lam)
    z[mask1]=lam
    z[mask2]=(gamma*lam-z[mask2])/(gamma-1)
    z[mask3]=0
    return z

def soft(a,b):
    c=abs(a)-b
    mask=(c<0)
    c[mask]=0
    return c*numpy.sign(a)

### old MR Jones
# def MRJones(bX,by,Theta,UU,UV,lamvec,tauvec,rho_theta=1,rho_delta=1,max_iter=20,max_eps=0.01,alpha=0):
#     # Theta: solve(R) where R is LD matrix
#     p=bX.shape[1]; n=bX.shape[0]
#     BTB=bX.T@Theta@bX
#     BTa=bX.T@Theta@by
#     BT=bX.T@Theta
#     Ta=Theta@by
#     TB=Theta@bX
#     Trhoinv=numpy.linalg.inv(Theta+rho_delta*numpy.eye(n)) # may be able to make faster, but not that slow already
#     TITa=Trhoinv@Ta
#     TITB=Trhoinv@TB
#     TC=numpy.linalg.cholesky(Theta).T
#     # byse=TC@byse.copy()
#     # bXse=TC@bXse.copy()
#     Suu=n*UU.copy()
#     Suv=n*UV.copy()
#     theta_ini=(numpy.linalg.inv(BTB-Suu)@(BTa-Suv)).reshape((p,1))
#     mask=abs(theta_ini)<0.001
#     theta_ini[mask]=0
#     delta_ini=by-bX@theta_ini
#     mad_delta=mad(delta_ini)
#     mask=abs(delta_ini)<(mad_delta*3)
#     delta_ini[mask]=0
    
#     Btheta=numpy.zeros((p,len(lamvec),len(tauvec)))
#     Bdelta=numpy.zeros((n,len(lamvec),len(tauvec)))
#     Bbic=numpy.zeros((len(lamvec),len(tauvec)))
#     for ss in range(0,len(lamvec)):
#         for sss in range(0,len(tauvec)):
#             theta=theta_ini.copy()
#             theta1=theta.copy(); mask1=abs(theta1)<lamvec[ss]; theta1[mask1]=0
#             delta=delta_ini.copy()
#             delta1=delta.copy(); mask2=abs(delta1)<tauvec[sss]; delta1[mask2]=0
#             u1=rho_theta*(theta-theta1)
#             u2=rho_delta*(delta-delta1)
#             # normtheta=sum(theta**2)**0.5
#             normtheta=numpy.linalg.norm(theta,ord=2)
#             theta2=theta1.copy()*0
#             error=1; errors=[]
#             iter_=1
#             while (error>max_eps) & (iter_<max_iter):
#                 theta2=theta1.copy()
#                 adjratio=numpy.sum(delta1==0)/n
#                 Hinv=numpy.linalg.inv(BTB-Suu*adjratio+(rho_theta+alpha*lamvec[ss])*makeI(p))
#                 theta=Hinv@(BTa-BT@delta-Suv*adjratio-(u1-rho_theta*theta1))
#                 weight1=dSCAD(theta,lamvec[ss]).reshape((p,1))
#                 theta1=soft(theta+u1/rho_theta,weight1/rho_theta)
#                 delta=TITa-TITB@theta-Trhoinv@(u2-rho_delta*delta1)
#                 weight2=dSCAD(delta,tauvec[sss],gamma=2.2).reshape((n,1))
#                 delta1=soft(delta+u2/rho_delta,weight2/rho_delta)
#                 u1=u1.copy()+rho_theta*(theta-theta1)
#                 u2=u2.copy()+rho_delta*(delta-delta1)
#                 iter_=iter_+1
#                 if iter_>3:
#                     if normtheta<=0: # sometimes normtheta can be 0 (consider if all coefficients are shrunken to 0)
#                         error=max_eps+1
#                         errors.append(error)
#                     else:
#                         error=numpy.linalg.norm(theta2-theta1,ord=2)/normtheta
#                         # error=((theta2.T@theta1)**0.5)/normtheta
#                         errors.append(error)
#             mask1=(theta1!=0)
#             mask2=(delta1!=0)
#             theta1[mask1]=theta[mask1]
#             delta1[mask2]=delta[mask2]
#             Btheta[:,ss,sss]=theta1.reshape((p,))
#             Bdelta[:,ss,sss]=delta1.reshape((n,))
#             residual=by-bX@theta1-delta1
#             # indtheta=(theta1!=0); indtheta=indtheta.reshape((indtheta.shape[0],))
#             # LinSm=BTB[indtheta,:][:,indtheta]-Suu[indtheta,:][:,indtheta]*adjratio+alpha*lamvec[ss]*makeI(sum(indtheta));
#             # LinSm=numpy.linalg.inv(LinSm)@(BTB[indtheta,:][:,indtheta]-Suu[indtheta,:][:,indtheta]*adjratio)
#             # df1=sum(numpy.diag(LinSm))
#             df1=numpy.sum(theta1!=0) # generally gives the same results as using df1 above, especially since lambda is usually taken to be extremely small
#             df2=numpy.sum(delta1!=0)
#             rss=TC@residual; rss=rss.reshape((rss.shape[0],1))
#             if (n-df1-df2)<1:
#                 Bbic[ss,sss]=math.inf
#             else:
#                 sig=(rss.T@rss)/(n-df1-df2)
#                 Bbic[ss,sss]=n*numpy.log(sig)+numpy.log(n)*(df1+df2)
#     A={}
#     A['thetamatrix']=Btheta
#     A['deltamatrix']=Bdelta
#     A['bic']=Bbic
#     return A

### new MR Jones
def MRJones(by,bX,LD,Rxx,rxy,lamvec=numpy.linspace(0.05,0.3,10),tauvec=numpy.linspace(5,25,10),rho_theta=1,rho_gamma=1,max_iter=30,max_eps=0.005,intercept=True,WIB=True,WIB_thres=5,bicfactor=1,normmax=3):
    warnings.filterwarnings('ignore')
    n=len(by.squeeze());
    if intercept:
        bX=numpy.column_stack(([1]*n,bX))
        Rxx=numpy.column_stack(([0]*Rxx.shape[0],Rxx)); Rxx=numpy.row_stack(([0]*(Rxx.shape[1]),Rxx))
        rxy=numpy.row_stack(([0],rxy))
    p=bX.shape[1]
    q=len(lamvec)
    w=len(tauvec)
    Theta=numpy.linalg.inv(LD)
    TC=numpy.linalg.cholesky(LD).T
    BT=bX.T@Theta
    BTB=BT@bX
    BTa=BT@by
    adjvec=numpy.diag(BTB-n*Rxx); adjvec.flags.writeable=True
    mask=adjvec<0
    adjvec[mask]=1/(n**3)
    adjvec=adjvec**0.5
    adjvec=adjvec/numpy.max(adjvec)
    
    if WIB:
        ss=numpy.diag(BTB)/n
        mask=ss<WIB_thres
        qq=numpy.diag(Rxx); qq.flags.writeable=True
        qq[mask]=(qq[mask]*ss[mask]/WIB_thres)**0.5
        Rxx=numpy.diag(qq)@Rxx@numpy.diag(qq)
    
    Thetarho=numpy.linalg.inv(LD+rho_gamma*numpy.diag([1]*n))
    by1=(TC@by).squeeze()
    if intercept:
        bX1=TC@bX[:,1:]
        fit0=QuantReg(by,bX1) 
        res=fit0.fit(q=0.5) # can give a warning if max iterations reached... how to silence?
        theta_ini=res.params; theta_ini=theta_ini.reshape((theta_ini.shape[0],1)); theta_ini.flags.writeable=True
        mask=numpy.abs(theta_ini)<0.001
        theta_ini[mask]=0
        theta_ini=numpy.row_stack((0,theta_ini))
        adjvec[0]=1
    else:
        bX1=TC@bX
        fit0=QuantReg(by,bX1)
        res=fit0.fit(q=0.5)
        theta_ini=res.params; theta_ini=theta_ini.reshape((theta_ini.shape[0],1)); theta_ini.flags.writeable=True
        mask=numpy.abs(theta_ini)<0.001
        theta_ini[mask]=0
        theta_ini=numpy.row_stack((0,theta_ini))
    
    normmax=max(numpy.sqrt(sum(theta_ini**2))*normmax); normmax=normmax if normmax>1 else 1
    r=(by-bX@theta_ini).squeeze()
    gamma_ini=(Theta@r).squeeze()
    mask=abs(gamma_ini/3)<mad(gamma_ini)
    gamma_ini[mask]=0
    
    lamvec=numpy.sort(lamvec)
    tauvec=numpy.sort(tauvec)
    Btheta=numpy.zeros((p,q,w))
    Bgamma=numpy.zeros((n,q,w))
    Bbic=numpy.zeros((q,w))
    
    for ss in range(0,q):
        for sss in reversed(range(0,w)):
            theta=theta_ini.copy().squeeze()
            theta1=theta.copy().squeeze()
            gamma=gamma_ini.copy().squeeze()
            gamma1=gamma.copy().squeeze()
            u1=rho_theta*(theta-theta1); u1=u1.squeeze()
            u2=rho_gamma*(gamma-gamma1); u2=u2.squeeze()
            theta2=theta1*0
            error=1; iter_=0
            while (error>max_eps) & (iter_<max_iter):
                iter_=iter_+1
                theta2=theta1.copy()
                effn=sum(gamma1==0)
                Hinv=numpy.linalg.inv(BTB-effn*Rxx+rho_theta*numpy.diag([1]*p))
                theta=Hinv@(BTa.squeeze()-bX.T@gamma-(rxy.squeeze())*effn-(u1.squeeze()-rho_theta*theta1.squeeze()))
                if numpy.sqrt(sum(theta.squeeze()**2))>normmax:
                    theta=theta/numpy.sqrt(sum(theta**2))*normmax
                
                theta1=mcp(theta+u1.squeeze()/rho_theta,lamvec[ss]/adjvec/rho_theta,ga=3)
                gamma=Thetarho@(by.squeeze()-bX@theta-u2.squeeze()+rho_gamma*gamma1)
                gamma1=mcp(gamma+u2.squeeze()/rho_gamma,tauvec[sss]/rho_gamma,ga=3)
                u1=u1+rho_theta*(theta-theta1)
                u2=u2+rho_gamma*(gamma-gamma1)
                error=sum((theta2.squeeze()-theta1.squeeze())**2) if iter_>3 else 1
                # end while
            
            Btheta[:,ss,sss]=theta1.squeeze()
            Bgamma[:,ss,sss]=gamma1.squeeze()
            df1=sum(theta1!=0)
            df2=sum(gamma1!=0)
            r=(by.squeeze()-bX@theta1-LD@gamma1)
            varr=mad(TC@r)**2
            if varr==0: # if MRJones wanted all to be pleiotropy
                rss=1e5 # something very large bc I never want this solution
            else:
                rss=sum(r*(Theta@r)/varr)
            Bbic[ss,sss]=n*numpy.log(rss)+(numpy.log(n)+numpy.log(p)*bicfactor)*df1+(numpy.log(n)+numpy.log(n)*bicfactor)*df2
            # end for sss
    # end for ss
    istar,jstar=numpy.where(Bbic==Bbic.min()); istar=istar[0]; jstar=jstar[0] # may be multiple minimums? Choose first
    theta=Btheta[:,istar,jstar]
    gamma=Bgamma[:,istar,jstar]
    # ignore everything intercept-related from the beginning
    if intercept:
        theta[0]=0 # instead of deleting bc I use bX[:,masktheta] below
        covg=numpy.zeros((p-1,p-1))
    else:
        covg=numpy.zeros((p,p))
    indtheta=[i for i in range(0,len(theta)) if theta[i]!=0]
    indgamma=[i for i in range(0,len(gamma)) if gamma[i]!=0]
    effn=n-len(indgamma)
    masktheta=numpy.zeros((len(theta)),dtype='bool'); masktheta[indtheta]=True
    maskgamma=numpy.ones((len(r)),dtype='bool'); maskgamma[indgamma]=False
    
    if (sum(indtheta)>0) & (sum(maskgamma)>2): # must have at least 2 non pleiotropic IVs
        Z=bX[:,masktheta][maskgamma,:]
        R0inv=numpy.linalg.inv(LD[maskgamma,:][:,maskgamma])
        r=by.squeeze()-bX@theta-LD@gamma.squeeze()
        r=r[maskgamma]
        D=numpy.diag(r.squeeze())
        V=Z.T@R0inv@D@(LD[maskgamma,:][:,maskgamma])@D@R0inv@Z
        Hinv=numpy.linalg.inv(Z.T@R0inv@Z-sum(maskgamma)*Rxx[masktheta,:][:,masktheta])
        covg=Hinv@V@Hinv
    
    # if intercept:
    #     intercept=theta[0]
    #     theta=theta[1:]
    #     # covg=covg[1:,:][:,1:] # don't need this bc I already made theta for intercept 0 above
    #     bX=bX[:,1:]
    if intercept:
        theta=theta[1:]
    
    dout={}
    dout['theta']=theta
    dout['gamma']=gamma
    dout['selectedgenes']=indtheta
    dout['pleiotropyIVs']=indgamma
    dout['covg']=covg
    dout['Bic']=Bbic
    dout['Btheta']=Btheta
    dout['Bgamma']=Bgamma
    dout['intercept']=intercept
    warnings.filterwarnings('default')
    return dout

def mcp(a,lam,ga=3):
    b=abs(a)
    z=soft(a,lam)/(1-1/ga)
    mask=b>(ga*lam)
    z[mask]=a[mask]
    return z

def bimin(mat):
    return numpy.where(mat==mat.min())

# search the BIC grid outputted by MR Jones and find optimal lambda1 and lambda2 (tau)
# def BICgridSearch(MRJonesOut,lamvec,tauvec):
#     # lambda2/tau (for pleiotropy) is in columns; lambda1 (for causal effects) is in rows
#     lamvec=lamvec.reshape((lamvec.shape[0],))
#     tauvec=tauvec.reshape((tauvec.shape[0],))
#     BICgrid=MRJonesOut['bic']
#     EstsGrid=MRJonesOut['thetamatrix']
#     PleioGrid=MRJonesOut['deltamatrix']
#     BICgrid=numpy.nan_to_num(BICgrid,nan=math.inf) # make large bc I'm searching for the minimum
#     inds=numpy.where(BICgrid==BICgrid.min())
#     i,j=inds[0][0], inds[1][0]
#     ests=EstsGrid[:,i,j]
#     deltas=PleioGrid[:,i,j]
#     isPleio=(deltas!=0)
#     A={'finalEsts': ests, 'finalDeltas': deltas, 'isPleio': isPleio, 'minLambda1': lamvec[i], 'minLambda2': tauvec[j], 'minBIC': BICgrid[i,j], 'i': i, 'j': j}
#     return A

# def MRJonesVariance(bX,by,deltas,ldR,UU,UV,finalEsts,isPleio,finalLambda,finalTau,alpha):
#     # Theta: inverse of LD matrix
#     # UU: SigmaUU; UV: SigmaUV
#     # finalEsts: optimized MRJones causal estimates
#     # isPleio: n-length vector of True/False values indicating if xi_j/gamma_j != 0 (ie is pleio after shrinkage)
#     # finalLambda: optimal lambda (lambda1) returned by MRJones; finalTau: optimal lambda 2 (for pleio) returned by MRJones
#     pleiomask=isPleio==False
#     thetamask=finalEsts!=0
#     X1=bX[pleiomask,:][:,thetamask]
#     y=by[pleiomask,:]
#     ld_=ldR[pleiomask,:][:,pleiomask]
#     ld_=numpy.linalg.inv(ld_)
#     UU_=UU[thetamask,:][:,thetamask]
#     UV_=UV[thetamask,:]; UV_=UV_.reshape((UV_.shape[0],))
#     m0=sum(pleiomask); p0=sum(thetamask)
#     H=X1.T@ld_@X1-m0*UU_+alpha*finalLambda*makeI(p0)
#     H=numpy.linalg.inv(H)
#     est0=finalEsts[thetamask].reshape((finalEsts[thetamask].shape[0],))
#     res=y.reshape((y.shape[0],))-X1@est0
#     ResMat=numpy.zeros((m0,m0))
#     numpy.fill_diagonal(ResMat, flatten_list(res))
#     ff=numpy.zeros((m0,p0))
#     dd=flatten_list([l.tolist() for l in UV_-UU_@est0])
#     for i in range(0,ff.shape[0]):
#         ff[i,:]=dd
#     E=ResMat@X1-ff
#     V=E.T@ld_@E
#     VB=H@V@H
#     return VB

def CD_BIC(x,y,UU,UV,VV,Rinv,lambda_min,lambda_max,alpha_,nLambda=100,eps=0.00001,max_iter=100,norm_stoppage=100,pen='scad'):
    # x: bx0, y: by, UU: SigmaUU, UV: SigmaUV, VV:SigmaVV, lambda_min/_max: minimum/maximum lambda, alpha_: alpha
    lambdas=numpy.linspace(lambda_min,lambda_max,nLambda)
    bic=[]
    for ll in range(0,lambdas.shape[0]):
        est_ll=CD(x,y,UU,UV,VV,Rinv,lambdas[ll],alpha_,eps,max_iter,norm_stoppage,pen)['est']
        mse=sum((y-x@est_ll)**2)[0]
        # df is trace of hat matrix, but that can be computationally costly since it is one large matrix inversion per lambda ... stick with #nonzero
        # Hat=x@numpy.linalg.inv(x.T@Rinv@x-m*UU)@x.T
        # HatRidge=x@numpy.linalg.inv(x.T@Rinv@x-m*UU+lambdas[ll]*makeI(p))@x.T
        # dfRidge=numpy.trace(HatRidge)
        # df=numpy.trace(Hat) # not even clear what df are for MRBEE
        effDF=len([x for x in range(0,est_ll.shape[0]) if est_ll[x][0]!=0])
        bic.append(numpy.log(mse)/m+effDF*numpy.log(m)/m)
    dd={'BICs': bic, 'lambda_min': lambdas[bic.index(min(bic))]}
    return dd

# MaxDet function for filling in correlation matrix used by MRBEE (requires the function from above: flatten_list())
def MaxDet(R,doToBadSolutions='drop'):
    if (doToBadSolutions=='drop' or doToBadSolutions=='zero')==False:
        raise ValueError('argument doToBadSolutions should be one of drop or zero')
    runningR=R.copy(); ns=[]; bad_variates=[]
    # for each row
    for i_ in range(0,runningR.shape[0]):
        # for each column
        for j_ in range(0,runningR.shape[0]):
            if i_>j_ or i_==j_: # working on upper-triangle only (all will get filled in)
                continue
            else:
                if math.isnan(R[i_,j_])==False:
                    continue
                else:
                    tempR=R.copy()
                    # reorder, putting target at (n-1)xn and nx(n-1) positions
                    inds_=list(numpy.array(range(0,tempR.shape[0])))
                    # (i_, j_) has a missing value
                    # steps
                    # 1) put j_ column last and fix rows
                    # find index location of j_ in inds_ and put j_ at the end
                    a=inds_[:(inds_.index(j_))]; b=inds_[(inds_.index(j_)+1):]
                    abc=[a,b,j_]; abc=flatten_list(abc)
                    tempR=tempR[:,abc][abc,:] # does not change the row that missing is on
                    # 2) put i_ row last and fix columns
                    a=inds_[:(inds_.index(i_))]; b=inds_[(inds_.index(i_)+1):]
                    abc=[a,b,i_]; abc=flatten_list(abc)
                    tempR=tempR[abc,:][:,abc]
                    # now missing value is in (n-1)xn and nx(n-1) spot
                    # remove other missing, after making target non-missing
                    n=tempR.shape[0]
                    tempR[n-1,n-2]=99; tempR[n-2,n-1]=99
                    # remove any rows/columns with missing
                    # first remove variables where the correlation with the target is missing
                    lastColMiss=numpy.isnan(tempR[:,-1])==False
                    tempR=tempR[:,lastColMiss][lastColMiss,:]
                    # now remove other missing
                    # the row with the target should no longer have any missing
                    # make upper triangle 0
                    tempRLower=numpy.tril(tempR)
                    mask=numpy.isnan(tempRLower).any(axis=1)==False
                    tempR=tempR[mask,:][:,mask]
                    # impute single missing with 0 and make positive definite (I hope this will ensure solutions stay in {-1,1}
                    # mask=numpy.isnan(tempR)
                    # tempR[mask]=0
                    tempR,alpha,deco=posDefifyCorrMat(tempR,0.001) # keep positive definite to avoid bad solutions
                    # identify sub-matrices used for MaxDet solution
                    n=tempR.shape[0]
                    A11=tempR[:(n-2),0:(n-2)]
                    B=tempR[(A11.shape[0]),0:A11.shape[1]]
                    C=tempR[tempR.shape[0]-1,0:A11.shape[0]]
                    sol=B.T@numpy.linalg.inv(A11)@C
                    ns.append(tempR.shape[0]-1)
                    if abs(sol)>1:
                        bad_variates.append((i_,j_))
                        # print(str(i_)+'and'+str(j_)+'sol='+str(sol))
                        if doToBadSolutions=='zero':
                            sol=0
                    runningR[i_,j_]=sol
                    runningR[j_,i_]=sol
    # It seems like one variable can cause multiple problems and the variable will always
    # be the first value in the tuples which are listed in bad_variates.
    # find the variables whose indices are listed first in each of the tuples in bad_variates
    bad_inds=set([bad_variates[x][0] for x in range(0,len(bad_variates))])
    bad_inds=list(bad_inds)
    goodMask=numpy.ones((runningR.shape[0],),dtype='bool')
    goodMask[bad_inds]=False
    nDropped=0 if doToBadSolutions=='zero' else len(bad_inds)
    if doToBadSolutions=='drop':
        runningR=runningR[goodMask,:][:,goodMask]
    return runningR, bad_inds, nDropped

def nuclearNorm(A):
    return sum(numpy.linalg.svd(A)[1])

def fullRankMatrixCompletion(A,lambda_,missImp=0,thresholding='soft',scad_a=3.7,eps=0.01,max_iter=100):
    nn=A.shape[0]*A.shape[1]; nns=[]; k=0
    Ac=A.copy(); Omega=numpy.isnan(Ac); 
    Ac[Omega]=missImp # we hypothesize the mean of missing to be 0 (these are just starting values; the program will add some error)
    lam=lambda_ # may affect results
    while (nn>eps) & (k<max_iter):
        k=k+1
        deco=numpy.linalg.svd(Ac,full_matrices=False)
        newd=Sop(deco[1],lam) # soft-thresholing
        # newd=scad(deco[1],scad_a,lam) # SCAD 
        lam=numpy.min(newd)/2 # kind of arbitrary, but don't want to set the smallest eigenvalue to 0
        AcPred=deco[0]@numpy.diag(newd)@deco[2] # no need to transpose the last
        AcPred[Omega==False]=Ac[Omega==False] # don't change observed values
        nn=nuclearNorm(Ac-AcPred); nns.append(nn)
        Ac=AcPred.copy()
    return Ac,k,nns

def taperLD(R,k):
    firstCol=numpy.linspace(1,R.shape[0],R.shape[0]); firstCol[:k]=1
    To=1/scipy.linalg.toeplitz(firstCol)
    Rtapered=R*To # will retain 1's on diagonal
    return Rtapered

def bandLD(R,k):
    if k>R.shape[0]:
        return R
    else:
        firstCol=numpy.zeros((R.shape[0],)); firstCol[:k]=1
        To=scipy.linalg.toeplitz(firstCol)
        Rbanded=R*To
        return Rbanded

def softThresholdLD(R,lambda_):
    Rsoft=soft(R,lambda_) # will NOT retain 1s on diagonal
    numpy.fill_diagonal(Rsoft,1) # retain correlation matrix structure
    return Rsoft

def regularizeLD(R,kTaper,kBand,lambda_,epsilon=1e-4,toCorr=True,justMCP=True):
    # epsilon: minimum eigenvalue final matrix is allowed to have
    # toCorr: If True, converts solution to a correlation matrix, else diagonal elements may be outside of {-1,1} 
    if justMCP:
        r3=mcp(R,lambda_)
        r4a,alpha,deco=posDefifyCorrMat(r3,epsilon)
    else:
        r1=taperLD(R,kTaper)
        r2=bandLD(r1,kBand)
        r3=softThresholdLD(r2,lambda_)
        r4a,alpha,deco=posDefifyCorrMat(r3,epsilon)
    if toCorr:
        D=1/numpy.diag(r4a)**0.5
        D=numpy.diag(D)
        r4a=D@r4a@D
    return r4a, alpha, deco

def posDefifyCorrMat(ndmatrix,epsilon=1e-4):
    # ndmatrix: pxp square negative-definite matrix of measurement error correlations/covariances
    # since it is a correlation matrix, abs max off-diagonal value is 1
    deco=numpy.linalg.eig(ndmatrix)
    eig0=min(numpy.real(deco[0])); eigP=max(numpy.real(deco[0]))
    if eig0<epsilon:
        mu=max([epsilon,(eig0+eigP)/2])
        alpha=(mu-epsilon)/(mu-eig0) # optimal choice given by Choi et al: https://doi.org/10.1016/j.jmva.2018.12.002
        Sigmahat=alpha*ndmatrix+(1-alpha)*mu*numpy.eye(ndmatrix.shape[0]) # the solution
        deco=numpy.linalg.eig(Sigmahat)
    else:
        alpha=1
        Sigmahat=ndmatrix.copy()
    deco=(numpy.real(deco[0]),numpy.real(deco[1]))
    return Sigmahat,alpha,deco

# function to calculate conditional F-statistics
# !!! at the moment, this function may be extremely memory-intensive (consider the kronecker products)
def conditionalFMRBEE(bx,UU,sparsePrecisionSqrt,tauLowerQuant=0,tauUpperQuant=0,opAlphaVal=0,method='mrbee'):
    # MRJones and MRBEE with IMRP adjustment give slightly different estimates of the conditional Fs
    # I trust MRJones more because it also controls for high correlations between exposures, which MRBEE cannot explicitly do
    bx_=sparsePrecisionSqrt@bx.copy(); Theta=makeI(bx.shape[0])
    Fs=[]
    for _ in range(0,bx.shape[1]):
        x_=numpy.delete(bx_,(_),axis=1); x_=x_.reshape((x_.shape[0],bx.shape[1]-1))
        y_=bx_[:,_]; y_=y_.reshape((y_.shape[0],1))
        UU_=numpy.delete(UU,(_),axis=0); UU_=numpy.delete(UU_,(_),axis=1)
        UV_=UU[:,_]; UV_=numpy.delete(UV_,(_),axis=0); UV_=UV_.reshape((UV_.shape[0],1))
        # if user just wants to use MRBEE (basically) to find conditional F-statistics (much faster)
        if method.lower()=='mrbee':
            ### using MRBEE-IMRP
            fit=imrbee(x_,y_,UU_,UV_,numpy.array([1]).reshape((1,1)),Theta,numpy.linalg.inv(Theta),PleioPThreshold=0.05/(x_.shape[0]**0.5),boot=False)
            finalEsts=fit[0] # est0, V0, Outliers, kR, k
            finalEsts=finalEsts[1:] # ignore intercept
            mask=numpy.ones((finalEsts.shape[0],),dtype='bool');
            estsVars=fit[1][1:,:][:,1:] # remove items corresponding to intercept (automatically added by MRBEE)
            df1=sum(mask);df2=0
        else:
            lamMax=(numpy.max(abs(x_.T@y_))/bx.shape[0])**0.5
            lamvec=numpy.linspace(0,lamMax,20)
            e0=numpy.linalg.inv(x_.T@x_-x_.shape[0]*UU_)@(x_.T@y_-x_.shape[0]*UV_) # initial unbiased estimate with no shrinkage or pleiotropy adjustment
            res=y_.reshape((y_.shape[0],))-(x_@e0).reshape((y_.shape[0],));
            kt1=max(abs(res))*tauLowerQuant; kt2=max(abs(res))*tauUpperQuant
            tauvec=numpy.linspace(kt1,kt2,20)
            ### using MRJones
            fit=MRJones(x_,y_,Theta,UU_,UV_,lamvec,tauvec,rho_theta=1,rho_delta=1,max_iter=10,max_eps=0.01)
            outSearched=BICgridSearch(fit,lamvec,tauvec) # finding optimal hyperparameters
            finalEsts=fit['thetamatrix'][:,outSearched['i'],outSearched['j']] # extracting final causal estimates from output
            mask=finalEsts!=0; 
            finalDeltas=fit['deltamatrix'] [:,outSearched['i'],outSearched['j']] # extracting final xi/horizontal pleiotropy estiamtes from output
            finalLambda=outSearched['minLambda1']; finalTau=outSearched['minLambda2'] # definining optimal lambda1 and lambda2/tau
            estsVars=MRJonesVariance(x_,y_,finalDeltas,Theta,UU_,UV_,finalEsts,finalDeltas!=0,finalLambda,finalTau,opAlphaVal) # estimating variance of causal estimates
            df1=sum(mask); df2=sum(finalDeltas!=0)
        
        # chi-square statistic (since I pre-multiplied by sparsePrecisionSqrt, all IVs are independent)
        finalEsts=finalEsts.reshape((finalEsts.shape[0],1))
        var_=1+float(finalEsts.T@finalEsts)+numpy.diag(x_[:,mask]@estsVars@x_[:,mask].T)
        res=y_.squeeze()-(x_[:,mask]@finalEsts[mask]).squeeze()
        chi=sum(res**2/var_)
        condF=chi/(x_.shape[0]-df1-df2-1)
        Fs.append(condF)
    return Fs

# organize causal estimate results
def organizeMetaResults(_dict_):
    dNames=list(_dict_)
    # in a running data frame with _dict_ elements indexed by group ID (ie gene group number)
    df0=pandas.DataFrame.from_dict(_dict_[dNames[0]]['meta'])
    df0['groupID']=dNames[0]
    for _ in range(1, len(dNames)):
        innerD=_dict_[dNames[_]]['meta']
        pp=pandas.DataFrame.from_dict(innerD)
        pp['groupID']=dNames[_]
        df0=pandas.concat([df0,pp])
    #df0=df0.reset_index()
    #df0=df0.drop(columns='index', axis=1)
    return df0

def concatUniMRRes(outerDict):
    keys=list(outerDict)
    start=outerDict[keys[0]]
    for II in range(1,len(keys)):
        start=pandas.concat((start,outerDict[keys[II]]))
    return start

# organize things to monitor
def organizeThingsMonitored(_dict_):
    df0=pandas.DataFrame.from_dict(_dict_).T
    df0['groupID']=df0.index
    return df0

# function to separate full eQTLGen GWAS data (all CHRs) into chromosome-specific files (speeds later operations up)
# def splitEQTLGenDataIntoCHRs(fulleqtlgenfp,diryouwant,verbose=True):
#     # fulleqtlgenfp: filepath to full eQTLGen data (with file extension)
#     # diryouwant: the directory in which you want me to make a new directory called chrSpecificEQTLs that I'll put the chr-specific data in
#     # load full data
#     if verbose:
#         print('loading full eQTLGen data for all chromosomes. please note that due to the size of the full eQTLGen data, this process may take a while (~20 minutes)')
#     fullData=pandas.read_table(fulleqtlgenfp,sep=None,engine='python')
#     allchrs=fullData['SNPChr'].unique()
#     for _ in range(0,len(allchrs)):
#         if verbose:
#             print('starting chromosome '+thischr)
#         datachr=fullData[fullData['SNPChr']==allchrs[_]]
#         fpOut=diryouwant+'/chr'+thischr+'eqtls.txt'
#         numpy.savetxt(fpOut, datachr)
#         cmd=["gzip", fpOut]
#         out=subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
#         if verbose:
#             print('finished chromosome '+thischr)

# Yihe's matrix completion
def poet(A,k):
    u,s,vh=numpy.linalg.svd(A)
    u=u[:,:k]
    return u

def ffimpute(X,K,dg,P,W,iter_max=25,eps=0.001):
    n,p=X.shape
    Del=numpy.isnan(X)
    # initialize
    X1=X.copy(); X1[Del]=0
    S=P@X1.T@W@X1@P
    Gamma=poet(S,K)
    Factor=X1@P@Gamma
    L=Factor@Gamma.T
    X1[Del]=L[Del]
    Gamma1=Gamma*0
    iter_=0
    while (numpy.linalg.norm(Gamma-Gamma1,ord='fro')>(eps*p**0.5)) & (iter_<iter_max):
        Gamma1=Gamma.copy()
        S=P@X1.T@W@X1@P
        Gamma=poet(S,K)
        Factor=X1@P@Gamma
        L=Factor@Gamma.T
        X1[Del]=L[Del]
        iter_=iter_+1
    return L, Factor, Gamma

def ffimputeBIC(X,Kmax,dg,dim_proj,W,iter_max=15,eps=0.001,tuning=0.001):
    if dim_proj>0:
        P=projectm(dg,k=dim_proj,tuning=tuning)
    if dim_proj==0:
        P=makeI(X.shape[1])
    n,p=X.shape
    Del=numpy.isnan(X);
    Del=Del+1-1 # converts to numeric (1=missing, 0=not)
    notDel=1-Del.copy()
    neff=numpy.sum(notDel)
    X1=X.copy(); X1[Del==1]=0
    bic=[]
    for kk in range(1,Kmax+1):
        L,Factor,Gamma=ffimpute(X=X,K=kk,P=P,W=W,dg=dg,iter_max=iter_max,eps=eps)
        df=n*kk+kk*p
        toAdd=numpy.log(numpy.linalg.norm(X1@P-notDel*L,ord="fro")**2)+numpy.log(neff/(n+p))*df/neff
        bic.append(toAdd)
    kopt=bic.index(min(bic))
    L,Factor,Gamma=ffimpute(X=X,K=kopt,P=P,W=W,dg=dg,iter_max=iter_max,eps=eps)
    return L, Factor, Gamma

def projectm(dg,bdeg=3,k=10,tuning=0):
    q1=numpy.quantile(dg,0.05)
    x=dg-numpy.min(dg)
    x=x/numpy.max(x)
    B,D=cosinebase(x,k)
    P=B@numpy.linalg.inv(B.T@B+tuning*D/p)@B.T
    return P

def cosinebase(t,k):
    p=len(t)
    A=numpy.ones((p,k+1))
    d=[0]
    for _ in range(1,k+1):
        # A=numpy.column_stack((A,numpy.sqrt(2)*numpy.cos(numpy.pi*_*t)))
        A[:,_]=numpy.sqrt(2)*numpy.cos(numpy.pi*_*t)
        d.append(_**4)
    D=numpy.zeros((k+1,k+1)); numpy.fill_diagonal(D,d)
    return A, D

# def specialToep(_): # _: number of rows/columns
#     hh=numpy.zeros((_,_))
#     for s_ in range(0,hh.shape[0]-1):
#         hh[s_:s_+2,s_:s_+2]=numpy.array([[2,-1],[-1,2]])
#     return hh

# def bspline(x, xl, xr, ndx, bdeg):
#     dx=(xr-xl)/ndx
#     knots=list(range(int(xl-bdeg*dx), int(xr+bdeg*dx), int(dx)))    
#     B=statsmodels.gam.smooth_basis.BSplines(x,df=[len(knots)],degree=[bdeg+1])
#     D=specialToep(B.shape[1]) # Yihe used in R: D=toeplitz(c(2,-1,rep(0,ncol(B)-2)))
#     D=D.T@D
#     return B, D

# writes out LD matrix to writableDir to be read in by another program (R program I will make) to make plots 
def LDLeadGeneSNPs(res, writableDir=os.getcwd(),ldRefDir='/mnt/rstor/SOM_EPBI_XXZ10/njl96/data/1000G/1kg_phase3_EUR_only'):
    # res: direct output from organizeMetaResults(outerDict)
    nondupres=res.drop_duplicates(subset='LeadSNP',keep='first')
    # note, some genes may share the same lead SNPs ... PLINK will only consider a SNP once. 
    outDir1=writableDir+'/myExtract.txt' # (for saving SNPs for LD calculations) write out to a directory that must be writable (since the user is there)
    outDir2=writableDir+'/tempOut' # (for saving LD matrix)
    res['LeadSNP'].to_csv(outDir1, header=False, index=False, sep=" ") # writing out list of SNPs to give PLINK to make LD calculations
    cmd=[callPlink(), "--r", "square", "--bfile", ldRefDir, "--extract", outDir1, "--out", outDir2] # the command for PLINK
    out=subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) # executing the command in the terminal
    ldmat=numpy.loadtxt(writableDir+'/tempOut.ld', dtype='f') # reading LD matrix back in; matrix format but with duplicates occupying only one row/column
    ### maybe for the plotting this does not matter actually ... 
    # ldmat is for unique SNPs only. write out correlation matrix for nondupres, which is all SNPs but none repeated >once, and the group start indices
    # if writing out for another program to read in
    os.chdir(writableDir)
    numpy.savetxt('ldmat.txt',ldmat)
    res.to_csv('res.txt',index=False)
    nondupres.to_csv('nondupres.txt',index=False)
    # note what this function does in the end - it doesnt return anything; it writes out a file

def ldscore(ld0):
    ell=[]; mask=numpy.zeros((ld0.shape[0],),dtype='bool')
    for _o in range(0,ld0.shape[0]):
        mask_o=mask.copy(); mask_o[_o]=True
        _v=ld0[mask_o,:][:,(mask_o==False)].reshape((ld0.shape[0]-1,))
        ell.append(sum(_v**2))
    ell=numpy.array(ell)
    return ell

def wmean(x,w):
    return sum(x*w)/sum(w)

def TwoMixtureEM(data,mu01,mu02,sig01,sig02,eps=0.001,k=100):
    p01=norm.pdf(data,mu01,sig01)
    p02=norm.pdf(data,mu02,sig02)
    error=1; k_=0; es=[]
    gg=numpy.ones((data.shape[0]),dtype='int')
    while (error>eps) & (k_<k):
        k_=k_+1
        mask=p01<p02 # True if predicted group 2, False if predicted group 1
        mu1est=wmean(data[mask==False],p01[mask==False]) # can change weights if you want
        mu2est=wmean(data[mask],p02[mask]) 
        sig1est=numpy.std(x[mask==False])
        sig2est=numpy.std(x[mask])
        p01=norm.pdf(data,mu1est,sig1est)
        p02=norm.pdf(data,mu2est,sig2est)
        error=(mu1est-mu01+mu2est-mu02+sig1est-sig01+sig2est-sig02)**2
        mu01=mu1est; mu02=mu2est; sig10=sig1est; sig20=sig2est # should all be floats so no need to use copies
        es.append(error)
        if k_>1:
            if abs(es[-2]-es[-1])<0.0001: # another stopping rule - if error did not change much at this iteration
                error=0
    gg[mask]=2; gg=gg-1
    return gg

# # testing the mixture
# n=1000;eps=0.7;mu1=0;mu2=0;sig1=(30000*0.2/100)**0.5;sig2=(1/30000)**0.5 # eps is P(group 1)
# mu1=0;mu2=1;sig1=1;sig2=1
# bj=numpy.random.binomial(1,eps,n)
# x=bj*numpy.random.normal(mu1,sig1,n)+(1-bj)*numpy.random.normal(mu2,sig2,n)
# mu01=0;mu02=-1;sig01=1;sig02=1
# bjhat=TwoMixtureEM(x,mu01,mu02,sig01,sig02,eps=0.01,k=100)
# import pandas
# df=pandas.DataFrame({'x': bj, 'xhat': bjhat})
# tt=numpy.array(pandas.crosstab(df['x'],df['xhat']))
# sum(numpy.diag(tt))/numpy.sum(tt)

# def bootOutcomeH2(x,H0target=x.shape[0],k=1000):
#     ests=[]
#     for k_ in range(0,k):
#         xstar=numpy.random.choice(x,size=x.shape[0],replace=True)
#         ests.append(sum(xstar**2))
#     ests=numpy.array(ests)
#     p_=1-sum(ests>H0target)/k
#     qs=numpy.quantile(ests,q=[0.05,0.95]) # 10% of data is outside of the range {0.05,0.95}
#     h2est=by.T@by-by.shape[0]
#     return h2est, p_, qs

def trace(A):
    out=numpy.sum(numpy.diag(A))
    return out

# single SNP MR
def singleSNP(bx,by,cn): # assumes Z-stat standardization
    m=bx.shape[0];p=bx.shape[1];by=by.squeeze()
    wald=[];waldse=[]
    for _ in range(0,p):
        bxp=bx[:,_].squeeze()
        maxix=numpy.argmax(abs(bxp))
        wald.append(by[maxix]/bxp[maxix])
        waldse.append(1/bxp[maxix]**2)
    wald=numpy.array(wald)
    waldse=numpy.array(waldse)
    ps=[2*norm.cdf(-abs(wald[_]/waldse[_]),0,1) for _ in range(0,p)]
    dfout=pandas.DataFrame({'Gene': cn, 'SingleSNPEst': wald, 'SingleSNP_P':ps})
    return dfout

# some measures of bias in SMR 
def perGeneSMRBias(bx,by,ld0,sparsePrecision,ogZs,UU,UV,cn,IVPThresh=5e-5,r2Thresh=0.25,minNIVs=10):
    # ogZs: original Z-stats in case a transformation of bx was made
    # this function uses gene groups defined by the MVMR code. For each gene, it estimates its total causal effect unconditional on
    # the other genes in the group. It then estimates the extent to which this causal estimate is biased by CHP.
    # this therefore doesn't require any additional LD calculations
    # mrbee estimates
    p=bx.shape[1]; m=bx.shape[0]; VV=numpy.array([1]).reshape((1,1))
    e0,v0,o0,k0,kiter0=imrbee(bx,by,UU,UV,VV,sparsePrecision,ld0,PleioPThreshold=0.05/m,boot=False) # strict-ish PleioPThreshold
    e0=e0[1:]; v0=v0[1:,:][:,1:] # take off intercept-relevant terms
    mask0=numpy.ones((bx.shape[0],),dtype='bool'); 
    if len(o0)>0:
        mask0[o0]=False
    _0=ld0[mask0,:][:,mask0]; _0sq=numpy.linalg.cholesky(_0).T; _0=numpy.linalg.pinv(_0)
    A1=numpy.linalg.pinv(bx[mask0,:].T@_0@bx[mask0,:]-m*UU)@bx[mask0,:].T@_0
    res0=by.reshape((by.shape[0],))-bx@e0.reshape((e0.shape[0],)); res0=res0[mask0]
    mjs=[]; jointStats=[]; jointPs=[]; pleioMeans=[]; pleioMeansP=[]; fronorm=[]; diffStats=[]; diffStatPs=[]
    uniests=[]; uniestses=[]
    q0=norm.ppf(1-IVPThresh)
    for _ in range(0,p):
        inds0=numpy.array(range(0,bx.shape[0]))
        bxj=bx[:,_].reshape((m,1))
        mask=numpy.abs(ogZs[:,_].reshape((m,1)))>q0
        mask=mask.reshape((m,))
        R11=ld0[mask,:][:,mask]
        R12=ld0[mask,:][:,mask==False]
        bxj=bxj[mask,:]; mj=bxj.shape[0]
        byj=by[mask,:]
        inds0=inds0[mask]
        if R12.shape[1]==0: # if all SNPs in bx are IVs for this gene (ie no other SNPs to compare these ones to)
            mjs.append(float('nan')); jointStats.append(float('nan')); jointPs.append(float('nan')); 
            pleioMeans.append(float('nan')); pleioMeansP.append(float('nan')); diffStats.append(float('nan')); diffStatPs.append(float('nan'))
            uniests.append(float('nan')); uniestses.append(float('nan'))
            continue
        # to make the MRBEE total causal effect estimate robust, I will remove SNPs that are in high LD with other SNPs in this region
        # recall that R12 has rows as this gene's IVs and columns as all OTHER gene's IVs
        mask2=[numpy.max((R12[x,:])**2)<r2Thresh for x in range(0,R12.shape[0])] # maks where IV not in LD with non-gene SNPs beyond r2>r2Thresh
        mask2=numpy.array(mask2)
        # if not enough data, cannot estimate anything reliably
        if sum(mask2)<minNIVs:
            if R12.shape[0]==0:
                fronorm.append(float('nan'))
            else:
                R11inv=numpy.linalg.pinv(R11)
                fronorm.append(numpy.linalg.norm(R12.T@R11inv@R11inv@R12,ord='fro'))
            mjs.append(float('nan')); jointStats.append(float('nan')); jointPs.append(float('nan')); 
            pleioMeans.append(float('nan')); pleioMeansP.append(float('nan')); diffStats.append(float('nan')); diffStatPs.append(float('nan'))
            uniests.append(float('nan')); uniestses.append(float('nan'))
            continue
        bxj=bxj[mask2,:]; mj=bxj.shape[0]
        byj=byj[mask2,:]
        UUj=numpy.array([1]).reshape((1,1))
        UVj=UV[_,:].reshape((1,1))
        R11=R11[mask2,:][:,mask2]
        R12=R12[mask2,:]
        inds0=inds0[mask2]
        # ld0, and thus R11, will always be positive definite because I already made it so in the code before this function is executed
        R11inv=numpy.linalg.pinv(R11)
        precisionjsqrt=numpy.linalg.cholesky(R11inv).T
        mjs.append(bxj.shape[0])
        # need to estimate total causal effect without bias (cannot use MVMR estimates - those are not total)
        # bootsel=True if bxj.shape[0]<30 else False
        bootsel=False # do not allow bootstrapping - only causes problems
        pleiot=0.05/bxj.shape[0]
        est0,V0,Outliers,kR,k=imrbee(bxj,byj,UUj,UVj,VV,R11inv,R11,PleioPThreshold=pleiot,boot=bootsel,niter=500) # conservative PleioPThreshold
        est0=est0[1:]; V0=V0[1:,:][:,1:] # take off intercept-relevant terms 
        uniests.append(float(est0)); uniestses.append(float(V0)**0.5)
        # Outliers are identified
        pleioPs=SpleioP(bxj,byj,UUj,UVj,VV,est0,V0)
        pleio=[(byj[x]-est0[0]*bxj[x,0])[0] for x in range(0,mj)]
        pleioPs=numpy.array(pleioPs); pleio=numpy.array(pleio)
        newmask=numpy.ones((bxj.shape[0],),dtype='bool')
        # potential problem here ... if all are outliers in IMRBEE (i.e., Outliers=list(range(0,m)), then I can't continue
        if len(Outliers)>0:
            newmask[Outliers]=False
        if sum(newmask)==0: # if all outliers
            quickAppend([jointStats,jointPs,pleioMeans,pleioMeansP,diffStats,diffStatPs])
            continue
        inds0=inds0[newmask]
        # joint test for no pleio
        V=R11+float(est0**2)*R11+float(V0)*bxj@bxj.T-2*float(est0*UVj)*R11
        Vinv=numpy.linalg.pinv(V)
        pleio=pleio.reshape((pleio.shape[0],1))
        chistat=numpy.sum(numpy.diag(Vinv@pleio@pleio.T))
        jointStats.append(chistat)
        jointPs.append(1-stats.chi2.cdf(chistat,pleio.shape[0]))
        # test for unbalanced HP
        jm=numpy.ones((bxj.shape[0],))
        Xi=R11*(numpy.eye(R11.shape[0])+float(est0*UUj)-2*float(est0*UVj))
        eta=1/len(pleio)**2*(jm.T@Xi@jm)
        stat=1/bxj.shape[0]*(jm.T@precisionjsqrt@pleio); stat=float(stat)
        statP=2*norm.cdf(-abs(stat),0,1/bxj.shape[0]**0.5)
        pleioMeans.append(float(stat)); 
        pleioMeansP.append(float(statP))
        # test for change in uni MR est vs MVMR est
        res1=pleio.copy().reshape((pleio.shape[0],)); res1=res1[newmask] # see above (non-pleio SNPs)
        # need Cov(y1|x1,y2|x2) after removing pleio
        _2=R11[newmask,:][:,newmask]; _2sq=numpy.linalg.cholesky(_2).T; _2=numpy.linalg.pinv(_2); 
        A2=numpy.linalg.pinv(bxj[newmask,:].T@_2@bxj[newmask,:]-sum(newmask)*UUj)@bxj[newmask,:].T@_2
        R01=ld0[mask0,:][:,inds0] # mask0: inds of nonpleio SNPs from MVMRBEE; inds0: inds of nonpleio SNPs from UniMRBEE
        Cov_=numpy.diag(res0)@R01@numpy.diag(res1)
        Cov_=A1@Cov_@A2.T
        Cov_=Cov_[_] # only selecting the covariance corresponding to this gene/exposure (don't care about Cov(thisGeneThetahat, otherGeneThetaHat)
        if v0.shape[0]==1:
            diff_=e0[_]-est0
            estvar=(V0+v0[_]-2*Cov_)
        else:
            diff_=e0[_]-est0
            estvar=(V0+v0[_,_]-2*Cov_)
        if estvar<0: # don't know why this happens but sometimes seems possible
            diffStats.append(float('nan'))
            diffStatPs.append(float('nan'))
            continue
        diffStat=diff_/(estvar**0.5)
        diffStatP=2*norm.cdf(-abs(diffStat),0,1)
        diffStats.append(float(diffStat))
        diffStatPs.append(float(diffStatP))
    data={'gene': cn, 'UniMRBEEEsts': uniests, 'UniMRBEESEs': uniestses, 'MRBEEEsts': flatten_list(e0.tolist()), 'MRBEESEs': (numpy.diag(v0)**0.5).tolist(),
          'SMRm': mjs, 'SMRpleioJointStats': jointStats, 'SMRpleioJointPs': jointPs, 'SMRHPmeans': pleioMeans, 'SMRHPmeansP': pleioMeansP,
         'OVBiasTestStat': diffStats, 'OVBiasTestP': diffStatPs}
    df=pandas.DataFrame(data)
    # some notes on the output
    # 'NaN' in all but last column where number of IVs<5; 'NaN' in last col if number of IVs was 0
    return df

def adjustInflation(res,infl):
    # want to adjust (Est/SE), (IVW_MVMR_Est/IVW_MVMR_SE), (IVW_UNIMR_Est/IVW_UNIMR_SE), (UniMRBEEEsts/UniMRBEESEs),
    # (MRBEEEsts/MRBEESEs), (MRBEERidgeEst/MRBEERidgeSE), (MRBEEPosSelEst/MRBEEPosSelSE)
    res2=res.copy()
    res2['MRJonesSE']=res2['MRJonesSE']*(infl**0.5)
    res2['IVW_MVMR_SE']=res2['IVW_MVMR_SE']*(infl**0.5)
    res2['IVW_UVMR_SE']=res2['IVW_UVMR_SE']*(infl**0.5)
    res2['MRBEE_UVMR_SE']=res2['MRBEE_UVMR_SE']*(infl**0.5)
    res2['MRBEE_MVMR_SE']=res2['MRBEE_MVMR_SE']*(infl**0.5)
    res2['MRBEERidge_MVMR_SE']=res2['MRBEERidge_MVMR_SE']*(infl**0.5)
    res2['MRBEEPostSelec_MVMR_SE']=res2['MRBEEPostSelec_MVMR_SE']*(infl**0.5)
    res2['MRBEEHuber_MVMR_SE']=res2['MRBEEHuber_MVMR_SE']*(infl**0.5)
    return res2

def callDelete():
    sys=platform.system()
    if (sys=='Linux') | (sys=='Darwin'):
        call_='rm'
    elif sys=='Windows':
        call_='del'
    return call_

def callPlink():
    sys=platform.system()
    if sys=='Linux':
        call_='./plinkdir/linux/plink'
    elif sys=='Darwin':
        call_='./plinkdir/mac/plink'
    elif sys=='Windows':
        call_='plinkdir/windows/plink.exe'
    return call_

def findOutcomeSignals(dataPheno,ldRefDir,writableDir,outcomeClumpingKBWindow,outcomeClumpingPthreshold,outcomeClumpingR2):
    # outcomeClumpingKBWindow, outcomeClumpingPthreshold, outcomeClumpingR2:
    #  - parameters to give PLINK to find clumps in the outcome GWAS (whole genome, not CHR-specific)
    dataPheno_=dataPheno.copy()[abs(dataPheno['phenoZ'])>norm.ppf(1-outcomeClumpingPthreshold)] # subset to speed things up
    dataPheno_['P']=numpy.array([norm.cdf(-abs(dataPheno_['phenoZ'].values[x])) for x in range(0,dataPheno_.shape[0])])
    dataPheno_[['phenoSNP','P']].rename(columns={'phenoSNP':'SNP'}).to_csv(writableDir+'/outcomePsout.txt', index=False, sep="\t")
    cmd=[callPlink(), "--bfile", ldRefDir, "--clump", writableDir+"/outcomePsout.txt", "--clump-kb", str(outcomeClumpingKBWindow), "--clump-p1", str(outcomeClumpingPthreshold), "--clump-p2", str(outcomeClumpingPthreshold), "--clump-r2", str(outcomeClumpingR2), "--out", writableDir+"/outcomeplinkout"]
    out=subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) # executing the command in the terminal
    outcomeClumps=pandas.read_fwf(writableDir+'/outcomeplinkout.clumped', sep='\t') # reading LD matrix back in
    outcomeClumps=outcomeClumps['SNP'].dropna().values # for some reason some NAs get appended to the end
    dataPheno['isOutcomeClump']=dataPheno['phenoSNP'].isin(outcomeClumps)
    return dataPheno

def geneGroupFinder(geneGroups, gene, isGtex=False):
    keynames=list(geneGroups)
    # if from GTEx, gene IDs will have a "." for some reason ... 
    found=0
    for s in range(0,len(keynames)):
        fl=geneGroups[keynames[s]]
        if isGtex:
            fl=[fl[i].split('.')[0] for i in range(0,len(fl))]
            gene=gene.split('.')[0]
        if gene in fl:
            found=1
            break
    if found==1:
        out=keynames[s]
        return out # returns keyname of group the gene is in
    else:
        return 'cannot find'

##################################################################################################################################################
# so-called upper-level functions

# function to check that files exist
def fileChecker(fp,ftype):
    if os.path.exists(fp)==False:
        raise ValueError('The {} does not exist. Are you sure you entered it correctly?'.format(ftype))

# function to load outcome data and determine the chromosome
def loadOutcomeGWASData(fp,effectAllele,z,rsid,ldRefDir):
    if '@' in z: # if no Z-stat column in phenotype data
        sep=z.split('@',1) # user must have passed 'BETA@SE' for Z-stat column name
        uc2=[rsid,sep[0],sep[1],effectAllele] # only pick out these columns from pheno GWAS data
        dataPheno=pandas.read_table(fp,sep=None,engine='python',usecols=uc2) # load pheno GWAS data
        dataPheno['phenoZ']=dataPheno[sep[0]]/dataPheno[sep[1]] # create my own Z-stat column
        dataPheno=dataPheno.rename(columns={uc2[0]: 'phenoSNP', uc2[3]: 'phenoEffectAllele'}) # rename to names I'll now later
        dataPheno=dataPheno[['phenoSNP','phenoZ','phenoEffectAllele']] # drop unneccessary columns
    else:
        uc2=[rsid,z,effectAllele]
        dataPheno=pandas.read_table(fp,sep=None,engine='python',usecols=uc2)
        dataPheno=dataPheno.rename(columns={uc2[0]: 'phenoSNP', uc2[1]: 'phenoZ', uc2[2]: 'phenoEffectAllele'})
    # filter to only LD ref dir SNPs
    onekgSNPs=pandas.read_table(ldRefDir+'.bim',sep=None,engine='python',header=None) # has column names CHR, SNP (SNP is for rsID)
    onekgSNPs=onekgSNPs[[1]] # the SNP column name is actually an index for .bim files 
    merged=pandas.merge(dataPheno,onekgSNPs,left_on='phenoSNP',right_on=1)
    merged=merged[['phenoSNP','phenoEffectAllele','phenoZ']]
    return merged # this is dataPeno and 1kg merged; NOT exposure and pheno/outcome merged

def rsidFinder(st):
    # x should be a string
    sp=st.split(':')
    if len(sp)==1:
        return st
    w=['rs' in sp[_].lower() for _ in range(0,len(sp))]
    if True in w:
        out=sp[w.index(True)]
    else:
        out=float('nan')
    return out

def prepRes(res):
    res=res[['Gene','geneBP','Est','MRBEEPosSelEst','MRBEEPosSelSE','MRJonesR2','geneZ','geneSNP','snpBP','conditionalF','groupID','CHR','infl','h2CisEstNaive','nClumps','SingleSNPEst','SingleSNP_P','IVW_MVMR_Est','IVW_MVMR_SE','IVW_UNIMR_Est','IVW_UNIMR_SE','UniMRBEEEsts','UniMRBEESEs','MRBEEEsts','MRBEESEs','MRBEERidgeEst','MRBEERidgeSE','MRBEEHuberEst','MRBEEHuberSE','SMRm','SMRpleioJointPs','SMRHPmeansP','OVBiasTestP']]
    res=res.rename(columns={'Est': 'MRJonesEst', 'MRJonesR2':'RsquaredMRJones','geneZ': 'LeadEQTLZstatistic','geneSNP':'leadEQTL','snpBP':'leadEQTLBP','groupID':'CHRspecificGroupID','CHR':'Chromosome','infl':'nullInflation','IVW_UNIMR_Est':'IVW_UVMR_Est','IVW_UNIMR_SE':'IVW_UVMR_SE','UniMRBEEEsts':'MRBEE_UVMR_Est','UniMRBEESEs':'MRBEE_UVMR_SE','MRBEEEsts':'MRBEE_MVMR_Est','MRBEESEs':'MRBEE_MVMR_SE','MRBEERidgeEst':'MRBEERidge_MVMR_Est','MRBEERidgeSE':'MRBEERidge_MVMR_SE','MRBEEHuberEst':'MRBEEHuber_MVMR_Est','MRBEEHuberSE':'MRBEEHuber_MVMR_SE','MRBEEPosSelEst':'MRBEEPostSelec_MVMR_Est','MRBEEPosSelSE':'MRBEEPostSelec_MVMR_SE','SMRm':'SMRNumIVs','SMRHPmeansP':'SMRUnbalancedHP_P','OVBiasTestP':'SMROVBias_P'})
    return res

# function to load exposure data and merge with outcome data (separated from outcome data bc many exposure files will be loaded)
# with specifying sample size column name
def loadExposureGWASData(fp,effectAllele='',z='',rsid='',snpBP='',geneLabel='',geneBPLabel='',nLabel='',isRawGzippedGTEx=False,mapwd='data/maps/1kgPhase3maps',nr=1e100): # outcome data already filtered to 1kg SNPs
    dataLoaded=False;
    # If this is raw GTEx data, I already know what to expect for column names 
    if isRawGzippedGTEx==True: # I take care of all column names here
        from io import BytesIO # Converting bytes to bytes input file
        import gzip
        import pyarrow # Fast reading of parquets
        file=gzip.open(fp, 'rb')
        # speed of loading full parquet data seems to be same as loading some subset
        # pf=pyarrow.parquet.ParquetFile(file);first_ten_rows=next(pf.iter_batches(batch_size=10));df=pa.Table.from_batches([first_ten_rows]).to_pandas()
        df=file.read()
        df=pandas.read_parquet(BytesIO(df)); dataLoaded=True
        dataGene=attachRsidToGTEx(df,mapwd=mapwd)
        dataGene['geneZ']=dataGene['slope']/dataGene['slope_se']
        dataGene=dataGene.rename(columns={'rsID': 'geneSNP', 'phenotype_id': 'Gene', 'ALT': 'geneEffectAllele'}) # already has snpBP column from attachRsidToGTEx()
        if is_number(nLabel):
            dataGene['geneN']=nLabel
        else: 
            dataGene=dataGene.rename(columns={nLabel: 'geneN'})
        dataGene['tempGeneBP']=dataGene['snpBP']-dataGene['tss_distance'] # snpBP is in hg38 since it is from GTEx
        aggd=dataGene.groupby('Gene').agg({'tempGeneBP': 'mean'})
        aggd=pandas.DataFrame({'geneBP':aggd.tempGeneBP.values,'Gene':list(aggd.index)})
        dataGene=pandas.merge(dataGene,aggd,left_on='Gene',right_on='Gene')
        dataGene=dataGene[['geneSNP','geneZ','Gene','geneEffectAllele','geneBP','snpBP','geneN']] # drop unneccessary columns
    else: # I need to take care of the column names now
        if '@' in z:
            sep=z.split('@',1) # user must have passed 'BETA@SE' for Z-stat column name
            uc1=[rsid,sep[0],sep[1],geneLabel,effectAllele,geneBPLabel,snpBP] # only pick out these columns from exposure GWAS data
            if dataLoaded==False:
                dataGene=pandas.read_table(fp,sep=None,engine='python',nrows=nr); # load exposure GWAS data
            dataGene['geneZ']=dataGene[sep[0]]/dataGene[sep[1]] # create my own Z-stat column
            dataGene=dataGene.rename(columns={uc1[0]: 'geneSNP', uc1[3]: 'Gene', uc1[4]: 'geneEffectAllele', uc1[5]: 'geneBP',uc1[6]: 'snpBP'})
            if is_number(nLabel):
                dataGene['geneN']=nLabel
            else: 
                dataGene=dataGene.rename(columns={nLabel: 'geneN'})
            dataGene=dataGene[['geneSNP','geneZ','Gene','geneEffectAllele','geneBP','snpBP','geneN']] # drop unneccessary columns
        else:
            uc1=[rsid,z,geneLabel,effectAllele,geneBPLabel,snpBP]
            if dataLoaded==False:
                dataGene=pandas.read_table(fp,sep=None,engine='python')
            dataGene=dataGene.rename(columns={uc1[0]: 'geneSNP', uc1[1]: 'geneZ', uc1[2]: 'Gene', uc1[3]: 'geneEffectAllele', uc1[4]: 'geneBP',uc1[5]: 'snpBP'})
            if is_number(nLabel):
                dataGene['geneN']=nLabel
            else: 
                dataGene=dataGene.rename(columns={nLabel: 'geneN'})
            dataGene=dataGene[['geneSNP','geneZ','Gene','geneEffectAllele','geneBP','snpBP','geneN']] # drop unneccessary columns
    # if data is from MetaBrain, change SNP column to only contain rsIDs
    dataGene['geneSNP']=dataGene['geneSNP'].apply(lambda x: rsidFinder(x))
    return dataGene

# merge exposure and outcome data
def mergeExposureAndOutcomeGWAS(dataGene,dataPheno):
    merged=pandas.merge(left=dataGene, right=dataPheno, left_on='geneSNP', right_on='phenoSNP')
    # harmonise exposure and outcome effect sizes by effect alleles
    mask=merged.geneEffectAllele.str.upper()!=merged.phenoEffectAllele.str.upper() # find where [uppercase] alleles do not match
    merged.loc[mask,'geneZ']=(-1*merged.loc[mask,'geneZ']) # change sign of Z-stat for gene expression if effect alleles do not match
    return merged

# function to attach rsIDs to GTEx data that is in hg38 and SNPs are named as eg chr9_10536_C_T_b38
def attachRsidToGTEx(dataGTEx,mapwd='data/maps/1kgPhase3maps'):
    # the mapfile is in a pre-determined location that will be set when the user clones our github repo:
    # alt allele in GTEx is effect allele
    # separate information in rsidcolumnname of dataGTEx
    splitted=dataGTEx['variant_id'].str.split('_',n=10)
    splitted=numpy.array(splitted)
    marker_=[]; chr_=[]; bp_=[]; alt_=[]; ref_=[]
    for _ in range(0,splitted.shape[0]):
        chr_.append(int(splitted[_][0].replace('chr','')))
        bp_.append(int(splitted[_][1]))
        ref_.append(splitted[_][2])
        alt_.append(splitted[_][3])
        marker_.append(str(chr_[_])+'_'+str(bp_[_]))
    dataGTEx['markername']=marker_
    dataGTEx['CHR']=chr_
    dataGTEx['REF']=ref_ # GTEx and merged_ will now have these columns, too
    dataGTEx['ALT']=alt_ # effect allele according to documentation
    dataGTEx['snpBP']=bp_
    # load mapfile corresponding to this chromosome (should all be one CHR at a time from dataGene so using any values in chr_ will work)
    # names of mapfiles do not change and are not set by the user
    mapfile=pandas.read_csv(mapwd+'/'+'1kgPhase3hg19Chrom'+str(chr_[0])+'MapToHg38.txt.gz',engine='python',sep=None) # uses chromosome-specific mapfiles. chromsome determined using data
    mapfile['markername']=mapfile['CHROM_b38'].astype(str)+'_'+mapfile['POS_b38'].astype(str) # make out of hg38 values bc that's what's in GTEx
    mapfile=mapfile[['rsID','markername']] # drop unnecessary columns
    # merge 
    merged_=pandas.merge(dataGTEx,mapfile,left_on='markername',right_on='markername')
    return merged_

### !!! NOTE !!! this function is pretty slow ... loading phenotype data always seems to take a long time 
def attachRsidToPhenoData(fpdata,chrom,bp,effectAllele,z,rsid,ldRefDir,build='NA',mapwd='data/maps/fullMaps/1kgPhase3MapAllSNPs.txt',nSNPsToTest=1000):
    if (is_number(build)==False) & (build.lower()!='na'):
        raise ValueError('please enter either 37 (hg19) or 38 (hg38) as the build parameter')
    data=pandas.read_csv(fpdata,sep=None,engine='python'); 
    mapdf=pandas.read_csv(mapwd,sep=None,engine='python')
    # all mapfiles should have rsID  CHROM_b37  CHROM_b38   POS_b37   POS_b38 (I know this because I will be the one to create them)
    if build.lower()=='na': # if the user does not know the build
        build=determineBuild(data,chrom,bp,mapdf,nSNPsToTest=nSNPsToTest)
    mapdf=mapdf.rename(columns={'CHROM_b'+str(build): 'usechr', 'POS_b'+str(build): 'usebp'}) # identifying columns in map corresponding to build
    mapdf['usemarker__']=mapdf['usechr'].astype(str)+':'+mapdf['usebp'].astype(str) # making markername column to merge map on
    data['usemarker__']=data[chrom].astype(str)+':'+data[bp].astype(str) # making markername column to merge data on
    data=data.drop_duplicates(subset=['usemarker__']) # cannot have duplicates when merging
    dataPheno=pandas.merge(data,mapdf[['rsID','usemarker__']],left_on='usemarker__',right_on='usemarker__') # merge data and map
    nMapped=dataPheno.shape[0] # number of mapped SNPs
    dataPheno=dataPheno.drop(columns=['usemarker__']) # dropping markername column I made
    del data # to save memory
    # note that any SNPs in data that are not in the map are dropped - you should therefore choose a map consistent with the LD ref panel since all 
    # SNPs will be subsetted to those in the LD ref panel anyway
    # 
    # instead of saving, I should now give it to the user in a cleaned way as I would if they used loadPhenoData()
    if '@' in z: # if no Z-stat column in phenotype data
        sep=z.split('@',1) # user must have passed (eg) 'BETA@SE' for Z-stat column name
        uc2=[rsid,sep[0],sep[1],effectAllele] # only pick out these columns from pheno GWAS data (rsid should be rsid in the map. if i made the map then it is 'rsID')
        dataPheno['phenoZ']=dataPheno[sep[0]]/dataPheno[sep[1]] # create my own Z-stat column
        dataPheno=dataPheno.rename(columns={uc2[0]: 'phenoSNP', uc2[3]: 'phenoEffectAllele'}) # rename to names I'll now later
        dataPheno=dataPheno.rename(columns={'rsID': 'phenoSNP'})[['phenoSNP','phenoZ','phenoEffectAllele']] # drop unneccessary columns. mapdf always attacheds 'rsID'
        dataPheno['phenoEffectAllele']=dataPheno['phenoEffectAllele'].str.upper() # make effect allele uppercase
    else:
        uc2=['rsID',z,effectAllele]
        dataPheno=dataPheno.rename(columns={uc2[0]: 'phenoSNP', uc2[1]: 'phenoZ', uc1[2]: 'phenoEffectAllele'})
        dataPheno['phenoEffectAllele']=dataPheno['phenoEffectAllele'].str.upper() # make effect allele uppercase
    # filter to only LD ref dir SNPs
    onekgSNPs=pandas.read_table(ldRefDir+'.bim',sep=None,engine='python',header=None) # has column names CHR, SNP (SNP is for rsID)
    onekgSNPs=onekgSNPs[[1]] # the SNP column name is actually an index for .bim files 
    merged=pandas.merge(dataPheno,onekgSNPs,left_on='phenoSNP',right_on=1)
    merged=merged.drop(columns=1)
    return merged

def determineBuild(data,chrom,bp,mapdf,nSNPsToTest=1000):
    # mapdf always has columns rsID  CHROM_b37  CHROM_b38   POS_b37   POS_b38
    mapdf['map__marker37']=mapdf['CHROM_b37'].astype(str)+':'+mapdf['POS_b37'].astype(str) # funny name so nobody else has it picked already
    mapdf['map__marker38']=mapdf['CHROM_b38'].astype(str)+':'+mapdf['POS_b38'].astype(str)
    mask=data[chrom].values==mapdf['CHROM_b37'].values[0]
    datacut=data[mask] # data will be for all CHRs since it is outcome data, but map data will be CHR-specific. This is fine, just subset outcome data this this CHR
    datacut=datacut.iloc[:(nSNPsToTest*100),:] # assuming only 1/100th of outcome GWAS SNPs will be in map
    datacut['data__marker']=datacut[chrom].astype(str)+':'+datacut[bp].astype(str)
    m1=pandas.merge(datacut,mapdf,left_on='data__marker', right_on='map__marker37') # for build 37
    m2=pandas.merge(datacut,mapdf,left_on='data__marker', right_on='map__marker38')  # for build 38
    if (m1.shape[0]/2)>m2.shape[0]:
        out=37
    elif m1.shape[0]<(m2.shape[0]/2):
        out=38
    else:
        out='cannot reliably determine build - please increase nSNPsToTest until this message disappears'
    return out

# identify gene groups for multivariable MR
def defGeneGroups(q0geneGroups, merged, allowSingletons=False):
    # merged=merged.sort_values(by='geneBP',axis=0) # sort by BP position of genes - not sure why this matters but it does
    genes=merged['Gene'].unique().tolist(); # genes=list(genes)
    mergedSig=merged[abs(merged['geneZ'])>q0geneGroups]
    usedGenes=[]; lens=[]; geneGroups={}
    for ii in range(0, len(genes)):
        if (genes[ii] in usedGenes)==True:
            continue
        else:
            qts=mergedSig[mergedSig['Gene']==genes[ii]]['geneSNP'].values # all SNPs that are cis-eQTLs for this gene (genes[ii])
        if len(qts)==0: # if no cis-eQTLs, skip gene
            continue
        qtsGenes=mergedSig[mergedSig['geneSNP'].isin(qts)]['Gene'].unique() # all genes for which these SNPs (qts) are cis-eQTLs
        if (len(qtsGenes)==1) & (allowSingletons==False): # I will not form a group of 1 ... why not?
            continue
        toAdd=list(qtsGenes) # genes sharing one of these eQTLs
        boo=[(toAdd[x] in usedGenes)==False for x in range(0,len(toAdd))] # boolean: True=gene has not been grouped; False=gene has been grouped already
        if sum(boo)==0: # If all genes sharing this eQTL are already in another group, next
            continue
        else:
            toAdd=(numpy.array(toAdd)[boo]).tolist() # else, only use the ungrouped genes to form a group with the current gene
        if (allowSingletons==False) & (sum(boo)==1):
            continue
        usedGenes.append(toAdd); usedGenes=flatten_list(usedGenes)
        geneGroups[genes[ii]]=toAdd # geneGroups is a dictionary with names that are the index genes in a group
        lens.append(len(toAdd))
    ggKeys=[list(geneGroups)[x] for x in range(0,len(geneGroups))]
    return geneGroups, ggKeys, lens, usedGenes

def defGeneGroupsByOutcome(q0geneGroups, merged, KbWindow=100,closestK=2):
    # assumes the "isOutcomeClump" column is in merged
    # KbWindow should be set to the same, or smaller, value as the value for outcomeClumpingKBWindow
    clumpBPs=merged[merged['isOutcomeClump']==True]['snpBP'].unique()
    usedGenes=[]; lens=[]; geneGroups={}
    mergedcut=merged.copy()[abs(merged['geneZ'])>q0geneGroups] # subsetting only to significant eQTLs
    for i in range(0,len(clumpBPs)):
        cutdf=mergedcut.copy()[abs(mergedcut['snpBP']-clumpBPs[i])<(KbWindow*1e3)] # subsetting to only those within KbWindow either side of outcome clump
        cutdf['absDist']=abs(cutdf['snpBP'].values-clumpBPs[i]) # calculating absolute distance from eQTL to outcome clump
        closestKgenes=(cutdf['Gene'].unique())[:closestK] # it puts them in order of appearance in data, which is now sorted ascending by distance from outcome clump
        # now find genes sharing eQTLs with any of these closestK genes closest to the clump
        theseeqtls=cutdf[cutdf['Gene'].isin(closestKgenes)]['geneSNP'] # find eQTLs for closest closestKgenes genes
        genestogroup=cutdf[cutdf['geneSNP'].isin(theseeqtls)]['Gene'].unique() # this is the group that is formed
        boo=[(genestogroup[x] in usedGenes)==False for x in range(0,len(genestogroup))] # True if not grouped, False otherwise
        if sum(boo)<1: # this will allow singletons, but will skip empty ones
            continue
        genestogroup=genestogroup[boo] # only group those not previously grouped 
        geneGroups['group'+str(i)]=genestogroup # key is a non-informative number - else there can be issues for some reason
        usedGenes.append(genestogroup.tolist()) # add newly grouped genes to running list of those already used
        usedGenes=flatten_list(usedGenes)
    ggKeys=list(geneGroups)
    lens=[len(geneGroups[ggKeys[x]]) for x in range(0,len(ggKeys))]
    return geneGroups,ggKeys,lens,usedGenes

def defineCandidateGeneGroups(merged,candidateGenes,MbWindow=2):
    geneGroups={}; lens=[]; usedGenes=[]
    for _ in range(0,len(candidateGenes)):
        gbp=merged[merged['Gene']==candidateGenes[_]]['geneBP'].values[0]
        aroundgenes=merged[abs(merged['geneBP']-gbp)<(MbWindow*1e6)]['Gene'].unique()
        geneGroups[candidateGenes[_]]=aroundgenes
        usedGenes.append(usedGenes)
        lens.append(len(aroundgenes))
    ggKeys=list(geneGroups)
    usedGenes=flatten_list(usedGenes)
    return geneGroups,ggKeys,lens,usedGenes

# perform analysis
def MVMRworkhorse(merged,geneGroups,ggKeys,writableDir,ldRefDir,isGtex=False,
                  analysisOnlyInOutcomeLoci=True,outcomeLociMbWindow=1,ldUpperLimit=0.5,
                  ldOtherLociOtherPt=0.0001,ldOtherLociR2=0.1,ldOtherLociWindow=1,q0Correls=3.290527,nMinCorrels=50,
                  jointChiGenesP=5e-8,opAlpha='dynamic',nMinIVs=50,hessMinScale=5,silence=False,
                  UniMRIVPThreshold=1e-8,candidateGenes='',assumeNoSampleOverlap=True,shrinkBiasCorrection=True,networkR2Thres=0.25,impute=True,saveData=False):
    # ldOtherLociOtherPt: only SNPs with P<this threshold will be considered when removing SNPs in IV set that are in LD with SNPs outside of the LD set
    # making sure the user hasn't passed anything crazy for ldOtherLociWindow
    if (ldOtherLociWindow>10) & (silence==False):
        raise ValueError("are you sure you gave `ldOtherLociWindow` the right value? This is {} megabases! If so, set the `-silence` flag to `true` to silence this message next time".format(str(ldOtherLociWindow*1e6)))
    # merged: full data set of merged gene expression and phenotype GWAS data - already went through QC (e.g., harmonizing alleles)
    # geneGroups: dictionary with list elements of genes that belong together in groups
    # ggKeys: list(geneGroups) - ie names of keys of geneGroups dictionary
    # writableDir: a directory that the user can write in - usually the current WD from which this software is executed
    # ldRefDir: filepath to LD reference panel (without extension - must have .bed, .fam, bim accompanying files)
    # analysisOnlyInOutcomeLoci: boolean, if true, will only perform MR within/near outcome clumps 
    # outcomeLociMbWindow: float; if analysisOnlyInOutcomeLoci==True, will only perform MR for groups with any genes within outcomeLociMbWindow Mb of outcome an clump
    # ldUpperLimit: for 2 SNPs with LD>this threshold, one will be randomly dropped
    # q0Correls: Only SNPs with P<this quantile will be considered when estimating correlations between measurement errors
    # nMinCorrels: If <this number of nonsig and nonmissing SNPs are avaiable to calcualte correls b/w measurement errors, set as missing
    # jointChiGenesP: only SNPs with P<this threshold in a p-degree of freedom joint test for p genes will be included as IVs
    # assumedMissingMean: assume that SNPs with missing association estimates (Z-stats) have this value as their mean
    # opAlpha: alpha to use in MR-Jones (note no grid searching is currently allowed so as to keep the program fast), or 'dynamic' (chosen from correlations; <1)
    # verbose: True/False - should progress be printed
    # nMinIVs: minimum number of IVs to use in MR; else group will be skipped
    # hessMinScale: diag(bx.T@LDinv@bx)/diag(m*Suu) must have be greater than hessMinScale
    # silence: a warning if it seems like the user set the ldOtherLociWindow to be extremely large
    # UniMRIVPThreshold: P-value threshold for univariable MR using the perGeneSMRBias() function
    ###
    ### want to record positions of outcome clumps up front
    # merged['Gene']=merged['Gene'].str.split(pat='.',n=0,expand=True).iloc[:,0] # dropping the '.<x>' suffixes from gene IDs (applies to GTEx, for example)
    outcomeclumpbps=merged[merged['isOutcomeClump']==True]['snpBP'].unique() # workhorse() receives CHR-specific data, so just need BP position
    ### first, if user just wants to test a set of candidate genes, make sure this CHR has some of them, else next
    if len(candidateGenes)>0:
        analysisOnlyInOutcomeLoci=False # perform analysis in candidateGenes loci
        tm=pandas.merge(merged,pandas.DataFrame({'Gene': candidateGenes}),left_on='Gene',right_on='Gene')
        if tm.shape[0]==0:
            del tm
            print('Candidate Genes not found in data')
            return {}, {}, {} # exit function and move to next CHR
    ###
    progs=list(numpy.linspace(0,len(geneGroups),10)); progs=[int(progs[x]) for x in range(0,len(progs))] # for progress printing
    thingsMonitored={}; outerDict={}; edgeDict={}
    for ogene in range(0, len(geneGroups)): # 0, len(geneGroups)
        thingsToMonitor={}
        ############################################################### start top 1/3rd (preparing data) (avg ~7 seconds per gene group)
        # 1) subsetting full data to only selected genes
        theseGenes=geneGroups[ggKeys[ogene]] # identify which genes are in this group
        if isinstance(theseGenes,str):
            theseGenes=[theseGenes] # if only one gene considered
        # if user only wants to test candidate genes, restrict to only those groups containing at least one of the candidate genes
        if len(candidateGenes)>0:
            tt=numpy.array([theseGenes[x] in candidateGenes for x in range(0,len(theseGenes))])
            if sum(tt)==0:
                continue
        
        thingsToMonitor['startingNumGenes']=len(theseGenes) # number of starting genes
        genedf=merged[merged['Gene'].isin(theseGenes)] # filter full data set to only these genes
        genedf=genedf.dropna(subset=['phenoZ']) # drop rows with missing values for outcome (I will not want to impute these)
        genedf=genedf.drop_duplicates(subset=['geneSNP','Gene'],keep='first') # cannot keep SNPs tested >once for a gene 
        sn=list(genedf['geneSNP'].unique()) # list of unique SNPs present in data for these genes (will be number of rows in wide data matrix)
        ### if the user wants to only look in significant outcome loci (see arguments analysisOnlyInOutcomeLoci and outcomeLociMbWindow)
        if (analysisOnlyInOutcomeLoci==True):
            sbps=genedf['snpBP'].unique()
            tt=[any(abs(sbps-outcomeclumpbps[x])<(outcomeLociMbWindow*1e6)) for x in range(0,len(outcomeclumpbps))]
            if sum(tt)==0: # if no genedf SNPs (candidate IVs) are within outcomeLociMbWindow of any outcome clumps, move to the next group
                thingsToMonitor['locusSkippedBCNoOutcomeSignal']=ggKeys[ogene] # defining the gene locus using the index gene
                continue
        
        ### may still want to impose a restriction on which SNPs can be included - some SNPs have >75% missing for genes - do we really want to keep these?
        gnKeep=genedf['Gene'].value_counts(); gnKeep=gnKeep[gnKeep>numpy.min([100,len(sn)*0.1])] # restricting genes to only those with >100 SNPs or at least 10% observed
        genedf=genedf[genedf['Gene'].isin(list(gnKeep.index))]; # restricting and noting which genes got removed
        newGenes=genedf['Gene'].unique() # list of theseGenes cut to those with a enough observed data to continue
        lostGenes=[theseGenes[x] for x in range(0,len(theseGenes)) if (theseGenes[x] in newGenes)==False] # which genes got removed?
        if len(newGenes)==1: # if all but one gene just got removed
            reshaped=genedf.copy(); reshaped=reshaped.rename(columns={'geneZ': newGenes[0]+'_Z'})
        else:
            reshaped=reshapeCorrelDf(genedf,newGenes,'Gene',joinType='outer') # reshaping data to be wide format (will have some NaN's if joinType='outer')
        
        ### identify lead SNP for each gene (its ok if some of these will not be included in MR, this is for later infererence and plotting)
        genedf['absZ']=numpy.abs(genedf['geneZ'])
        leadGeneSNPDf=genedf.loc[genedf.groupby('Gene')['absZ'].idxmax()][['geneSNP','Gene','geneZ']]
        genedf=genedf.drop(columns=['absZ'],axis=1)
        cn=[newGenes[x]+'_Z' for x in range(0,len(newGenes))] # column names for Z-statistics 
        thingsToMonitor['genesDroppedBCMissing']=len(theseGenes)-len(cn) # how many genes got dropped because of missing?
        thingsToMonitor['startingNSNPs']=reshaped.shape[0] # starting number of SNPs
        geneLocs=genedf[['Gene','geneBP']].drop_duplicates() # record locations of all genes in this group now before any get dropped
        ############################################################### end top 1/3rd (preparing data) (avg ~7 seconds per gene group)
        ############################################################### start middle 1/3rd (imputation, measurement error correlations, IV selection, LD estimation) (avg ~7 seconds per gene group)
        ### Imputation of missing in actual gene exposure data matrix
        # Yihe's reduced rank matrix completion (imputation with 0s or GAM predictions can both be problematic)
        X=reshaped[cn].values; # actually a lot of missing if the number of genes in a group is very large
        if len(cn)==1: # if only one gene, still want this to be a matrix
            X=X.reshape((X.shape[0],1)) 
        # cn and columns of X are in the same order. How to find geneBPs that are in this same order:
        cnstripped=[cn[x].split('_',1)[0] for x in range(0,len(cn))]
        gg=pandas.merge(pandas.DataFrame({'hi': cnstripped}), geneLocs, left_on='hi', right_on='Gene')['geneBP'].values
        mask=numpy.isnan(X); Del=mask+1-1; newX=X.copy()
        if (X.shape[1]>1) & (impute==True):
            L,Factor,Gamma=ffimputeBIC(X,Kmax=15,dg=gg,dim_proj=0,W=makeI(X.shape[0]),iter_max=15,eps=0.001,tuning=0)
            newX[mask]=L[mask]
        elif (X.shape[1]>1) & (impute==False):
            newX[mask]=0
        
        # Omega=numpy.isnan(X); Xc=X.copy(); Xc[Omega]=assumedMissingMean; d0=numpy.linalg.svd(Xc,compute_uv=False)
        # newX,out2,out3=fullRankMatrixCompletion(X.copy(),missImp=0,lambda_=min(d0)/2,eps=0.01*X.shape[1]) # lambda_ chosen somewhat arbitrarily
        mm=numpy.isnan(X).sum(axis=0)/X.shape[0]
        thingsToMonitor['minPropMissingZs']=numpy.min(mm)
        thingsToMonitor['maxPropMissingZs']=numpy.max(mm)
        thingsToMonitor['minNumNonMissing']=X.shape[0]-numpy.max(numpy.isnan(X).sum(axis=0))
        # probably need to do some more work on the method in fullRankMatrixCompletion(), but I think generally the idea is appropriate.
        ### calculation of R=DSD where S=SigmaUU
        ccDf=reshaped[cn]; # `reshaped` with only columns with names in `cn`. Must reset index
        copycat=(ccDf.copy()).reset_index(); 
        ccDf=ccDf.assign(myInd=list(copycat.index)) # add column to end of ccDf indicating row number (indices get all messed up)
        cond=numpy.array(abs(ccDf)>q0Correls);
        # for each pair, find all SNPs where cond==False for both, which may ==True for others
        p=ccDf.shape[1]-1; cc=makeI(p); nn=numpy.zeros((p,p))
        for i in range(0,p):
            for j in range(0,p):
                if i>j or i==j:
                    continue
                else:
                    condi=flatten_list(cond[:,i]); condj=flatten_list(cond[:,j]) # finding logical for column i and j
                    keepijs=[x for x in range(0,len(condi)) if (condi[x]==False and condj[x]==False)] # want where both nonsignificant
                    subCCDF=ccDf.copy(); subCCDF=subCCDF.iloc[keepijs,[i,j,p]].dropna() # and remove missing; p+1 is to keep the row number column I made
                    mask=numpy.zeros(len(condi),dtype='bool'); mask[subCCDF['myInd'].values]=True
                    subCCDF=subCCDF.drop('myInd', axis=1)
                    nn[i,j]=subCCDF.shape[0]
                    if nn[i,j]<nMinCorrels: # I previously chose 50 somehwat arbitrarily
                        # print('problem with '+str(i)+' and '+str(j)) # probably want to comment this out
                        cc[i,j]=float('nan') # make missing - to be filled in with MaxDet imputation
                        continue
                    # correlation using only nonsig
                    # recall indexing by rows and columns must happen sequentially not simultaneously
                    cc[i,j]=numpy.corrcoef(subCCDF.values,rowvar=False)[0,1] # no need to premultiply by LD matrix bc since nonsig assume not affected by LD
                    # correlation using all with resampling
                    # cc[i,j]=rs,vstotal=subSamCalcR(subCCDF,k=subCCDF.shape[0])
        
        cc=cc+cc.T; numpy.fill_diagonal(cc,1)
        del ccDf
        Omega=numpy.isnan(cc)
        thingsToMonitor['numMissingMECorrelations']=sum(Omega.reshape((Omega.shape[0]*Omega.shape[1],)))/2 # number of missing values in correlation matrix
        ### imputation of missing in R=DSD where S=SigmaUU
        # using MaxDet to fill in missing & drop any genes that you need to (where MaxDet sol not in {-1,-1} (if no missing this code can still run and nothing will happen)
        ccComplete,ccDropped,ccNDropped=MaxDet(cc)
        mask=numpy.ones((cc.shape[0],),dtype='bool') # creating a mask to drop genes/exposures (see below)
        mask[ccDropped]=False # if gene was dropped bc of a bad MaxDet solution, do not consider it below
        Omega=Omega[mask,:][:,mask] # find indices of missing for non-bad MaxDet solutions (ie, find which were successfully imputed)
        cnMask=numpy.ones(len(cn),dtype='bool')
        cnMask[ccDropped]=False
        newX=newX[:,mask]
        cn=numpy.array(cn)[cnMask].tolist()
        thingsToMonitor['numBadMaxDetSolutions']=cc.shape[0]-ccComplete.shape[0] # number of genes dropped bc of bad MaxDet solutions
        # do not drop any SNPs - this may change in the future, but for now do not
        ccComplete,alpha,deco=posDefifyCorrMat(ccComplete,epsilon=1e-5) # making ccComplete positive definite via shrinkage if it is not
        thingsToMonitor['alphaToMakeMECorrelMatPosDef']=alpha # factor that off-diagonal elements of ccComplete needed to multiplied by to make it posdef
        # remove genes with norm(newX[:,j],2)^2<5m
        toKeep=(numpy.diag(newX.T@newX)>(hessMinScale*newX.shape[0]/2)) # technically >hessMinScale*newX.shape[0]*diag(Suu), but diag(Suu)=(1)
        if sum(toKeep)==0:
            continue
        newX=newX[:,toKeep]
        ccComplete=ccComplete[toKeep,:][:,toKeep]
        droppedGenes=numpy.array(cn)[toKeep==False].tolist()
        cn=numpy.array(cn)[toKeep].tolist()
        reshaped=reshaped.drop(labels=droppedGenes,axis=1)
        # perform p-degree of freedom chi-square test to test SNP assoc with at least one gene at P<tau. drop SNPs not meeting this criteria
        # Note I need SigmaUU, but not the LD matrix, for this so it happens here and not earlier
        tau=jointChiGenesP; # no need to pre-multiply by sparsePrecisionSqrt bc I am considering one SNP at a time
        chi_=numpy.diag(newX@numpy.linalg.inv(ccComplete.copy())@newX.T) # the inverse matrix is the inverse of SigmaUU
        chiPs=[1-stats.chi2.cdf(chi_[x],newX.shape[1]) for x in range(0,newX.shape[0])]
        mask=numpy.array(chiPs)<tau; mask=mask.reshape((mask.shape[0],))
        reshaped=reshaped[mask];
        newX=newX[mask,:]
        if newX.shape[0]<nMinIVs:
            continue # also need to go to next if not many SNPs here
        
        # 3) retrieving LD matrix for selected genes
        # `reshaped` contains the SNPs I want to use in MR (`ccDf` was only used for calculating the correlation matrix)
        # generating SNP LD matrix
        reshaped=reshaped.reset_index()
        reshaped=reshaped.drop_duplicates(subset=['geneSNP'],keep=False) # drop duplicate SNPs (cannot appear in file given to PLINK)
        mask=numpy.zeros((newX.shape[0]),dtype='bool'); mask[reshaped.index.tolist()]=True
        newX=newX[mask,:] # in case any duplicates got removed
        outDir1=writableDir+'/myExtract.txt' # (for saving SNPs for LD calculations) write out to a directory that must be writable (since the user is there)
        outDir2=writableDir+'/tempOut' # (for saving LD matrix)
        reshaped['index']=list(range(0,reshaped.shape[0])) # this is the current ordering of rows of newX
        reshaped=reshaped.sort_values(by='snpBP',ascending=True) # need to order by BP so I know order of SNPs in ldmat from PLINK 
        newX=newX[reshaped['index'].values];mstart=newX.shape[0] # Need to make sure the data I already created (newX/bX) matches this order
        minBP=min(reshaped['snpBP'].values); maxBP=max(reshaped['snpBP'].values)
        thingsToMonitor['BPSpanOfIVSet']=maxBP/1e6-minBP/1e6
        
        if ldOtherLociWindow>0:
            # defaults: ldOtherLociOtherPt=5e-5,ldOtherLociR2=0.1
            mask1=merged[(merged['snpBP']<minBP) & (merged['snpBP']>(minBP-ldOtherLociWindow*1e6))]
            mask2=merged[(merged['snpBP']>maxBP) & (merged['snpBP']<(maxBP+ldOtherLociWindow*1e6))]
            mask1=mask1[abs(mask1['geneZ'])>norm.ppf(1-ldOtherLociOtherPt)]
            mask2=mask2[abs(mask2['geneZ'])>norm.ppf(1-ldOtherLociOtherPt)]
            mask1=mask1[['geneSNP','snpBP']].drop_duplicates(keep='first').sort_values(by='snpBP',ascending=True)
            mask2=mask2[['geneSNP','snpBP']].drop_duplicates(keep='first').sort_values(by='snpBP',ascending=True)
            # write out big list in this order: mask1, reshaped, mask2
            big_list=pandas.concat([mask1['geneSNP'],reshaped['geneSNP'],mask2['geneSNP']])
            big_list.to_csv(outDir1, header=False, index=False, sep=" ") # writing out list of SNPs to give PLINK to make LD calculations
            cmd=[callPlink(), "--r", "square", "--bfile", ldRefDir, "--extract", outDir1, "--out", outDir2] # the command for PLINK
            out=subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) # executing the command in the terminal
            ldmat=numpy.genfromtxt(outDir2+'.ld', dtype='f') # reading LD matrix back in
            # now need to know if any SNPs in reshaped are in LD>ldOtherLociR2(=0.1 be default) with SNPs outside of reshaped
            maskreshaped=numpy.ones((ldmat.shape[0],),dtype='bool')
            maskreshaped[:mask1.shape[0]]=False; maskreshaped[(mask1.shape[0]+reshaped.shape[0]):]=False
            inds1=ldmat[mask1.shape[0]:(mask1.shape[0]+reshaped.shape[0]),:][:,:mask1.shape[0]]; inds1=numpy.nan_to_num(inds1,nan=0)
            inds2=ldmat[mask1.shape[0]:(mask1.shape[0]+reshaped.shape[0]),:][:,(mask1.shape[0]+reshaped.shape[0]):]; inds2=numpy.nan_to_num(inds2,nan=0)
            ldmask1=(((inds1**2)>ldOtherLociR2).sum(axis=1))>0
            ldmask2=(((inds2**2)>ldOtherLociR2).sum(axis=1))>0
            ldmask=(ldmask1+ldmask2)==0
            reshaped=reshaped[numpy.array(ldmask).astype(bool)]
            newX=newX[ldmask]
            ldmat=ldmat[maskreshaped,:][:,maskreshaped] # subset to only SNPs in reshaped
            ldmat=ldmat[ldmask,:][:,ldmask] # subset again to only SNPs in reshaped uncorrelated with SNPs in surrounding window
            # some LD correlations are totally missing for a SNP (rare, but possible). remove these. (possible why? observed AF=0 or 1 so no variation?)
            mask=((numpy.isnan(ldmat).sum(axis=0))==ldmat.shape[0])==False # create mask: True if non-missing LD; False if all missing LD (I've never seen partial missing)
            reshaped=reshaped[mask]; newX=newX[mask,:]; ldmat=ldmat[mask,:][:,mask] # drop SNPs with missing LD
            del mask1,mask2,ldmask,inds1,inds2,ldmask1,ldmask2
        else:
            # personal check (should be all True's): [sum(newX[:,xx]==reshaped[cn[xx]].values)==newX.shape[0] for xx in range(0,newX.shape[1])]
            reshaped['geneSNP'].to_csv(outDir1, header=False, index=False, sep=" ") # writing out list of SNPs to give PLINK to make LD calculations
            cmd=[callPlink(), "--r", "square", "--bfile", ldRefDir, "--extract", outDir1, "--out", outDir2] # the command for PLINK
            out=subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) # executing the command in the terminal
            ldmat=numpy.genfromtxt(outDir2+'.ld', dtype='f') # reading LD matrix back in
            # ldmat is ordered by BP position, NOT by order of SNPs in myExtract.txt. That is why I sorted reshaped by snpBP
            # some LD correlations are totally missing for a SNP (rare, but possible). remove these. (possible why? observed AF=0 or 1 so no variation?)
            mask=((numpy.isnan(ldmat).sum(axis=0))==ldmat.shape[0])==False # create mask: True if non-missing LD; False if all missing LD (I've never seen partial missing)
            reshaped=reshaped[mask]; newX=newX[mask,:]; ldmat=ldmat[mask,:][:,mask] # drop SNPs with missing LD
        
        thingsToMonitor['snpsDroppedBCCHPCorrectionByLD']=mstart-newX.shape[0]
        ldc=abs(ldmat.copy()); ss=(ldc>ldUpperLimit).sum(axis=0).sum(axis=0)
        if (ldmat.shape[0]<nMinIVs) or (ss==ldmat.shape[0]**2): # if not enough SNPs or all are in near-perfect LD with each other
            continue
        thingsToMonitor['snpsDroppedBCMissingLD']=sum(mask==False) # save what happened here
        snps=reshaped['geneSNP'].values # find SNPs currently considered
        someDropped,snpsKept=LDShrinkUpper(ldmat,ldUpperLimit,snps); ld0=someDropped.copy() # find those not in ~perfect LD (>tau) with any others (makes ldmat pos def)
        mask=reshaped['geneSNP'].isin(snpsKept)
        reshaped=reshaped[mask] # restrict to only SNPs from line above    
        newX=newX[mask]
        thingsToMonitor['snpsDroppedBCHighLD']=ldmat.shape[0]-ld0.shape[0] # number of SNPs dropped because of high LD correlations
        ######################### if the number of IVs is very small, we can't really trust any results
        if ld0.shape[0]<nMinIVs: 
            continue
        #########################
        # R=ld0.copy();kTaper=int(ld0.shape[0]/4);kBand=ld0.shape[0]+1;lambda_=0.01;epsilon=1e-4;toCorr=True
        ld0,alpha,deco=regularizeLD(ld0.copy(),int(ld0.shape[0]/4),ld0.shape[0]+1,0.1,epsilon=1e-4,toCorr=True,justMCP=False) # int rounds down
        w,v=deco[0],deco[1]
        sparsePrecision=v@numpy.diag(1/w)@v.T
        sparsePrecisionSqrt=v@numpy.diag(1/w**0.5)@v.T
        thingsToMonitor['alphaToMakeLDMatPosDef']=alpha
        ############################################################### end middle 1/3rd (correlations and imputation)
        
        ############################################################### start final 1/3rd (MR and saving data) (avg ~2 seconds per gene group)
        # organizing data
        bx=numpy.array(newX); by=reshaped['phenoZ'].values; bx0=bx.copy(); by0=by.copy()
        p=1 if len(bx.shape)==1 else bx.shape[1]
        q=1 if len(by.shape)==1 else by.shape[1]
        m=bx.shape[0]
        bx=bx.reshape((m,p)); by=by.reshape((m,q))
        bx0=bx0.reshape((m,p)); by0=by0.reshape((m,q))
        thingsToMonitor['finalNumIVs']=bx.shape[0]
        thingsToMonitor['outcomeEstH2InLocus']=(by.T@by)-by.shape[0]
        ########################################## 
        # remove LD structure at the outset - seems to cause problems in standard error estimation - I currently have not found the problem or another solution
        bx=sparsePrecisionSqrt@bx; by=sparsePrecisionSqrt@by
        # if isGtex: # if GTEx, all effect estimates are already joint
        #     bx=bx0.copy()
        #     by=by0.copy()
        #     ld0=numpy.eye(m)
        #     sparsePrecision=ld0.copy()
        #     sparsePrecisionSqrt=ld0.copy()
        ld0og=ld0.copy(); sparsePrecisionog=sparsePrecision.copy()
        ld0=numpy.eye(bx.shape[0]); sparsePrecision=ld0.copy(); sparsePrecisionSqrt=ld0.copy() 
        # this has implications for how nclumps are estimated - that's why I keep a copy of the original LD matrix
        ##########################################
        ########################################## (yes I am doing this for a second time - works well in both cases)
        # remove genes with norm(newX[:,j],2)^2<5m
        toKeep=(numpy.diag(bx.T@bx)>(hessMinScale*bx.shape[0]/2)) # can change factor if you want; default is hessMinScale=5
        if sum(toKeep)==0:
            continue
        
        bx=bx[:,toKeep]; bx0=bx0[:,toKeep]
        ccComplete=ccComplete[toKeep,:][:,toKeep]
        droppedGenes=numpy.array(cn)[toKeep==False].tolist()
        cn=numpy.array(cn)[toKeep].tolist()
        reshaped=reshaped.drop(labels=droppedGenes,axis=1)
        p=sum(toKeep)
        ##########################################
        # ells=ldscore(ld0); ells=numpy.column_stack(([1]*ells.shape[0],ells)) # only approximate LD scores (only considers LD in ld0)
        # (numpy.linalg.inv(ells.T@ells)@ells.T@bx).shape
        # variance-covariance matrices of measurement errors
        VV=numpy.array(1).reshape((1,1))
        UU=ccComplete.copy() # maybe not memory friendly to make a copy but both object names work well in context
        # now to calculate UV
        if assumeNoSampleOverlap: # if user wants to assume there is no sample overlap
            UV=numpy.zeros((p,1))
        else:
            bybx=numpy.hstack((by,bx)) # already no missing because I removed it from newX (bx)
            mask=abs(bybx>q0Correls).sum(axis=1) == 0 # True values in mask are the SNPs to use in correls calculations
            mask=numpy.where(mask)[0] # for some reason returns a tuple so I add the [0]
            if mask.shape[0]<nMinCorrels:
                uv0=numpy.zeros((p,q))*float('nan')
                a_=numpy.hstack((VV,uv0.T))
                b_=numpy.hstack((uv0,UU))
                c_=numpy.vstack((a_,b_))
                md=MaxDet(c_,doToBadSolutions='drop')[0]
                mask=abs(md)>1
                md[mask]=0
                UV=md.copy(); UV=UV[:,0]; UV=UV[1:]; UV=UV.reshape((p,1))
            else:
                UV=numpy.corrcoef(sparsePrecisionSqrt[mask,:][:,mask]@bybx[mask,:],rowvar=False)[:,0] # Corr(by,bx) will be in first column
                UV=numpy.delete(UV,0); UV=UV.reshape((p,1)) # drop the 1 for Corr(by,by)
        
        # add intercept-relevant terms for MRJones
        bx_=numpy.column_stack(([1]*m,bx)); bx0_=numpy.column_stack(([1]*m,bx0)); 
        UU_=numpy.column_stack(([0]*p,UU)); UU_=numpy.row_stack(([0]*(p+1),UU_))
        UV_=numpy.row_stack((0,UV))
        ### shrinking of UU and UV to make sure they don't introduce big problems in causal estimation
        if shrinkBiasCorrection:
            UU=UU*0.5; UU_=UU_*0.5
            UV=UV*0.5; UV_=UV_*0.5
        
        if (opAlpha=="dynamic") & (bx.shape[0]>10):
            c_c=numpy.corrcoef(bx,rowvar=False); numpy.fill_diagonal(c_c.reshape((p,p)), 0)
            opAlphaVal=float(numpy.max(abs(c_c)))
        elif is_number(opAlpha):
            opAlphaVal=opAlpha
        else:
            opAlphaVal=0
        
        thingsToMonitor['MRJonesAlpha']=opAlphaVal
        ## first defining hyperparameter lambda1, lambda2/tau grids
        # multiple lambda vectors
        lamvec=numpy.linspace(0.01,1,40)
        t1=numpy.linspace(0.01,1/2,10)
        t2=numpy.linspace(1,3,10)
        t3=numpy.linspace(5,10,10)
        t4=numpy.linspace(20,40,20)
        tauvec=numpy.concatenate((t1,t2,t3,t4))
        n=reshaped['geneN'].values; n=numpy.nanmedian(n)
        ### if you want to fix alpha at a certain value to make it faster
        out=MRJones(by0,bx0_,ld0og,UU_,UV_,lamvec=lamvec,tauvec=tauvec,rho_theta=3/2)
        # out=MRJones(by,bx,ld0,UU,UV,lamvec=lamvec,tauvec=tauvec)
        deltas=out['gamma']
        finalEsts=out['theta']; 
        # variance explained
        if all(finalEsts==0):
            r2mr=0
        else:
            em=numpy.column_stack((by0,bx0_@finalEsts)); 
            if numpy.all(em[:,0]==em[0,0]) | numpy.all(em[:,1]==em[1,1]):
                r2mr=0
            else:
                r2mr=numpy.corrcoef(em,rowvar=False)[0,1]**2
        estsVars=out['covg']; estsVars=estsVars.reshape((estsVars.shape[0],estsVars.shape[0]))
        # The MRJones() function already drops the intercept-relevant term from the variance-covariance matrix
        finalEsts=finalEsts[1:]
        estDf=pandas.DataFrame.from_dict({'gg': numpy.array(cn), 'Est': finalEsts})
        thingsToMonitor['nNonzeroMRJonesDeltas']=sum(deltas!=0)
        # if all(deltas!=0):
        #     continue
        # intercept terms not included in new MRJones function
        # isNotZero=(finalEsts!=0) # finding which genes we need a variance estimate for (ie non-shrunken-to-0 ones)
        # don't even consider MR-Jones SEs - we rely on post selection later anyway
        paraDf=estDf.copy()
        paraDf['SE']=math.inf
        # SEs=(numpy.diag(estsVars)**0.5).squeeze()
        # # if all estimates shrunken to 0
        # if len(SEs)==0:
        #     seDf=estDf.copy(); seDf=seDf.rename(columns={'Est': 'SE'}); seDf['SE']=math.inf
        # elif all(finalEsts==0): # MRJones will return all SEs as 0s for all genes, evenif coefficients are 0
        #     seDf=pandas.DataFrame.from_dict({'gg': numpy.array(cn), 'SE': numpy.ones((len(cn),))*float('nan')})
        # else:
        #     seDf=pandas.DataFrame.from_dict({'gg': numpy.array(cn)[isNotZero], 'SE': SEs})
        
        # paraDf=pandas.merge(estDf,seDf,how='outer',left_on='gg', right_on='gg')
        # paraDf['SE']=paraDf['SE'].fillna(math.inf)
        # add lead SNPs for each gene to use in later plotting
        leadGeneSNPDf['Gene']=leadGeneSNPDf['Gene']+'_Z' # for merging
        paraDf=pandas.merge(paraDf,leadGeneSNPDf,left_on='gg',right_on='Gene') # merge
        Ps=[2*norm.cdf(-abs(paraDf['Est'].values[x]/paraDf['SE'].values[x]),0,1) if paraDf['SE'].values[x]>0 else 1 for x in range(0,paraDf.shape[0])]
        # conditional F-statistics (MRJones recommended for this; lamvec determined automatically from within the function)
        if bx.shape[1]<2:
            condFs=float('nan')
        else:
            condFs=conditionalFMRBEE(bx,UU,sparsePrecisionSqrt,tauLowerQuant=0.75,tauUpperQuant=0.95,opAlphaVal=opAlphaVal,method='mrbee') # 'notfree': MRBEE; 'free': MR-Jones
        paraDf['conditionalF']=condFs
        # adding gene BP to results
        gl=[]; cn_=[cn[x].split('_',1)[0] for x in range(0,len(cn))] # cn_ drops '_Z' from cn (all cn==paraDf['gg'])
        for _ in range(0,len(cn_)):
            mask=cn_[_]==geneLocs['Gene'].values
            ta=geneLocs[mask]['geneBP'].values
            gl.append(float(ta))
        
        # finally, perform analysis using MVMR and uni MR with MR-Egger (this is equivalent to SMR and the Porcu et al methods)
        ivwout=IVW(bx,by,sparsePrecision,newX,uniIVPThreshold=1e-5) # only IVs with P<uniIVPThreshold will be considered in univariable MR
        ivwout['gg']=cn # cn is the order of genes in paraDf (check: [cn[x]==paraDf['gg'].values[x] for x in range(0,len(cn))])
        paraDf=pandas.merge(paraDf,ivwout,left_on='gg',right_on='gg')
        # add BP position for lead SNP to paraDf
        paraDf=pandas.merge(paraDf,genedf[['geneSNP','snpBP']].drop_duplicates(keep='first'),left_on='geneSNP',right_on='geneSNP') # add snpBP to paraDf for lead SNP
        paraDf['Gene']=[paraDf['Gene'][_].split('_',1)[0] for _ in range(0,paraDf.shape[0])] # for merging
        paraDf=pandas.merge(paraDf,genedf[['Gene','geneBP']].drop_duplicates(keep='first'), left_on='Gene',right_on='Gene') # add gene BP to paraDf for each gene
        paraDf=paraDf.sort_values('geneBP') # order by gene BP
        # approximate bias in using SMR
        smrbias=perGeneSMRBias(bx,by,ld0,sparsePrecision,newX,UU,UV,cn_,IVPThresh=UniMRIVPThreshold,r2Thresh=0.25,minNIVs=5)
        paraDf=pandas.merge(paraDf,smrbias,left_on='Gene',right_on='gene')
        # MRBEE with ridge
        mrbeeRidge=MRBEERidge(bx,by,UU,UV,VV,sparsePrecision,nLambda=20)
        ridgeDf=pandas.DataFrame({'Gene': cn_, 'MRBEERidgeEst': mrbeeRidge[0].squeeze(), 'MRBEERidgeSE': (numpy.diag(mrbeeRidge[1])**0.5).squeeze()})
        paraDf=pandas.merge(paraDf,ridgeDf,left_on='Gene',right_on='Gene')
        # MRBEE with Huber loss
        gam=sum((by.squeeze()-numpy.median(by.squeeze()))**2)/by.shape[0]/2
        huberE,huberV,huberW=MRBEEHuber(bx,by,ld0,UU,UV,gamma=gam,initial="bee",eps=0.0001*p,max_iter=20,boot=False) # keep gamma small (as gamma->inf, huber->quantile)
        nim=numpy.ones((p+1,),dtype='bool'); nim[0]=False # non-intercept mask
        huberdf=pandas.DataFrame({'Gene': cn_,'MRBEEHuberEst': huberE.squeeze()[nim], 'MRBEEHuberSE': ((numpy.diag(huberV)**0.5).squeeze())[nim]})
        paraDf=pandas.merge(paraDf,huberdf,left_on='Gene',right_on='Gene')
        ### MRBEE Post-Selection
        mask=finalEsts!=0 # True where gene should be included, False where not
        if sum(mask)==0:
            toa=numpy.array([0]*len(mask))*float('nan') # vector of NaN's
            toAdd=pandas.DataFrame({'Gene': cn_, 'MRBEEPosSelEst': toa,'MRBEEPosSelSE': toa})
            paraDf=pandas.merge(paraDf,toAdd,left_on='Gene',right_on='Gene')
        else:
            bxkeep=bx[:,mask]; UUkeep=UU[mask,:][:,mask]; UVkeep=UV[mask]
            p1,p2,p3,p4,p5=imrbee(bxkeep,by,UUkeep,UVkeep,VV,sparsePrecision,ld0,0.05/bxkeep.shape[0]**0.5,boot=False)
            p1=p1[1:]; p2=p2[1:,:][:,1:] # remove intercept-relevant terms 
            toAdd=pandas.DataFrame({'Gene': numpy.array(cn_)[mask].tolist(), 'MRBEEPosSelEst': p1.squeeze(), 'MRBEEPosSelSE': numpy.diag(p2)**0.5})
            paraDf=pandas.merge(paraDf,toAdd,how="outer",left_on='Gene',right_on='Gene')
        # single SNP
        singlesnpdf=singleSNP(bx,by,cn_)
        paraDf=pandas.merge(paraDf,singlesnpdf,how='outer',left_on='Gene',right_on='Gene')
        # estimate gene expression heritabilities and numbers of SNPs with nonzero associations with gene expression (from theory)
        # eps=[sum(abs(bx[:,__])>2) for __ in range(0,bx.shape[1])] # estimates of beta mixture proportions (see next line) (proportions of nonzero true associations with gene exp)
        # eps=numpy.array(eps)/bx.shape[0]
        # medn=numpy.median(reshaped['geneN'].values) # median sample size for these genes in this tissue
        ### this may not be the true eps ... recall I assumed the true effects are not in LD
        # estimate heritability
        # ... how to do it?
        # from theory of the mixture (see Supplement)
        # from LDSC (approximate ldscores)
        # ldsc=ldscore(ld0og); ldsc=sparsePrecisionog@ldsc; ldsc=numpy.column_stack(([1]*ldsc.shape[0],ldsc))
        # hat=(numpy.linalg.inv(ldsc.T@ldsc)@ldsc.T@sparsePrecision@bx)
        # hat=hat[1,:]
        # problems with ldsc (estimates outside range {-1,1}
        # basically performing clumping without window size restriction within the SNP set now: for each gene, number of P(Z)<0.05 and independent with other SNPs
        nclumps=[]; r2=0.01 # like r2 in PLINK
        ld0=ld0og.copy() # since all MR calculations are done, resetting LD matrix to be actual LD matrix (ie not identity matrix bc I made a transformation)
        for p_ in range(0,bx.shape[1]):
            mask=numpy.abs(bx0[:,p_])>norm.ppf(1-1e-5) # (1e-5)th quantile; only considered LD between P(Z)<0.05 SNPs
            ld0mask=ld0[mask,:][:,mask]**2 # squaring all values for convenience
            ninds=[]
            for s_ in range(0,ld0mask.shape[0]):
                m1=numpy.zeros((ld0mask.shape[0],),dtype='bool'); m1[s_]=True # creating mask to index this specific SNP
                ldm1=ld0mask[m1,:][:,m1==False].squeeze();
                ss=numpy.sum(ldm1>(r2)) # SNPs with squared LD < r2 will be considered independent
                if ss==0: # if no other SNPs in LD with this one
                    ninds.append(1) # put 1 (I will sum up later)
            if len(ninds)==0:
                nclumps.append(1) # if can't find any for this gene, just put 1 - assume there is at least 1
            else:
                nclumps.append(sum(ninds)) # else, put total number I counted
        # nclumps is an estimate of the number of independent causal SNPs in this locus
        nclumps=numpy.array(nclumps)
        eps=numpy.array(nclumps)/bx.shape[0] # an estimate of the proportion of causal SNPs in this locus
        medns=reshaped['geneN'].values; N=numpy.diag(1/medns);
        medns=numpy.nanmedian(medns)
        cisEstNaive=((numpy.diag(bx.T@bx)-1/medns)/bx.shape[0]-1)*nclumps/medns # medns is the median sample sizes for each gene
        toAdd=pandas.DataFrame({'Gene': cn_, 'nClumps': nclumps, 'h2CisEstNaive': cisEstNaive})
        paraDf=pandas.merge(paraDf,toAdd,left_on='Gene',right_on='Gene')
        paraDf['MRJonesR2']=r2mr
        # pretty sure h2 cannot be reliably estimated from the data because bx=bhat*nk/sigk where I do not know sigk and there is evidence that it is not 1.        
        # d1['data']=reshaped # probably don't want to save the data ... will take up a lot of memory
        d1={}
        d1['data']={'Gene': cn_, 'GeneBP': gl, 'IVs': reshaped['geneSNP'].values.tolist()}
        if saveData: # currently not used after this point. data gets outputted by MVMRWorkhorse() but is not saved after that
            d1['bX']=bx; d1['by']=by; d1['UU']=UU; d1['UV']=UV; d1['regLD']=ld0
        paraDf.index=paraDf['Gene']; paraDf=paraDf.reindex(cn_)
        d1['meta']=paraDf # put whole data frame
        # outerDict[ggKeys[ogene]]=d1
        outerDict['group'+str(ogene)]=d1
        thingsMonitored[ggKeys[ogene]]=thingsToMonitor # the key is the key for the group in geneGroups
        ### MN networks
        try:
        # back-tabbed the if: else: below twice
          if (r2mr>networkR2Thres) & (bx.shape[0]>25) & (bx.shape[1]>1):
              XX=numpy.column_stack((by,bx))
              XXSE=numpy.ones(XX.shape)
              Rnoise=numpy.eye(bx.shape[1]+1); Rnoise[1:,1:]=UU; Rnoise=posDefifyCorrMat(Rnoise)[0]
              lv1=lamvec=numpy.linspace(1,5,10)/100; lv2=lamvec=numpy.linspace(7,12,10)/100; lv3=lamvec=numpy.linspace(15,20,10)/100; lv=numpy.concatenate((lv1,lv2,lv3))
              net=entropy_mcp_spearman_sampling(XX,XXSE,Rnoise/5,lamvec=lv) # make Rnoise small
              c1=sum(sum(numpy.isnan(net['Theta'])))>0
              c2=sum(sum(numpy.isnan(net['Theta_Alt'])))>0
              if (c1==True) & (c2==False):
                  edges=(net['Theta_Alt']!=0).astype(int)
              elif (c1==False) & (c2==True):
                  edges=(net['Theta']!=0).astype(int)
              elif (c1==False) & (c2==False):
                  edges=(net['Theta']!=0).astype(int)
              else: # could not conpute network
                  edges=numpy.zeros(XX.shape)
              numpy.fill_diagonal(edges,0)
          else:
              edges=numpy.zeros((bx.shape[1]+1,bx.shape[1]+1))
        except:
            edges=numpy.zeros((bx.shape[1]+1,bx.shape[1]+1))
        edgeDict['group'+str(ogene)]=edges 
        ############################################################### final 1/3rd (MR and saving data)
    return thingsMonitored,outerDict,edgeDict

def UVMRworkhorse(merged,writableDir,ldRefDir,LDWindowMb=1,LDMaxInWindow=0.5,eQTLInOtherWindowP=5e-5,IVP=0.025,minNumIVs=10,LDWithinIVSetMax=0.9,
                 restrictToOutcomeSignals=True,outcomeSignalsWithinXMb=1):
    # LDWindowMb: remove IVs that have LD>LDMaxInWindow with other eQTLs within LDWindowMb
    # LDMaxInWindow: IVs with LD>LDMaxInWindow with other eQTLs within LDWindowMb are removed
    # eQTLInOtherWindowP: from two above, only consider SNPs in nearby windows to be eQTLs for other genes if P<eQTLInOtherWindowP
    # IVP: P-value threshold for defining IVs in UniMR
    # minNumIVs: self-explanatory
    genes=merged['Gene'].unique()
    outcomeClumpBPs=merged[merged['isOutcomeClump']==True]['snpBP'].values
    mergedsig=merged[abs(merged['geneZ'])>norm.ppf(1-eQTLInOtherWindowP)]
    thingsToMonitor={}; outerDict={}
    for ogene in range(0,len(genes)):
        mcut=merged.copy()[merged['Gene']==genes[ogene]]
        mcut=mcut[abs(mcut['geneZ'])>norm.ppf(1-IVP)] 
        if mcut.shape[0]<minNumIVs: # if not enough eQTLs, next gene
            continue
        # if user only wants to do the analysis near outcome signals, make sure there are some here
        if restrictToOutcomeSignals==True:
            diffs=[min(abs((mcut['snpBP'].values)[x]-outcomeClumpBPs)) for x in range(0,mcut.shape[0])]
            if any(numpy.array(diffs)<(1e6*outcomeSignalsWithinXMb))==False:
                continue
        # find LD between these SNPs
        mcut=mcut.drop_duplicates(subset=['geneSNP'],keep=False) # drop duplicate SNPs (cannot appear in file given to PLINK)
        mcut=mcut.sort_values(by='snpBP',ascending=True)
        if mcut.shape[0]<minNumIVs: # if no non-duplicated eQTLs, next gene
            continue
        # make sure not in LD with eQTLs with any genes within LDWindowMb
        minBP,maxBP=min(mcut['snpBP']),max(mcut['snpBP'])
        mask1=mergedsig[(mergedsig['snpBP']<minBP) & (mergedsig['snpBP']>(minBP-LDWindowMb*1e6)) & (mergedsig['Gene']!=genes[ogene])]
        mask2=mergedsig[(mergedsig['snpBP']>maxBP) & (mergedsig['snpBP']<(maxBP+LDWindowMb*1e6)) & (mergedsig['Gene']!=genes[ogene])]
        mask1=mask1[['geneSNP','snpBP']].drop_duplicates(keep='first').sort_values(by='snpBP',ascending=True)
        mask2=mask2[['geneSNP','snpBP']].drop_duplicates(keep='first').sort_values(by='snpBP',ascending=True)
        # write out big list in this order: mask1, reshaped, mask2
        big_list=pandas.concat([mask1['geneSNP'],mcut['geneSNP'],mask2['geneSNP']])
        outDir1=writableDir+"/myExtract.txt"
        outDir2=writableDir+"/ld"
        big_list.to_csv(outDir1, header=False, index=False, sep=" ") # writing out list of SNPs to give PLINK to make LD calculations
        cmd=[callPlink(), "--r", "square", "--bfile", ldRefDir, "--extract", outDir1, "--out", outDir2] # the command for PLINK
        out=subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) # executing the command in the terminal
        ldmat=numpy.loadtxt(outDir2+'.ld', dtype='f') # reading LD matrix back in
        # now need to know if any SNPs in reshaped are in LD>ldOtherLociR2(=0.1 be default) with SNPs outside of reshaped
        maskreshaped=numpy.ones((ldmat.shape[0],),dtype='bool')
        maskreshaped[:mask1.shape[0]]=False; maskreshaped[(mask1.shape[0]+mcut.shape[0]):]=False
        inds1=ldmat[mask1.shape[0]:(mask1.shape[0]+mcut.shape[0]),:][:,:mask1.shape[0]]; inds1=numpy.nan_to_num(inds1,nan=0)
        inds2=ldmat[mask1.shape[0]:(mask1.shape[0]+mcut.shape[0]),:][:,(mask1.shape[0]+mcut.shape[0]):]; inds2=numpy.nan_to_num(inds2,nan=0)
        ldmask1=(((inds1**2)>(LDMaxInWindow**2)).sum(axis=1))>0
        ldmask2=(((inds2**2)>(LDMaxInWindow**2)).sum(axis=1))>0
        ldmask=(ldmask1+ldmask2)==0
        if sum(ldmask)<minNumIVs: # if no SNPs that are valid IVs for univariable MR
            continue
        mcut=mcut[ldmask]
        ldmat=ldmat[maskreshaped,:][:,maskreshaped] # subset to only SNPs in mcut
        ldmat=ldmat[ldmask,:][:,ldmask] # subset again to only SNPs in reshaped uncorrelated with SNPs in surrounding window
        # some LD correlations are totally missing for a SNP (rare, but possible). remove these. (possible why? observed AF=0 or 1 so no variation?)
        mask=((numpy.isnan(ldmat).sum(axis=0))==ldmat.shape[0])==False # create mask: True if non-missing LD; False if all missing LD (I've never seen partial missing)
        mcut=mcut[mask];
        del mask1,mask2,ldmask,inds1,inds2,ldmask1,ldmask2
        someDropped,snpsKept=LDShrinkUpper(ldmat,LDWithinIVSetMax,mcut['geneSNP'].values); ld0=someDropped.copy() # find those not in LD too highly
        mcut=mcut[mcut['geneSNP'].isin(snpsKept)]
        # make sure LD between remaining SNPs isn't too high to cause inflation
        if ld0.shape[0]<minNumIVs:
            continue
        # perform MR (bias-correction for SigmaUV assumed to be 0 since very unlikely gene expression samples overlap with any phenotype samples)
        # first step is to transform by sqrt of inverse of LD matrix
        bxog,byog=mcut['geneZ'].values,mcut['phenoZ'].values
        m=bxog.shape[0]
        bxog,byog=bxog.reshape((m,1)), byog.reshape((m,1))
        # but first make LD matrix positive definite
        ld0,alpha,deco=regularizeLD(ld0,int(ld0.shape[0]/4),ld0.shape[0]+1,0.1,epsilon=1e-4,toCorr=True) # int() rounds down
        ld0inv=deco[1]@numpy.diag(1/deco[0])@deco[1].T
        ld0invsqrt=deco[1]@numpy.diag(1/deco[0]**0.5)@deco[1].T
        bx=ld0invsqrt@bxog; by=ld0invsqrt@byog; ld0=numpy.eye(m); ld0inv=ld0.copy(); ld0invsqrt=ld0.copy()
        bx_=numpy.column_stack(([1]*m,bx))
        UU=numpy.ones((1,1)) # under transformation, variance of measurement error is 1
        UV=numpy.zeros((1,1))
        VV=numpy.ones((1,1))
        UU_=numpy.zeros((2,2)); UU_[1,1]=1
        UV_=numpy.zeros((2,1))
        # MRBEE
        est0,V0,Outliers,kR,k=imrbee(bx,by,UU,UV,VV,ld0inv,ld0,PleioPThreshold=0.05/bx.shape[0]**0.5,boot=False,niter=1000,initial="bee",max_iter=15)
        pleiop=SpleioP(bx_,by,UU_,UV_,VV,est0,V0)
        # IVW (MR-Egger actually)
        bxx=numpy.linalg.inv(bx_.T@bx_)
        ivwest=bxx@bx_.T@by
        s2=sum((by.squeeze()-bx_@ivwest.squeeze())**2)/(m-2)
        ivwvar=s2*bxx
        # huber loss
        huberE,huberV,huberW=MRBEEHuber(bx,by,ld0,UU,UV,gamma=0.5,initial="bee",eps=0.001,max_iter=15,boot=False) # keep gamma small (as gamma-> huber->quantile)
        # end
        paraDict={'Gene': genes[ogene], 'GeneBP': int(mcut['geneBP'].unique()),
                  'MRBEEEst': float(est0[1]), 'MRBEESE': float(V0[1,1]**0.5), 
                  'HuberEst': float(huberE[1]), 'HuberSE': float(huberV[1,1]**0.5),
                  'IVWEst': float(ivwest[1]), 'IVWSE': float(ivwvar[1,1]**0.5), 'nIVs': m}
        outerDict[genes[ogene]]=pandas.DataFrame(paraDict, index=[0])
    return outerDict

def huberWeight(e,gamma):
    e_=numpy.array(e)
    mask=abs(e_)<=gamma
    w=e_*0
    w[mask]=1
    w[mask==False]=gamma/abs(e[mask==False])
    return w

def bootVarWLS(x,y,w,UU,k=1000):
    m=len(y);p=x.shape[1]
    ests=numpy.zeros((k,p))
    for ii in range(0,k):
        ind=numpy.random.randint(0,m,m)
        xi=x[ind];yi=y[ind];wi=w[ind];D=numpy.diag(wi.squeeze())
        ei=numpy.linalg.inv(xi.T@D@xi-sum(wi)*UU)@xi.T@D@yi
        ests[ii,:]=ei.squeeze()
    return numpy.cov(ests,rowvar=False)

def MRBEEHuber(bx,by,R,UU,UV,gamma=0.5,initial="ivw",eps=0.001,max_iter=15,boot=False):
    m=bx.shape[0];p=bx.shape[1]
    UU_=numpy.column_stack(([0]*p,UU)); UU_=numpy.row_stack(([0]*(p+1),UU_))
    UV_=numpy.row_stack((0,UV))
    w,v=numpy.linalg.eig(R); Rinvsqrt=v@numpy.diag(1/w**0.5)@v.T
    bx_=numpy.column_stack(([1]*m,Rinvsqrt@bx)) # transforming up front
    by_=Rinvsqrt@by # transforming up front
    if initial=="robust":
        mod=QuantReg(by_,bx_) # includes an intercept
        res=mod.fit(q=0.5)
        est0=res.params.reshape((p+1,1))
        res0=by_-bx_@est0
    else:
        est0=numpy.linalg.inv(bx_.T@bx_)@bx_.T@by_
        res0=by_-bx_@est0
    error=1+eps; k=0
    while (error>eps) & (k<max_iter):
        k=k+1
        w0=huberWeight(res0,gamma)
        D=numpy.diag(w0.squeeze())
        adjfactor=sum(w0)
        estk=numpy.linalg.inv(bx_.T@D@bx_-adjfactor*UU_)@bx_.T@D@by_
        error=(est0-estk).T@(est0-estk)
        est0=estk.copy()
        res0=by_-bx_@est0
    adjfactor=sum(w0)
    if boot:
        varest=bootVarWLS(bx_,by_,w0,UU)
    else:
        s2=sum((w0*res0)**2)/(adjfactor-2) # this may underestimate variance (assumes fixed weights)
        Hinv=numpy.linalg.inv(bx_.T@D@bx_-adjfactor*UU_)
        XWX=bx_.T@(D**2)@bx_
        varest=s2*Hinv@XWX.T@Hinv
    return est0,varest,w0

# calculating inflation using areas where I know there is no outcome signal. need to use the actual data for this
# the data I will use for this is in the `merged` object and I'll want to do this before running the analyses
def calcInflationFactor(merged,ldRefDir,writableDir,mbWindow=5,ldUpperLimit=0.5,outcomeP=0.01,nGenesToTry=50,geneDistMb=1,minM=25):
    # 1) find loci in which there is no outcome signal
    # 2) select a gene with no eQTLs
    # 3) calculate LD for SNPs available for gene
    # 4) prune LD to only |r|<ldUpperLimit (default=0.5)
    # 5) calculate test statistic and store
    # 6) repeat for many genes that are at least geneDistMb Mb away from each other (default is 5)
    # 7) store inflation
    outDir1=writableDir+'/myExtract.txt' # (for saving SNPs for LD calculations) write out to a directory that must be writable (since the user is there)
    outDir2=writableDir+'/tempOut' # (for saving LD matrix)
    # 1)
    merged_=merged.copy().reset_index()
    clumpBPs=merged_[merged_['isOutcomeClump']==True]['snpBP'].unique()
    mask=numpy.ones((merged_.shape[0]),dtype='bool')
    for _ in range(0,len(clumpBPs)):
        m_=merged[abs(merged['snpBP']-clumpBPs[_])<(mbWindow*1000000)] # the 1Mb region around clumpBPs
        mask[m_.index]=False
    merged_=merged_.sort_values(by='snpBP',ascending=True)
    merged_=merged_[mask]
    # 2)
    # first finding genes that are at least geneDistMb Mb away from each other
    genedf=merged_[['Gene','geneBP']].drop_duplicates().sort_values(by='geneBP',ascending=True)
    if genedf.shape[0]==0:
        return 1,0,1,[],0,float('nan')
    mask=numpy.zeros((genedf.shape[0],),dtype='bool'); mask[0]=True
    for _ in range(1,len(mask)): # starting with second gene
        tind=numpy.where(mask==True)[0][-1]
        toadd=True if abs(genedf.iloc[_,1]-genedf.iloc[tind,1])>(geneDistMb*1e6) else False
        mask[_]=toadd
        if sum(mask)==nGenesToTry: # if I already have enough, don't need to continue
            break
    
    ggs=genedf[mask]['Gene'].unique()
    cstats=list(); ccs=list(); ccps=list(); counter=0; usedGenes=list(); ms=list()
    for _ in range(0,len(ggs)):
        dat=merged_[merged_['Gene']==ggs[_]]
        dat=dat[abs(dat['phenoZ'])<norm.ppf(1-outcomeP,0,1)]
        # dat=dat[abs(dat['geneZ'])<norm.ppf(1-outcomeP,0,1)] # also make sure no eQTLs to rule out horizontal pleiotropy causing inflation
        if dat.shape[0]<(minM):
            continue
        # 3) only use every kth SNP to speed things up and ensure I don't need to do as much removing bc of high LD later.
        # I want approximately 100 SNPs before LD pruning so the k in kth is dynamically chosen based on how much data is originally available
        roundFactor=math.floor(dat.shape[0]/100)
        if roundFactor==0:
            roundFactor=1
        dat=dat.iloc[::roundFactor,:].drop_duplicates('geneSNP')
        dat['geneSNP'].to_csv(outDir1,header=False, index=False, sep=" ") # writing out list of SNPs to give PLINK to make LD calculations
        cmd=[callPlink(), "--r", "square", "--bfile", ldRefDir, "--extract", outDir1, "--out", outDir2] # the command for PLINK
        out=subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) # executing the command in the terminal
        ldmat=numpy.genfromtxt(outDir2+'.ld', dtype='f') # reading LD matrix back in
        mask=(numpy.isnan(ldmat).sum(axis=0)==ldmat.shape[0])==False # if any missing (bc MAF==0), will have all nan's for that SNP
        ldmat=ldmat[mask,:][:,mask]
        dat=dat[numpy.array(mask).astype(bool)]
        # 4) prune
        prunedLD,snpsKept=LDShrinkUpper(ldmat,ldUpperLimit,dat['geneSNP'].values)
        dat=dat[dat['geneSNP'].isin(snpsKept)]
        # 5) MR + test stats
        m=dat.shape[0]
        bx=numpy.column_stack(([1]*m,dat['geneZ'].values))
        by=dat['phenoZ'].values
        res=scipy.stats.spearmanr(dat[['geneZ','phenoZ']])
        ccor=res[0]
        pcorr=res[1]
        # if pcorr<0.05:
        #     continue
        if m<minM: # if small sample size or the correlation between bx and by is ~significant; changed to be 0 bc I actually dont
            continue
        counter=counter+1
        usedGenes.append(ggs[_])
        ms.append(m)
        ccs.append(ccor)
        ccps.append(pcorr)
        # print(_)
        if counter>nGenesToTry:
            break
        prunedLD,alpha,deco=regularizeLD(prunedLD,m+1,m+1,0.1,epsilon=1e-4,toCorr=True,justMCP=True)
        Theta=deco[1]@numpy.diag(1/deco[0])@deco[1].T
        BTB=numpy.linalg.inv(bx.T@Theta@bx)
        BTa=bx.T@Theta@by
        est=BTB@BTa
        D=numpy.diag((by-bx@est).squeeze())
        BT=bx.T@Theta
        rob=BTB@BT@D@prunedLD@D@BT.T@BTB
        # 6)
        se=numpy.diag(rob)**0.5
        stat=((est/se)**2)[1] # only care about causal estimate, not intercept
        cstats.append(stat)
    # 7)
    infl=numpy.median(cstats)/stats.chi2.ppf(0.5,1) if len(cstats)>10 else 1
    # finding BP positions of genes used in inflation calculation
    bps=genedf[genedf['Gene'].isin(usedGenes)]['geneBP'].values
    return infl,ccs,ccps,usedGenes,ms,bps

