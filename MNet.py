import numpy
import scipy
import networkx as nx
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

def entropyloss(A,B,eps=1e-8):
    C=A@B
    a=numpy.trace(C)
    b,o=numpy.linalg.eig(C)
    b=numpy.real(b)
    mask=b>eps
    b=sum(numpy.log(b[mask]))
    return a-b-A.shape[1]

def Penamatrix(H2):
    S=numpy.eye(len(H2))*0
    for i in range(0,len(H2)):
        for j in range(0,len(H2)):
            S[i,j]=1/(H2[i]*H2[j])**0.5
            S[j,i]=S[i,j]
    numpy.fill_diagonal(S,0)
    return S

def spearmancov(A):
    p=A.shape[1]
    s=numpy.arange(0,p)
    # s=scipy.stats.median_abs_deviation(A,axis=0)
    s=numpy.std(A,axis=0)
    s=numpy.array(s)
    R=scipy.stats.spearmanr(A,axis=0).correlation
    R=2*numpy.sin(R*numpy.pi/6)
    S=numpy.diag(s)@R@numpy.diag(s)
    return S

def soft(a,b):
    c=abs(a)-b
    mask=(c<0)
    c[mask]=0
    return c*numpy.sign(a)

def mcp(a,lam,ga=3):
    b=abs(a)
    z=soft(a,lam)/(1-1/ga)
    mask=b>(ga*lam)
    z[mask]=a[mask]
    return z

def MCPthreshold(S,lam,ga=3):
    p=S.shape[0]
    d=numpy.diag(S)**0.5
    R=cov2cor(S)
    s=R.reshape((p*p,))
    s1=abs(s)
    s2=mcp(s1,lam,ga)
    s2=s2*numpy.sign(s)
    S1=s2.reshape((p,p))
    numpy.fill_diagonal(S1,1)
    D=numpy.diag(d)
    S2=D@S1@D
    return S2

def cov2cor(A):
    A_=A
    d=numpy.diag(A_); 
    mask=d<0; d[mask]=1
    D=numpy.diag(1/d**0.5)
    return D@A@D

def matrixsqrt(A):
    d,v=numpy.linalg.eig(A)
    d1=d*0
    mask=d>0
    d1[mask]=1/d[mask]
    d=d**0.5
    d1=d1**0.5
    A=v@numpy.diag(d)@v.T
    B=v@numpy.diag(d1)@v.T
    return A,B

def positiveadj(A,min_eps=0.001):
    d,v=numpy.linalg.eig(A)
    mask=d<min_eps
    d[mask]=0
    B=v@numpy.diag(d)@v.T
    return B

def rowMedian(A):
    g=[numpy.median(A[_,:]) for _ in range(0,A.shape[0])]
    return g

def entropy_mcp_spearman_sampling(BETA,SE,Rnoise,lamvec=numpy.linspace(2,15,10)/100,
                                  max_eps=0.001,max_iter=25,rho=0.25,
                                  min_eig=0.01,subtime=100,subfrac=0.66,subthres=0.95):
    m,p=BETA.shape
    #########################tuning parameter selection############################## 
    subvecerror=numpy.zeros((len(lamvec),subtime))
    Thetalist=numpy.zeros((len(lamvec),subtime,p,p))
    
    for j in range(0,subtime):
        indsub=numpy.random.choice(m,size=round(subfrac*m),replace=False)
        mask=numpy.zeros((m,),dtype='bool'); mask[indsub]=True
        BETAS=BETA[mask]
        SES=SE[indsub]
        S=spearmancov(BETAS)*BETAS.shape[0]
        M=S*0
        for ii in range(0,BETAS.shape[0]):
            M=M+SES[ii,:]*(Rnoise*SES[ii,:]).T
        S1=cov2cor(S-M)
        Theta0=numpy.linalg.inv(MCPthreshold(S1,0.1))
        BETAS2=BETA[mask==False,:]
        SES2=SE[mask==False,:]
        S2=spearmancov(BETAS2)*BETAS2.shape[0]
        M2=S2*0
        for ii in range(0,BETAS2.shape[0]):
            M2=M2+SES[ii,:]*(Rnoise*SES2[ii,:]).T
        S2=cov2cor(S2-M2)
        
        for i in range(0,len(lamvec)):
            Theta=Theta0.copy()
            Theta1=Theta0*0
            Delta1=Theta.copy()
            Gamma1=Delta1*0
            Delta2=Theta.copy()
            Gamma2=Delta2*0
            error=numpy.linalg.norm(Theta-Theta1,ord='fro')/p**0.5
            iter=0
            while (error>max_eps) & (iter<max_iter):
                Theta1=Theta.copy()
                Q=S1+Gamma1-rho*Delta1
                Theta=(-Q+matrixsqrt(Q@Q+4*rho*numpy.eye(p))[0])/2/rho
                Delta1=mcp((Theta+Gamma1/rho).reshape((p*p,)),lamvec[i]/rho,ga=3).reshape((p,p))
                Gamma1=Gamma1+rho*(Theta-Delta1)
                if (iter%5==0) & (min(numpy.linalg.eig(Theta)[0])<min_eig):
                    Q=S1+Gamma1+Gamma2-rho*(Delta1+Delta2)
                    Theta=(-Q+matrixsqrt(Q@Q+8*rho*numpy.eye(p))[0])/4/rho
                    Delta1=mcp((Theta+Gamma1/rho).reshape((p*p,)),lamvec[i]/rho,ga=3).reshape((p,p))
                    Gamma1=Gamma1+rho*(Theta-Delta1)
                    Delta2=positiveadj(Theta+Gamma2/rho,min_eps=min_eig)
                    Gamma2=Gamma2+rho*(Theta-Delta2)
                # end if
                iter=iter+1
                error=numpy.linalg.norm(Theta-Theta1,ord='fro')/p**0.5
            # end while
            df=sum((sum(Delta1!=0)+p)/2)
            subvecerror[i,j]=entropyloss(S2,Delta1)+numpy.log(m)/m*df
            Theta0=Delta1.copy()
            Thetalist[i,j,:,:]=Theta0
        # end for lamvec
    # end for subtime
    istar=numpy.argmin(rowMedian(subvecerror))
    Thetalist=Thetalist[istar,:,:,:]
    K=Thetalist[1,:,:]*0
    Thetasub=K.copy()
    for i in range(0,subtime):
        K=K+sum(sum(Thetalist[i,:,:]!=0))/subtime
        Thetasub=Thetasub+Thetalist[i,:,:]/subtime
    S=spearmancov(BETA)*BETA.shape[0]
    M=S*0
    for i in range(0,BETA.shape[0]):
        M=M+SE[i,:]*(Rnoise*SE[i,:]).T
    S1=cov2cor(S-M)
    Theta0=numpy.linalg.inv(MCPthreshold(S1,0.1))
    Theta=Theta0.copy()
    Theta1=Theta0*0
    Delta1=Theta.copy()
    Gamma1=Delta1*0
    Delta2=Theta.copy()
    Gamma2=Delta2*0
    error=numpy.linalg.norm(Theta-Theta1,ord='fro')/p**0.5
    iter=0
    while (error>max_eps) & (iter<(2*max_iter)):
        Theta1=Theta.copy()
        Q=S1+Gamma1-rho*Delta1
        Theta=(-Q+matrixsqrt(Q@Q+4*rho*numpy.eye(p))[0])/2/rho
        Delta1=mcp((Theta+Gamma1/rho).reshape((p*p,)),lamvec[istar]/rho,ga=3).reshape((p,p))
        Gamma1=Gamma1+rho*(Theta-Delta1)
        if (iter%5==0) & (min(numpy.linalg.eig(Theta)[0])<min_eig):
            Q=S1+Gamma1+Gamma2-rho*(Delta1+Delta2)
            Theta=(-Q+matrixsqrt(Q@Q+8*rho*numpy.eye(p))[0])/4/rho
            Delta1=mcp((Theta+Gamma1/rho).reshape((p*p,)),lamvec[istar]/rho,ga=3).reshape((p,p))
            Gamma1=Gamma1+rho*(Theta-Delta1)*(K>subthres)
            Delta2=positiveadj(Theta+Gamma2/rho,min_eps=min_eig)
            Gamma2=Gamma2+rho*(Theta-Delta2)
        # end for
        iter=iter+1
        error=numpy.linalg.norm(Theta-Theta1,ord='fro')/p**0.5
    # end while
    Theta0=Delta1.copy()
    A=dict()  
    A['Theta']=Theta0
    A['Theta_Alt']=Thetasub
    A['K']=K
    A['cv_error']=subvecerror
    return A

def addEdges(edges):
    p=edges.shape[0]
    G=nx.Graph()
    for i in range(0,p):
        for j in range(0,p):
            cond=(j>=i) | (edges[i,j]==0)
            if cond:
                continue
            G.add_edges_from([(j,i)])
    return G        

def colorCoreGenes(edges):
    color_map=['powderblue']; p=edges.shape[0]
    for node in range(1,p):
        nodek=edges[:,node]
        if nodek[0]==1:
            color_map.append('tomato')
        else: 
            color_map.append('lightgrey') 
    return color_map

def mapNodeNames(edges,newlabs):
    p=edges.shape[0]
    if len(newlabs)!=p:
        raise ValueError('size of edge matrix does not equal length of new node names')
    old=range(0,p)
    new={}
    for _ in range(0,p):
        new[old[_]]=newlabs[_]
    return new

















