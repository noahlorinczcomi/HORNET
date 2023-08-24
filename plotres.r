
############################################################################################################################################
# ensuring necessary package exists (ggplot2) and confirm working directory
############################################################################################################################################
if(!require('ggplot2')) install.packages('ggplot2')
if(!require('ggrepel')) install.packages('ggrepel')
source("modCMplot.r") # also using CMplot
# data must be present like this /HORNET/res.csv  --- I am going to make a cpres.py file to ensure this
############################################################################################################################################
# loading data
############################################################################################################################################
pat=normalizePath('res.csv')
if(!file.exists(pat)) stop('A file named res.csv of HORNET results is expected to be in the HORNET/ directory. This file is automatically saved when you use hornet.py, so it got deleted somewhow. You can copy any of your results files into the HORNET/ directory, rename the copy as res.csv, then run this program.')
data=read.csv(pat,sep='\t');
data=as.data.frame(data)
if(all.equal(data[,1],data[,2])) data=data[,-2]
data$Gene=sapply(data$Gene,function(h) unlist(strsplit(h,'[.]'))[1])
# I'm going to want to use gene names, not Ensembl IDs, throughout, so load the genekey data
# but only do this if the user used ensembl ids, if they didn't then skip this step
isens=sum(sapply(data$Gene,function(h) tolower(substr(h,start=1,stop=4))=='ensg'))/nrow(data)>0.5
if(isens) {
  genekey=read.csv('data/biomart_gene_legend.csv')
  colnames(genekey)[1]='Gene'
  genekey=as.data.frame(genekey)
  genekey=genekey[,c('Gene','GeneName')]
  data=merge(data,genekey,by='Gene',all.x=TRUE) # left join
  # remove ensembl id column and use gene name in its place
  data=data[,-which(colnames(data)=='Gene')]
  colnames(data)[which(colnames(data)=='GeneName')]='Gene'
}
############################################################################################################################################
# manhattan plots using CMplot (annotated P-values and rsquared)
############################################################################################################################################
mandf=data[,c('Gene','Chromosome','geneBP','RsquaredMRJones','MRBEEPostSelec_MVMR_Est','MRBEEPostSelec_MVMR_SE')]
mandf=na.omit(mandf)
if(nrow(mandf)==0) {
  mess='Could not create Manhattan plot because of insufficient data available'
} else {
  mandf$GeneP=1-pchisq((mandf$MRBEEPostSelec_MVMR_Est/mandf$MRBEEPostSelec_MVMR_SE)^2,1)
  traitdf=mandf[,c('Gene','Chromosome','geneBP','RsquaredMRJones','GeneP')]; colnames(traitdf)[3]='Position'
  traitdf$Gene=sapply(traitdf$Gene,function(h) unlist(strsplit(h,'[.]'))[1])
  ### Circular manhattan of MR-Jones R-squared (inner) and gene P-value (outer)
  # actually just rsquared
  thres=0.2
  topkgenes=traitdf[order(traitdf$RsquaredMRJones,decreasing=TRUE),]
  topkgenes=topkgenes[topkgenes$RsquaredMRJones>thres,]
  genes=c()
  uv=unique(topkgenes$RsquaredMRJones)[1:50] # top 50
  for(i in 1:length(uv)) {
      dv=topkgenes[round(topkgenes$RsquaredMRJones,3)==round(uv[i],3),]
      genes=c(genes,dv$Gene[which.min(dv$GeneP)])
  }
  setwd('plots')
  mess=paste0('Manhattan plot saved to ',normalizePath(getwd()),'/manhattan.png \n')
  suppressWarnings(
    CMplot(traitdf[,1:4], plot.type="m", band=0.5, LOG10=FALSE, ylab="Genetic variance explained by gene locus",threshold=thres,
        threshold.lty=2, threshold.lwd=2, threshold.col="black", amplify=TRUE, width=14,height=6,
        signal.col=NULL, chr.den.col=NULL, file="png",file.name="manhattan",dpi=300,file.output=TRUE,signal.line=NULL,
        cex=0.8,verbose=FALSE,
        highlight=genes,
        highlight.text.col=rep("black",length(genes)),
        highlight.col='black',col=c('cornflowerblue','gray60'),
        highlight.cex=1, highlight.text=genes,
        highlight.text.cex=3/4)
  )
}
cat(mess)
# rsquared and gene P
# mandf$GeneP=-log10(mandf$GeneP)
# nextmax=mandf$GeneP[which.max(mandf$GeneP[mandf$GeneP<Inf])]
# mandf$GeneP=ifelse(mandf$GeneP>nextmax,nextmax+5,mandf$GeneP)
# CMplot(traitdf,type="p",plot.type="c",r=0.4,col=c("grey30","grey60"),LOG10=FALSE,
#  threshold=-log10(5e-5),,cir.chr.h=1.5,amplify=TRUE,threshold.lty=1,threshold.col="black",
#  signal.line=NULL,signal.col=c("indianred","indianred"),chr.den.col=c("darkgreen","yellow","red"),
#  bin.size=1e6,outward=TRUE,file="png",file.name="manhattan",dpi=300,file.output=TRUE,verbose=FALSE,width=10,height=10)

############################################################################################################################################
# volcano plot
############################################################################################################################################
prattcut=0.2
qcut=qnorm(1-5e-5)
vdf=data[,c('Gene','Chromosome','geneBP','RsquaredMRJones','MRBEEPostSelec_MVMR_Est','MRBEEPostSelec_MVMR_SE','MRBEE_UVMR_Est')]
vdf=na.omit(vdf)
vdf$Gene=sapply(vdf$Gene,function(h) unlist(strsplit(h,'[.]'))[1])
vdf$z=vdf$MRBEEPostSelec_MVMR_Est/vdf$MRBEEPostSelec_MVMR_SE
vdf$pratt=vdf$MRBEEPostSelec_MVMR_Est*vdf$MRBEE_UVMR_Est
vdf$shapevar=NA
vdf$shapevar[(abs(vdf$z)>qcut) & (vdf$pratt<prattcut)]=1
vdf$shapevar[(abs(vdf$z)>qcut) & (vdf$pratt>prattcut)]=2
vdf$shapevar[(abs(vdf$z)<qcut) & (vdf$pratt>prattcut)]=3
vdf$shapevar[(abs(vdf$z)<qcut) & (vdf$pratt<prattcut)]=4
topkgenes=vdf[order(vdf$pratt,decreasing=TRUE),]
topkgenes=topkgenes[topkgenes$shapevar==2,]
topkgenes=topkgenes[topkgenes$pratt<1,]
topkgenes=topkgenes[order(topkgenes$pratt,decreasing=TRUE),]
genes=topkgenes$Gene[1:50] # top 50
vdf$genelabel=ifelse(vdf$Gene %in% genes, vdf$Gene,NA)

if(nrow(vdf)==0) {
  mess='Could not create Volcano plot because of insufficient data available'
} else {
  mess=paste0('Volcano plot saved to ',normalizePath(getwd()),'/volcano.png \n')
  mp=ggplot(vdf,aes(x=z,y=pratt,color=factor(shapevar))) +
    geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=0,ymax=prattcut),inherit.aes=FALSE,fill='gray93') +
    geom_rect(aes(xmin=-qcut,xmax=qcut,ymin=-Inf,ymax=Inf),inherit.aes=FALSE,fill='gray93') +
    lims(x=c(-20,20)) +
    theme_classic() +
    geom_vline(xintercept=0,linetype='dashed') +
    scale_color_manual(values=c('#47CCDE','#47DE6E','#DE4747','black')) + # blue, green, red, black
    #geom_hline(yintercept=prattcut) +
    #geom_hline(yintercept=0) +
    guides(shape='none',color='none',size='none') +
    theme(legend.position='bottom',legend.text=element_text(size=11)) +
    labs(x='Z-statistic for causal estimate',y='Pratt index') +
    geom_point() +
    geom_text_repel(aes(label=genelabel),data=vdf[vdf$shapevar==2,],size=2,angle=45,color='black',force=10,ylim=c(0.2,1)) +
    scale_y_continuous(breaks=seq(0,1,0.2),labels=seq(0,1,0.2),limits=c(-0.05,1))

  suppressWarnings(ggsave('volcano.png',mp,width=8,height=6))
}
cat(mess)
############################################################################################################################################
# clean up
############################################################################################################################################
sys=Sys.info()['sysname']
delo=ifelse(sys %in% c('Darwin','Linux'), 'rm', 'del')
setwd('..') # was in plots/, need to back up to HORNET/ where res.csv is stored
oo=system(paste0(delo,' ','res.csv'),intern=TRUE)





