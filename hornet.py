#!/usr/bin/env python
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from functions import *
from MNet import * # this is a source file in the HORNET directory
import argparse

parser=argparse.ArgumentParser(prog='HORNET: Genome-wide robust multivariable MR with gene expression',
                    description='This program performs univariable and multivariable MR with MR-Jones to identify causal genes and their regulatory networks.')

### data
parser.add_argument('--eQTLGWAS',action='store',type=str,help='(Required) path to directory/folder containing multiple files with .gz extensions that are the eQTL GWAS summary statistics for gene expression, where each file corresponds to each chromosome. There should be nothing else in this directory/folder other than these files. If you only want to perform the analysis on certain chromosomes, create a directory/folder containg the eQTL GWAS summary statistics for only those chromosomes.')
parser.add_argument('--phenoGWAS',action='store',type=str,help='(Required) filepath with .gz extension to GWAS for outcome phenotype')
parser.add_argument('--LDRef',action='store',type=str,default=os.getcwd()+'/data/ldref/1kg3TRANS',help='(Required) filepath without extension to the LD reference panel. The .bed, .bim, .fam files must each be present.')
parser.add_argument('--isRawGTEx',action='store',type=str,help='(Required) If the eQTL association estimates are directly from the "eQTL Tissue-Specific All SNP Gene Associations" data sets from the GTEx Consortium, put True; else, put False.')
parser.add_argument('--writableDir',action='store',type=str,default='tempfiles',help='(Optional) A directory/folder in which temporary files (eg, LD matrices) can be written and removed.')
### data labels
parser.add_argument('--snpLabels',action='store',type=str,help='(Required,Required) comma-separated list of rsID column names for gene expression (first) and outcome phenotype (second) GWAS, respectively. If either of the data do not contain an rsID column, please put NA in the list.')
parser.add_argument('--eQTLSNPBPLabel',action='store',type=str,help='(Required) column name of SNP base pair position values in the eQTL GWAS data sets')
parser.add_argument('--zLabels',action='store',type=str,help='(Required,Required) comma-separated list of Z-statistic column names for gene expression (first) and outcome trait (second) GWAS, respectively. If Z-statistics are not present in one or either data set, please input BETA@SE where BETA is the name of the effect size column and SE is the name of the standard error column,')
parser.add_argument('--effectAlleles',action='store',type=str,help='(Required, Required) column names in the eQTL GWAS (first) and phenotype GWAS (second) representing the effect alleles, respectively')
parser.add_argument('--eQTLGeneLabel',action='store',type=str,help='(Required) the name of the column in the eQTL GWAS that indicates the gene name/label')
parser.add_argument('--eQTLGeneBP',action='store',type=str,help='(Required) the name of the column in the eQTL GWAS that indicates the base pair location of the centers of genes or their transcription start sites')
parser.add_argument('--eQTLGWASn',action='store',type=int,default=500,help='(Optional) the sample size of the eQTL GWAS')
parser.add_argument('--phenotypeGWASn',action='store',type=int,default=30000,help='(Optional) the sample size of the outcome phenotype GWAS')
### type of analysis
parser.add_argument('--candidateGenes',action='store',type=str,default='',help='(Optional) If you only want to perform the analysis for a list of candidate genes, please specify that comma-separated list here. If your list of candidate genes is large, you can place the directory to that file here. Each row of this file should contain gene IDs that are present in the eQTL GWAS data.')
parser.add_argument('--analysisInPhenotypeLoci',action='store',type=str,default='true',help='(Optional) One of True or False. Should the the analyses only be performed in loci that are significantly associated with the outcome phenotype? See -phenoLociPvalue and -phenoLociMbWindow for how these loci are defined. The default is True.')
parser.add_argument('--phenoLociPvalue',action='store',type=float,default=5e-5,help='(Optional) P-value threshold for considering a locus to be significantly associated with the outcome phenotype. The default is 5e-5.')
parser.add_argument('--phenoLociMbWindow',action='store',type=float,default=0.5,help='(Optional) Mb window size for considering two phenotype loci as distinct. The default is 0.5.')
parser.add_argument('--phenoLociR2',action='store',type=float,default=0.1,help='(Optional) SNPs associated with the outcome at significance level -phenoLociPvalue and in squared LD beyond this thresold will be clumped together using PLINKs clumping procedure. The default is 0.1.')
### forming gene groups
parser.add_argument('--geneGroupPvalueThreshold',action='store',type=float,default=5e-8,help='(Optional) only genes sharing cis-eQTLs with P<this threshold will be grouped together for multivariable MR. The default is 5e-8, and we recommend a threshold at or very close to genome-wide significance (5e-8) because this ensures that gene groups do not become exceedingly large. If they become exceedingly large, causal estimation may be unstable and there is a greater dependence on reliable missing data imputation.')
parser.add_argument('--onlyWithinMbOfPhenoLocus',action='store',type=float,default=1,help='(Optional) Only form gene groups within this many Mb of an outcome locus if -analysisInPhenotypeLoci is True. The default is 1 (1Mb).')
parser.add_argument('--numIndexGenesToFormGroup',action='store',type=int,default=4,help='(Optional) If -analysisInPhenotypeLoci is True, the -numIndexGenesToFormGroup closest to each lead SNP indexing a phenotype locus will be used to form a group of genes. Note, these genes must also be within -onlyWithinMbOfPhenoLocus of the lead phenotype SNP for the locus. The default is 4.')
### forming IV sets
parser.add_argument('--ldMax',action='store',type=float,default=0.5,help='(Optional) two cis-eQTLs with an LD correlation greater than this threshold will be considered as the same cis-eQTL and the one with the largest eQTL GWAS P-value will be excluded')
parser.add_argument('--LDOtherLoci',action='store',type=float,default=0.2,help='(Optional) IVs that are in absolute LD beyond this threshold with SNPs in other loci may be excluded (see -POtherLoci). The default is 0.2.')
parser.add_argument('--POtherLoci',action='store',type=float,default=0.05,help='(Optional) IVs that are in absolute LD beyond -LDOtherLoci with SNPs in other loci that have P-values less than this threshold may be excluded (see -OtherLociWindow). The default is 0.05.')
parser.add_argument('--otherLociMbWindow',action='store',type=float,default=0.05,help='(Optional) IVs that meet the conditions of -LDOtherLoci and -POtherLoci and are within this Mb window will be excluded. The default is 1, corresponding to 1Mb.')
parser.add_argument('--MVMRIVPThreshold',action='store',type=float,default=5e-8,help='(Optional) IVs must have P<this threshold in a joint test for multiple genes to be considered as IVs in multivariable MR. The default is 5e-8.')
parser.add_argument('--UVMRIVPThreshold',action='store',type=float,default=5e-5,help='(Optional) IVs must have P<this threshold in the eQTL GWAS to be considered as IVs in univariable MR. The default is 5e-5.')
### estimation
parser.add_argument('--imputeMissingEQTLs',action='store',type=str,default='yes',help='(Optional) Should missing values in the MVMR IV set (which is a union set of gene-specific IV sets) be imputed using the matrix completion method described in Lorincz-Comi & Yang et al (2023)? If you put "no" or "False" here, missing values will still be imputed, but with all 0s. The default is "yes".')
parser.add_argument('--nSNPsBiasCorrection',action='store',type=int,default=50,help='(Optional) The minimum number of non-significant SNPs to be used when calculating bias-correction terms used by MR-Jones, following the methods of Lorincz-Comi & Yang et al. (2023). The default is 50.')
parser.add_argument('--PvalueThresholdBiasCorrection',action='store',type=float,default=0.01,help='(Optional) Only SNPs with P<this threshold will be used to calculatee the bias-correction terms in MR-Jones, following the methods of Lorincz-Comi & Yang et al. (2023). The default is 0.01.')
parser.add_argument('--assumeNoSampleOverlap',action='store',type=str,default='yes',help='(Optional) Should it be assumed that the eQTL and disease GWAS do not share any of the same participants? If this assumption should not be made, bias-correction for sample overlap will be made following the methods of Lorincz-Comi & Yang et al. (2023). The default is "yes".')
parser.add_argument('--minMVMRIVs',action='store',type=int,default=30,help='(Optional) The minimum number of IVs that may be used in multivariable MR with MR-Jones. If less than this threshold are available for the locus, the locus will be skipped. The default is 30.')
parser.add_argument('--hessMinScale',action='store',type=float,default=2.5,help='(Optional) The signal-to-noise ratio for associations between SNPs and gene expression. Larger values will restrict the analysis to only genes with the strongest signals. Values equal to or less than 1/2 may lead to a Hessian matrix with diagonal elements that are not the correct sign, which is introduced by the MRBEE/MR-Jones correction for weak instruments. Note that the bias-correction term in the Hessian matrix of MR-Jones is automatically multiplied by the factor 1/2. The default is 2.5, corresponding to an effective signal-to-noise ratio of 5.')
### post-estimation stuff
parser.add_argument('--adjustSEsForInflation',action='store',type=str,default='yes',help='(Optional) One of True or False. Should the standard errors of causal effect estimates be adjusted for inflation due to misspecified LD structure? The default is "yes".')
parser.add_argument('--saveRawData',action='store',type=str,default='false',help='(Optional) One of True or False. Should the raw data used in multivariable MR be saved? This includes association estimates with each gene in a group, their standard errors, correspondingly the same for the outcome phenotype, and the estimated LD matrix.')
parser.add_argument('--out',action='store',type=str,default='results',help='(Optional) filename without extension of location in which results should be written, in tab-separated .txt format. The default is results. Note that a tab-separated file named diagnostics.txt will automatically be written. This file contains information about the performance of multivariable MR.')
parser.add_argument('--iterativelySave',action='store',type=str,default='true',help='(Optional) Should results be saved iteratively as each chromosome is completed? This is helpful if you anticipate the analysis may take a relatively long time and you do not want to lose progress in case your access to the machine it is running on expires. The default is True.')
### creating networks
parser.add_argument('--graphNetworks',action='store',type=str,default='yes',help='(Optional) Should graph images of causal networks be stored? If "yes", they will go into the "plots" subfolder. The default is "yes".')
parser.add_argument('--graphNetworksFilePrefix',action='store',type=str,default='',help='(Optional) The default prefix for the filenames of graphed regulatory networks. Note that all network graph files will contain the lead gene in the filename automatically. What you put here will come before that. The default is "".')
parser.add_argument('--networkR2Thres',action='store',type=float,default=0.25,help='(Optional) Only loci with at least this much genetic variance in the disease phenotype explained will be considered for graph analysis using MNet. The default is 0.25.')
parser.add_argument('--networkDiseaseLabel',action='store',type=str,default='disease',help='(Optional) In graphs of causal networks, the disease node will have this label. This default is "disease"')
parser.add_argument('--networkGraphsOut',action='store',type=str,default='plots',help='(Optional) The directory/folder in which you would like plots of causal network graphs to be saved. The default is the "plots/" subdirectory/subfolder of the "hornet/" subdirectory/subfolder.')
### flags related to what is printed or not
parser.add_argument('--silence',action='store',type=str,default='no',help='(Optional) Should warnings about the size of the CHP window outside of the target locus be ignored? Put "true" or "yes". The default is "no".')
### done
print(' ')
print('HORNET started '+time.ctime())

args=parser.parse_args()

### first checking if all of the files/direcetories they gave actually exist
fileChecker(args.eQTLGWAS, 'eQTL GWAS directory')
fileChecker(args.phenoGWAS, 'phenotype GWAS filepath')
fileChecker(args.LDRef+'.bed', 'LD reference panel filepath (w/o extension)')
fileChecker(args.writableDir, 'temporing working directory')
od=(args.out.split('/'))[:-1]
fileChecker('/'.join(od), 'output/results directory')
# if no files in eQTL GWAS directory, there's a problem
if (os.listdir(args.eQTLGWAS))==0:
    raise ValueError('It looks like there are no files in the eQTL GWAS directory ({})'.format(args.eQTLGWAS))
# output directory for graphs
graphoutdir=args.networkGraphsOut+'/' if args.networkGraphsOut!='plots' else os.getcwd()+'/plots/'
if os.path.isdir(graphoutdir)==False:
    raise ValueError('\n It looks like the directory you want to save network graph plots to does not exist({})'.format(graphoutdir))

## 
candidateGenes=args.candidateGenes.split(',') if len(args.candidateGenes)>0 else []
snplabs=args.snpLabels.split(',')
zlabs=args.zLabels.split(',')
ealabs=args.effectAlleles.split(',')
startwd=os.getcwd() # starting working directory - will be the one containing the HORNET software
impute=(args.imputeMissingEQTLs.lower()=='true') | (args.imputeMissingEQTLs.lower()=='yes')

### exposure parameters
dirGene=args.eQTLGWAS 
isRawGTEx=(args.isRawGTEx.lower()=='true') | (args.isRawGTEx.lower()=='yes')
effectAlleleGene=ealabs[0]
zGene=zlabs[0]
rsidGene=snplabs[0]
geneLabelGene=args.eQTLGeneLabel
geneBPLabelGene=args.eQTLGeneBP
snpBPGene=args.eQTLSNPBPLabel
nLabel=args.eQTLGWASn
# if GTEx, you need to choose these carefully or else no analysis can be done
ldUpperLimit=args.ldMax
ldOtherLociOtherPt=args.POtherLoci
ldOtherLociR2=(args.LDOtherLoci)**2;
ldOtherLociWindow=args.otherLociMbWindow;
q0Correls=1-norm.ppf(args.PvalueThresholdBiasCorrection)
nMinCorrels=args.nSNPsBiasCorrection
jointChiGenesP=args.MVMRIVPThreshold
assumedMissingMean=0;
opAlpha=0
verbose=True
nMinIVs=args.minMVMRIVs # default is 30
hessMinScale=args.hessMinScale # default is 2.5
silence=(args.silence.lower()=='true') | (args.silence.lower()=='yes')
UniMRIVPThreshold=args.UVMRIVPThreshold
adjustSEsForInflation=(args.adjustSEsForInflation.lower()=='true') | (args.adjustSEsForInflation.lower()=='yes')
analysisInPhenotypeLoci=(args.analysisInPhenotypeLoci.lower()=='true') | (args.analysisInPhenotypeLoci.lower()=='yes')
outcomeLociMbWindow=args.onlyWithinMbOfPhenoLocus
q0geneGroups=norm.ppf(1-args.geneGroupPvalueThreshold)
assumeNoSampleOverlap=(args.assumeNoSampleOverlap.lower()=='true') | (args.assumeNoSampleOverlap.lower()=='yes')
shrinkBiasCorrection=True # I don't think we want this adjusted; keep it simple
saveData=(args.saveRawData.lower()=='true') | (args.saveRawData.lower()=='yes')
numIndexGenesToFormGroup=args.numIndexGenesToFormGroup

### phenotype parameters
fpPheno=args.phenoGWAS
writableDir=args.writableDir
ldRefDir=args.LDRef
mapwd=os.path.dirname(os.path.dirname(ldRefDir))+'/maps/1kgPhase3maps'
rsidPheno=snplabs[1]
effectAllelePheno=ealabs[1]
zPheno=zlabs[1]
nLabel=args.phenotypeGWASn
outcomeClumpingKBWindow=args.phenoLociMbWindow
outcomeClumpingPthreshold=args.phenoLociPvalue
outcomeClumpingR2=args.phenoLociR2
### loading phenotype data
print('Loading phenotype GWAS data')
dataPheno=loadOutcomeGWASData(fpPheno,effectAllelePheno,zPheno,rsidPheno,ldRefDir)
dataPheno=findOutcomeSignals(dataPheno,ldRefDir,writableDir,outcomeClumpingKBWindow,outcomeClumpingPthreshold,outcomeClumpingR2)
### loading key/dictionary/lookup/reference data
bim=pandas.read_csv(ldRefDir+'.bim',sep='\t',names=['chr','rsid','x','bp','a1','a2']) # for figuring out which chromosome is being used later
genekey=pandas.read_csv(os.getcwd()+'/data/biomart_gene_legend.csv.gz',engine='python') # to be used only if candidateGenes are specified

### running
# os.chdir(dirGene)
runningres=pandas.DataFrame() # to be filled in
runningdiagnostics=pandas.DataFrame() # to be filled in
infls=list()
edgeLists={}
for _ in range(0, len(os.listdir(dirGene))):
    t0=time.perf_counter() # start timer
    # os.chdir(dirGene)
    fpGene=os.listdir(dirGene)[_] # define filepath for exposure data
    if len(candidateGenes)>0:
        ch=loadExposureGWASData(dirGene+'/'+fpGene,effectAlleleGene,zGene,rsidGene,snpBPGene,geneLabelGene,geneBPLabelGene,nLabel=nLabel,isRawGzippedGTEx=isRawGTEx,mapwd=mapwd,nr=10)
        mm=pandas.merge(bim,ch,left_on='rsid',right_on='geneSNP') # find chromosome
        mchr=mm['chr'][0]
        mm=mm[mm['Gene'].isin(candidateGenes)]
        del ch
        if mm.shape[0]==0:
            print('skipping chromsome '+str(mchr)+' because there are no candidate genes on this chromosome')
            continue
    dataGene=loadExposureGWASData(dirGene+'/'+fpGene,effectAlleleGene,zGene,rsidGene,snpBPGene,geneLabelGene,geneBPLabelGene,nLabel=nLabel,isRawGzippedGTEx=isRawGTEx,mapwd=mapwd) # load genes
    merged=mergeExposureAndOutcomeGWAS(dataGene,dataPheno) # merge with already loaded (outside of the loop) outcome GWAS data
    mm=pandas.merge(bim,merged['geneSNP'].drop_duplicates(),left_on='rsid',right_on='geneSNP') # find chromosome
    chromosome=mm['chr'].values[0]
    del dataGene; # delete data we no longer need
    print('Starting chromosome '+str(chromosome))
    if (analysisInPhenotypeLoci==True) & (len(candidateGenes)==0): # form gene groups
        geneGroups,ggKeys,lens,usedGenes=defGeneGroupsByOutcome(q0geneGroups, merged, KbWindow=outcomeClumpingKBWindow,closestK=numIndexGenesToFormGroup)
    elif len(candidateGenes)>0:
        if any(merged['Gene'].isin(candidateGenes))==False:
            continue
        geneGroups,ggKeys,lens,usedGenes=defineCandidateGeneGroups(merged,candidateGenes,MbWindow=2)
    else:
        geneGroups,ggKeys,lens,usedGenes=defGeneGroups(q0geneGroups,merged)
    # [print(geneGroupFinder(geneGroups,candidateGenes[i],isGtex=True if isRawGTEx else False)) for i in range(0,len(candidateGenes))]
    thingsMonitored,outerDict,edgeD=MVMRworkhorse(merged,geneGroups,ggKeys,writableDir=writableDir,ldRefDir=ldRefDir,isGtex=isRawGTEx,
                                                  analysisOnlyInOutcomeLoci=analysisInPhenotypeLoci,outcomeLociMbWindow=outcomeLociMbWindow,
                                                  ldUpperLimit=ldUpperLimit, ldOtherLociOtherPt=ldOtherLociOtherPt,ldOtherLociR2=ldOtherLociR2,
                                                  ldOtherLociWindow=ldOtherLociWindow, q0Correls=q0Correls,nMinCorrels=nMinCorrels,jointChiGenesP=jointChiGenesP,
                                                  opAlpha=opAlpha,nMinIVs=nMinIVs,
                                                  hessMinScale=hessMinScale,silence=silence, UniMRIVPThreshold=UniMRIVPThreshold, candidateGenes=candidateGenes,
                                                  assumeNoSampleOverlap=assumeNoSampleOverlap,shrinkBiasCorrection=shrinkBiasCorrection,networkR2Thres=args.networkR2Thres,impute=impute,saveData=saveData)
    edgeLists['Chromosome'+str(chromosome)]=edgeD # save edges
    if len(outerDict)==0:
        print('  Skipping chromosome '+str(chromosome)+' because of insufficient phenotype and/or eQTL signals')
        continue
    res=organizeMetaResults(outerDict); res['CHR']=chromosome; res['infl']=1; 
    res=prepRes(res) # store results
    diagnostics=organizeThingsMonitored(thingsMonitored) # store diagnostics
    diagnostics['Chromosome']=chromosome
    runningdiagnostics=pandas.concat([runningdiagnostics,diagnostics.copy()])
    if adjustSEsForInflation: # should SEs be adjusted for inflation?
        infl,ccs,ccps,ug,ms,bps=calcInflationFactor(merged,ldRefDir,writableDir,mbWindow=1,ldUpperLimit=ldUpperLimit)
        infls.append(infl)
        res['nullInflation']=numpy.max((1,infl))
        res=adjustInflation(res,numpy.max((1,infl)))
    runningres=pandas.concat([runningres,res.copy()])
    # LDLeadGeneSNPs(res,thischr,writableDir,ldRefDir) # write out results to be read in by an R program that will make figures etc.    
    t1=time.perf_counter()-t0
    print('  Chromosome '+str(chromosome)+' complete '+str(round(t1/60,1))+' minutes later')
    # iteratively save so I don't lose any progress between genes
    if (args.iterativelySave.lower()=='true') | (args.iterativelySave.lower()=='yes'):
        runningres.to_csv(args.out+'_tempresults.txt',sep='\t')
        runningdiagnostics.to_csv(args.out+'_diagnostics.txt',sep='\t')

fp1=args.out+'_results.txt'
fp2=args.out+'_diagnostics.txt'
runningres.to_csv(fp1,sep='\t')
runningdiagnostics.to_csv(fp2,sep='\t')
# delete iteratively saven files if user chose to iteratively save results
delcmd='del' if platform.system()=='Windows' else 'rm'
if (args.iterativelySave.lower()=='true') | (args.iterativelySave.lower()=='yes'):
    out=subprocess.call([delcmd, args.out+'_tempresults.txt'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
print('Results are written to '+fp1)
print('Analysis diagnostics are written to '+fp2)

# remove unnecessary files
# now I want to remove these files I just wrote out because they may be large
delcall=callDelete()
out=subprocess.call([delcall, writableDir+'/tempOut.ld'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
out=subprocess.call([delcall, writableDir+'/tempOut.log'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
out=subprocess.call([delcall, writableDir+"/myExtract.txt"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
out=subprocess.call([delcall, writableDir+"/outcomeplinkout.clumped"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
out=subprocess.call([delcall, writableDir+"/outcomeplinkout.log"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
out=subprocess.call([delcall, writableDir+"/outcomePsout.txt"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

### network construction (done in Python)
# first find top args.networksInTopKLoci loci
graphNetworks=(args.graphNetworks.lower()=='yes') | (args.graphNetworks.lower()=='true')
hasSaved=False
if (runningres.shape[1]>1) & (graphNetworks):
    vals=numpy.sort(runningres['RsquaredMRJones'].unique())[::-1]
    for _ in range(0,len(vals)):
        dcut=runningres[runningres['RsquaredMRJones']==vals[_]]
        genes=dcut['Gene'].values.tolist(); genes=[genes[_].split(',')[0] for _ in range(0,len(genes))]
        if dcut.shape[0]<1:
            continue
        chrom,grp=int(dcut['Chromosome'][0]),dcut['CHRspecificGroupID'][0]
        edges=edgeLists['Chromosome'+str(chrom)][grp]
        if sum(sum(edges))==0:
            continue
        G=addEdges(edges)
        fl=flatten_list([args.networkDiseaseLabel,genes])
        if len(fl)!=edges.shape[0]:
            print(dcut.head())
            print(dcut.shape)
            print(genes)
            print('length of names: '+str(len(fl))+', dimension of edge matrix: '+str(edges.shape[0]))
        G=nx.relabel_nodes(G,mapping=mapNodeNames(edges,fl)) # add new labels
        # color_map=colorCoreGenes(edges) # add colors for core and peripheral genes
        fig=matplotlib.pyplot.figure(figsize=(12,8))
        warnings.filterwarnings("ignore")
        nx.draw(G, with_labels=True,ax=fig.add_subplot(),alpha=0.5,font_size=10,node_color='skyblue') # all same color - less likely to produce an error
        leadgene=dcut.loc[abs(dcut['MRJonesEst']).idxmax()]['Gene']
        leadgene=leadgene.split('.')[0]
        fig.savefig(graphoutdir+args.graphNetworksFilePrefix+'_'+leadgene+'_graph.png',dpi=500)
        hasSaved=True
        ### run if using as an example
        # edges=numpy.array([[0,1,0],[1,0,1],[0,1,0]])
        # genes=['geneA','geneB']
        # leadgene='geneA'

print(' ')
if graphNetworks & hasSaved:
    print('Network graphs are saved in the '+graphoutdir+' folder')
print(' ')
print('HORNET ended '+time.ctime())

