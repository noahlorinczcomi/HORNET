#!/usr/bin/env python
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from functions import *
import argparse
import os
import time

parser=argparse.ArgumentParser(prog='HORNET: Genome-wide robust multivariable MR with gene expression',
                    description='This program performs univariable and multivariable MR with MR-Jones to identify causal genes and their regulatory networks.')

### data
parser.add_argument('--eQTLGWAS',action='store',type=str,help='(Required) path to directory/folder containing multiple files with .gz extensions that are the eQTL GWAS summary statistics for gene expression, where each file corresponds to each chromosome. There should be nothing else in this directory/folder other than these files. If you only want to perform the analysis on certain chromosomes, create a directory/folder containg the eQTL GWAS summary statistics for only those chromosomes.')
parser.add_argument('--phenoGWAS',action='store',type=str,help='(Required) filepath with .gz extension to GWAS for outcome phenotype')
parser.add_argument('--LDRef',action='store',type=str,help='(Required) filepath without extension to the LD reference panel. The .bed, .bim, .fam files must each be present.')
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
parser.add_argument('--POtherLoci',action='store',type=float,default=0.05,help='(Optional) IVs that are in absolute LD beyond -LDOhterLoci with SNPs in other loci that have P-values less than this threshold may be excluded (see -OtherLociWindow). The default is 0.05.')
parser.add_argument('--otherLociMbWindow',action='store',type=float,default=0.05,help='(Optional) IVs that meet the conditions of -LDOtherLoci and -POtherLoci and are within this Mb window will be excluded. The default is 1, corresponding to 1Mb.')
parser.add_argument('--MVMRIVPThreshold',action='store',type=float,default=5e-8,help='(Optional) IVs must have P<this threshold in a joint test for multiple genes to be considered as IVs in multivariable MR. The default is 5e-8.')
parser.add_argument('--UVMRIVPThreshold',action='store',type=float,default=5e-5,help='(Optional) IVs must have P<this threshold in the eQTL GWAS to be considered as IVs in univariable MR. The default is 5e-5.')
### estimation
parser.add_argument('--nSNPsBiasCorrection',action='store',type=int,default=50,help='(Optional) The minimum number of non-significant SNPs to be used when calculating bias-correction terms used by MR-Jones, following the methods of Lorincz-Comi et al. (2023). The default is 50.')
parser.add_argument('--PvalueThresholdBiasCorrection',action='store',type=float,default=0.01,help='(Optional) Only SNPs with P<this threshold will be used to calculatee the bias-correction terms in MR-Jones, following the methods of Lorincz-Comi et al. (2023). The default is 0.01.')
parser.add_argument('--assumeNoSampleOverlap',action='store',type=str,default='yes',help='(Optional) Should it be assumed that the eQTL and disease GWAS do not share any of the same participants? If this assumption should not be made, bias-correction for sample overlap will be made following the methods of Lorincz-Comi et al. (2023). The default is "yes".')
parser.add_argument('--minMVMRIVs',action='store',type=int,default=30,help='(Optional) The minimum number of IVs that may be used in multivariable MR with MR-Jones. If less than this threshold are available for the locus, the locus will be skipped. The default is 30.')
parser.add_argument('--hessMinScale',action='store',type=float,default=2.5,help='(Optional) The signal-to-noise ratio for associations between SNPs and gene expression. Larger values will restrict the analysis to only genes with the strongest signals. Values equal to or less than 1/2 may lead to a Hessian matrix with diagonal elements that are not the correct sign, which is introduced by the MRBEE/MR-Jones correction for weak instruments. Note that the bias-correction term in the Hessian matrix of MR-Jones is automatically multiplied by the factor 1/2. The default is 2.5, corresponding to an effective signal-to-noise ratio of 5.')
### post-estimation stuff
parser.add_argument('--adjustSEsForInflation',action='store',type=str,default='true',help='(Optional) One of True or False. Should the standard errors of causal effect estimates be adjusted for inflation due to misspecified LD structure? The default is True.')
parser.add_argument('--saveRawData',action='store',type=str,default='false',help='(Optional) One of True or False. Should the raw data used in multivariable MR be saved? This includes association estimates with each gene in a group, their standard errors, correspondingly the same for the outcome phenotype, and the estimated LD matrix.')
parser.add_argument('--out',action='store',type=str,default='results',help='(Optional) filepath without extension of location in which results should be written, in tab-separated .txt format. The default is results. Note that a tab-separated file named diagnostics.txt will automatically be written. This file contains information about the performance of multivariable MR.')
parser.add_argument('--iterativelySave',action='store',type=str,default='true',help='(Optional) Should results be saved iteratively as each chromosome is completed? This is helpful if you anticipate the analysis may take a relatively long time and you do not want to lose progress in case your access to the machine it is running on expires. The default is True.')
### flags related to what is printed or not
parser.add_argument('--silence',action='store',type=str,default='no',help='(Optional) Should warnings about the size of the CHP window outside of the target locus be ignored? Put True or yes. The default is "no".')
### done
print(' ')
print('HORNET started '+time.ctime())

args=parser.parse_args()

### first checking if all of the files they gave actually exist
fileChecker(args.eQTLGWAS, 'eQTL GWAS directory')
fileChecker(args.phenoGWAS, 'phenotype GWAS filepath')
fileChecker(args.LDRef+'.bed', 'LD reference panel filepath (w/o extension)')
fileChecker(args.writableDir, 'temporing working directory')
od=(args.out.split('/'))[:-1]
fileChecker('/'.join(od), 'output/results directory')
# if no files in eQTL GWAS directory, there's a problem
if (os.listdir(args.eQTLGWAS))==0:
    raise ValueError('It looks like there are no files in the eQTL GWAS directory ({})'.format(args.eQTLGWAS))

## 
candidateGenes=args.candidateGenes.split(',') if len(args.candidateGenes)>0 else []
snplabs=args.snpLabels.split(',')
zlabs=args.zLabels.split(',')
ealabs=args.effectAlleles.split(',')
startwd=os.getcwd() # starting working directory - will be the one containing the HORNET software

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
mapwd='data/maps/1kgPhase3maps'
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

### running
# os.chdir(dirGene)
runningres=pandas.DataFrame() # to be filled in
runningdiagnostics=pandas.DataFrame() # to be filled in
bim=pandas.read_csv(ldRefDir+'.bim',sep='\t',names=['chr','rsid','x','bp','a1','a2']) # for figuring out which chromosome is being used later
infls=list()
print('Beginning the analysis')
for _ in range(0, len(os.listdir(dirGene))):
    t0=time.perf_counter() # start timer
    # os.chdir(dirGene)
    fpGene=os.listdir(dirGene)[_] # define filepath for exposure data
    dataGene=loadExposureGWASData(dirGene+'/'+fpGene,effectAlleleGene,zGene,rsidGene,snpBPGene,geneLabelGene,geneBPLabelGene,nLabel=nLabel,isRawGzippedGTEx=isRawGTEx,mapwd=mapwd) # load genes
    merged=mergeExposureAndOutcomeGWAS(dataGene,dataPheno) # merge with already loaded (outside of the loop) outcome GWAS data
    del dataGene; # delete data we no longer need
    mm=pandas.merge(bim,merged['geneSNP'].drop_duplicates(),left_on='rsid',right_on='geneSNP') # find chromosome
    chromosome=mm['chr'].values[0]
    print('Starting the analysis with chromosome '+str(chromosome))
    if analysisInPhenotypeLoci==True: # form gene groups
        geneGroups,ggKeys,lens,usedGenes=defGeneGroupsByOutcome(q0geneGroups, merged, KbWindow=outcomeClumpingKBWindow,closestK=numIndexGenesToFormGroup)
    else:
        geneGroups,ggKeys,lens,usedGenes=defGeneGroups(q0geneGroups,merged)
    # [print(geneGroupFinder(geneGroups,candidateGenes[i],isGtex=True if isRawGTEx else False)) for i in range(0,len(candidateGenes))]
    thingsMonitored,outerDict=MVMRworkhorse(merged,geneGroups,ggKeys,writableDir=writableDir,ldRefDir=ldRefDir,isGtex=isRawGTEx,
                                            analysisOnlyInOutcomeLoci=analysisInPhenotypeLoci,outcomeLociMbWindow=outcomeLociMbWindow,
                                            ldUpperLimit=ldUpperLimit, ldOtherLociOtherPt=ldOtherLociOtherPt,ldOtherLociR2=ldOtherLociR2,
                                            ldOtherLociWindow=ldOtherLociWindow, q0Correls=q0Correls,nMinCorrels=nMinCorrels,jointChiGenesP=jointChiGenesP,
                                            assumedMissingMean=assumedMissingMean,opAlpha=opAlpha,verbose=verbose,nMinIVs=nMinIVs,
                                            hessMinScale=hessMinScale,silence=silence, UniMRIVPThreshold=UniMRIVPThreshold, candidateGenes=candidateGenes,
                                            assumeNoSampleOverlap=assumeNoSampleOverlap,shrinkBiasCorrection=shrinkBiasCorrection,saveData=saveData)
    if len(outerDict)==0:
        continue
    res=organizeMetaResults(outerDict); res['CHR']=chromosome # store results
    diagnostics=organizeThingsMonitored(thingsMonitored) # store diagnostics
    runningdiagnostics=runningdiagnostics.append(diagnostics.copy())
    if adjustSEsForInflation: # should SEs be adjusted for inflation?
        infl,ccs,ccps,ug,ms,bps=calcInflationFactor(merged,ldRefDir,writableDir,mbWindow=1,ldUpperLimit=ldUpperLimit,outcomeP=0.01,nGenesToTry=100,geneDistMb=1/2,minM=50)
        infls.append(infl)
        res['infl']=numpy.max((1,infl))
        res=adjustInflation(res,numpy.max((1,infl)))
    runningres=runningres.append(res.copy())
    # LDLeadGeneSNPs(res,thischr,writableDir,ldRefDir) # write out results to be read in by an R program that will make figures etc.    
    t1=time.perf_counter()-t0
    print('\t Chromosome '+str(chromosome)+' complete '+str(round(t1/60,1))+' minutes later')
    # iteratively save so I don't lose any progress between genes
    if (args.iterativelySave.lower()=='true') | (args.iterativelySave.lower()=='yes'):
        runningres.to_csv(args.out+'.txt',sep='\t')
        runningres.to_csv(args.out+'_diagnostics.txt',sep='\t')

fp1=args.out+'.txt'
fp2=args.out+'_diagnostics.txt'
runningres.to_csv(fp1,sep='\t')
runningres.to_csv(fp2,sep='\t')
print('Results are written to '+fp1)
print('Analysis diagnostics are written to '+fp2)

# remove unnecessary files
# now I want to remove these files I just wrote out because they may be large
out=subprocess.call(["rm", writableDir+'/tempOut.ld'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
out=subprocess.call(["rm", writableDir+'/tempOut.log'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
out=subprocess.call(["rm", writableDir+"/myExtract.txt"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
out=subprocess.call(["rm", writableDir+"/outcomeplinkout.clumped"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
out=subprocess.call(["rm", writableDir+"/outcomeplinkout.log"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
out=subprocess.call(["rm", writableDir+"/outcomePsout.txt"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

# network construction
## steps:
# 1) write out file to be read in by R
# 2) execute R program
# 3) 