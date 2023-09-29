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
parser.add_argument('--geneScreener',action='store',type=str,default='mrjones',help='(Optional) The tool to use to perform gene selection in the screening step of penalized eQTL-MVMR. The options are "gscreen" for the G-Screen method and "mrjones" for the MR-Jones method. The default is "mrjones"')
### post-estimation stuff
parser.add_argument('--adjustSEsForInflation',action='store',type=str,default='yes',help='(Optional) One of True or False. Should the standard errors of causal effect estimates be adjusted for inflation due to misspecified LD structure? The default is "yes".')
parser.add_argument('--saveRawData',action='store',type=str,default='false',help='(Optional) One of "true" "yes", "false" or "no" . Should the raw data used in multivariable MR be saved? This includes association estimates with each gene in a group, their standard errors, correspondingly the same for the outcome phenotype, and the estimated LD matrix. NOTE!!! This really only works if you are exploring a single candidate gene/locus')
parser.add_argument('--whereSaveRawData',action='store',type=str,default='',help='(Optional) If you gave --saveRawData "yes" or "true", this should be the directory/folder in which you want the raw data saved. The files will have the names "bybx.txt" (where -by- is in the first column),"UU.txt", and "LD.txt"')
parser.add_argument('--out',action='store',type=str,default='results',help='(Optional) filename without extension of location in which results should be written, in tab-separated .txt format. The default is results. Note that a tab-separated file named diagnostics.txt will automatically be written. This file contains information about the performance of multivariable MR.')
parser.add_argument('--iterativelySave',action='store',type=str,default='true',help='(Optional) Should results be saved iteratively as each chromosome is completed? This is helpful if you anticipate the analysis may take a relatively long time and you do not want to lose progress in case your access to the machine it is running on expires. The default is True.')
parser.add_argument('--cleanResults',action='store',type=str,default='yes',help='(Optional) Should the minimum interested results be saved? Put "no" or "false" if you want the complete results, including causal estiamtes from other MR methods. The default is "yes", to onnly output the most important results')
### creating networks
parser.add_argument('--graphNetworks',action='store',type=str,default='yes',help='(Optional) Should graph images of causal networks be stored? If "yes", they will go into the "plots" subfolder. The default is "yes".')
parser.add_argument('--graphNetworksFilePrefix',action='store',type=str,default='',help='(Optional) The default prefix for the filenames of graphed regulatory networks. Note that all network graph files will contain the lead gene in the filename automatically. What you put here will come before that. The default is "".')
parser.add_argument('--networkR2Thres',action='store',type=float,default=0.25,help='(Optional) Only loci with at least this much genetic variance in the disease phenotype explained will be considered for graph analysis using MNet. The default is 0.25.')
parser.add_argument('--networkDiseaseLabel',action='store',type=str,default='disease',help='(Optional) In graphs of causal networks, the disease node will have this label. This default is "disease"')
parser.add_argument('--networkGraphsOut',action='store',type=str,default='plots',help='(Optional) The directory/folder in which you would like plots of causal network graphs to be saved. The default is the "plots/" subdirectory/subfolder of the "hornet/" subdirectory/subfolder.')
### flags related to what is printed or not
parser.add_argument('--silence',action='store',type=str,default='no',help='(Optional) Should warnings about the size of the CHP window outside of the target locus be ignored? Put "true" or "yes". The default is "no".')
parser.add_argument('--hideVersion',action='store',type=str,default='yes',help='(Optional) Should the version number of HORNET be printed at the beginning? If not, put "false" or "no". The default is "yes".')
### done
args=parser.parse_args()
dontHide=(args.hideVersion.lower()=='yes') | (args.hideVersion.lower()=='true')
if dontHide:
    f=open('hornet.txt','r')
    con=f.read()
    print(con)
    f.close()

print(' ')
print('HORNET started '+time.ctime())

### first checking if all of the files/direcetories they gave actually exist
fileChecker(os.path.abspath(args.eQTLGWAS), 'eQTL GWAS directory')
fileChecker(os.path.abspath(args.phenoGWAS), 'phenotype GWAS filepath')
fileChecker(os.path.abspath(args.LDRef+'.bed'), 'LD reference panel filepath (w/o extension)')
fileChecker(os.path.abspath(args.writableDir), 'temporing working directory')
if args.whereSaveRawData!='':
    fileChecker(os.path.abspath(args.whereSaveRawData),'saved raw data directory')
od=(args.out.split('/'))[:-1]
fileChecker('/'.join(od), 'output/results directory')
# if no files in eQTL GWAS directory, there's a problem
if (os.listdir(args.eQTLGWAS))==0:
    raise ValueError('It looks like there are no files in the eQTL GWAS directory ({})'.format(args.eQTLGWAS))
# output directory for graphs
# graphoutdir=args.networkGraphsOut+'/' if args.networkGraphsOut!='plots' else os.getcwd()+'/plots/'
if os.path.isdir(os.path.abspath(args.networkGraphsOut))==False:
    raise ValueError('\n It looks like the directory you want to save network graph plots to does not exist({})'.format(args.networkGraphsOut))
# check other flags
if (args.geneScreener in ['gscreen','mrjones'])==False:
    raise ValueError('the "--geneScreener" flag accepts either "gscreen" or "mrjones", but you gave it {}'.format(args.geneScreener.lower()))

##
candidateGenesIsFile=os.path.isfile(args.candidateGenes)
if candidateGenesIsFile==False:
    candidateGenes=args.candidateGenes.split(',') if len(args.candidateGenes)>0 else []
else:
    candidateGenes=pandas.read_table(args.candidateGenes,header=None).iloc[:,0].values.tolist()

# drop Ensembl version type if it is present (if it is not, this code will have no effect)
if len(candidateGenes)>0:
    candidateGenes=[candidateGenes[i].split('.')[0] for i in range(0,len(candidateGenes))]

snplabs=args.snpLabels.split(',')
zlabs=args.zLabels.split(',')
ealabs=args.effectAlleles.split(',')
startwd=os.getcwd() # starting working directory - will be the absolute path of the HORNET software
impute=(args.imputeMissingEQTLs.lower()=='true') | (args.imputeMissingEQTLs.lower()=='yes')

### exposure parameters
dirGene=os.path.abspath(args.eQTLGWAS)
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
fpPheno=os.path.abspath(args.phenoGWAS)
writableDir=os.path.abspath(args.writableDir)
ldRefDir=os.path.abspath(args.LDRef)
# mapwd is fixed - users must use one of the maps I created
mapwd=os.path.abspath('data/maps/1kgPhase3maps')
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
dataPheno[dataPheno['isOutcomeClump']==True].to_csv(os.path.abspath('results/outcomeloci.csv'))

### loading key/dictionary/lookup/reference data
bim=pandas.read_csv(ldRefDir+'.bim',sep='\t',names=['chr','rsid','x','bp','a1','a2']) # for figuring out which chromosome is being used later
genekey=pandas.read_csv(os.getcwd()+'/data/biomart_gene_legend.csv.gz',engine='python') # to be used only if candidateGenes are specified

### running
# os.chdir(dirGene)
runningres=pandas.DataFrame() # to be filled in
runningdiagnostics=pandas.DataFrame() # to be filled in
infls=list()
edgeLists={}
usedCandidateGenes=[]
for _ in range(0, len(os.listdir(dirGene))):
    t0=time.perf_counter() # start timer
    # make sure I haven't already used all candidate genes
    mask=[candidateGenes[i] in usedCandidateGenes for i in range(0,len(candidateGenes))]
    if all(mask) & (len(candidateGenes)>0): # if all candidate genes have already been used
        print('All available candidate genes have been tested')
        break # don't need to consider additional chromosomes
    fpGene=os.listdir(dirGene)[_] # define filepath for exposure data
    dataGene=loadExposureGWASData(dirGene+'/'+fpGene,effectAlleleGene,zGene,rsidGene,snpBPGene,geneLabelGene,geneBPLabelGene,nLabel=nLabel,isRawGzippedGTEx=isRawGTEx,mapwd=mapwd) # load genes
    merged=mergeExposureAndOutcomeGWAS(dataGene,dataPheno) # merge with outcome/phenotype data
    mm=pandas.merge(bim,merged['geneSNP'].drop_duplicates(),left_on='rsid',right_on='geneSNP') # find chromosome
    chromosome=mm['chr'].values[0]
    del dataGene; # delete data we no longer need
    # check merged data for candidate gene(s) if user provided any (not much speed advantage checking merged vs dataGene bc merging is very fast)
    merged['Gene']=merged['Gene'].apply(lambda x: x.split('.')[0])
    if len(candidateGenes)>0:
        ch=merged.copy()
        ch=ch[ch['Gene'].isin(candidateGenes)] # candidateGenes already got splitted by '.'
        usedCandidateGenes.append(ch['Gene'].unique().tolist()) # keeping track of which candidate genes I've used so I know when to stop
        usedCandidateGenes=flatten_list(usedCandidateGenes)
        if ch.shape[0]==0:
            print('skipping chromsome '+str(chromosome)+' because there are no candidate genes on this chromosome')
            del ch
            continue
    
    print('Starting chromosome '+str(chromosome))
    if (analysisInPhenotypeLoci==True) & (len(candidateGenes)==0): # form gene groups
        geneGroups,ggKeys,lens,usedGenes=defGeneGroupsByOutcome(q0geneGroups, merged, KbWindow=outcomeClumpingKBWindow,closestK=numIndexGenesToFormGroup)
    else:
        geneGroups,ggKeys,lens,usedGenes=defGeneGroups(q0geneGroups,merged)
    # [print(geneGroupFinder(geneGroups,candidateGenes[i],isGtex=True if isRawGTEx else False)) for i in range(0,len(candidateGenes))]
    thingsMonitored,outerDict,edgeD=MVMRworkhorse(merged,geneGroups,ggKeys,writableDir=writableDir,ldRefDir=ldRefDir,isGtex=isRawGTEx,
                                                  analysisOnlyInOutcomeLoci=analysisInPhenotypeLoci,outcomeLociMbWindow=outcomeLociMbWindow,
                                                  ldUpperLimit=ldUpperLimit, ldOtherLociOtherPt=ldOtherLociOtherPt,ldOtherLociR2=ldOtherLociR2,
                                                  ldOtherLociWindow=ldOtherLociWindow, q0Correls=q0Correls,nMinCorrels=nMinCorrels,jointChiGenesP=jointChiGenesP,
                                                  opAlpha=opAlpha,nMinIVs=nMinIVs,
                                                  hessMinScale=hessMinScale,silence=silence, UniMRIVPThreshold=UniMRIVPThreshold, candidateGenes=candidateGenes,
                                                  assumeNoSampleOverlap=assumeNoSampleOverlap,shrinkBiasCorrection=shrinkBiasCorrection,networkR2Thres=args.networkR2Thres,impute=impute,saveData=saveData,
                                                  geneSelector=args.geneScreener.lower())
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
    # if user wanted to save raw data, do that now inside of the chosen directory
    if saveData: 
        fpout=os.path.abspath(args.whereSaveRawData)+'/'
        print('raw data are saved in the "bybx.txt", "UU.txt", and "LD.txt" files in the {} directory'.format(fpout))
        ln=list(outerDict)[0]
        numpy.savetxt(fpout+'bybx.txt',numpy.column_stack((outerDict[ln]['by'],outerDict[ln]['bX'])))
        numpy.savetxt(fpout+'UU.txt',outerDict[ln]['UU'])
        numpy.savetxt(fpout+'LD.txt',outerDict[ln]['regLD'])
        numpy.savetxt(fpout+'genes.txt',outerDict[ln]['data']['Gene'])
        numpy.savetxt(fpout+'snps.txt',outerDict[ln]['data']['IVs'])
        numpy.savetxt(fpout+'geneBPs.txt',outerDict[ln]['data']['GeneBP'])

    # LDLeadGeneSNPs(res,thischr,writableDir,ldRefDir) # write out results to be read in by an R program that will make figures etc.    
    t1=time.perf_counter()-t0
    print('  Chromosome '+str(chromosome)+' complete '+str(round(t1/60,1))+' minutes later')
    # iteratively save so I don't lose any progress between genes
    if (args.iterativelySave.lower()=='true') | (args.iterativelySave.lower()=='yes'):
        runningres.to_csv(args.out+'_tempresults.txt',sep='\t')
        runningdiagnostics.to_csv(args.out+'_diagnostics.txt',sep='\t')

# if no chromosomes could be analyzed, there is nothing to save - tell the user and exit
if runningres.shape[0]==0:
    raise ValueError('\n\n No genes could be analyzed in these data.\n Consider adjusting your criteria for forming \ngene groups, IV sets, etc., or the data itself\nas the cause of this result\n\n')

fp1=os.path.abspath(args.out)+'_results.txt'
fp2=os.path.abspath(args.out)+'_diagnostics.txt'
# did the user want us to clean up the results?
cleanres=(args.cleanResults.lower()=='yes') | (args.cleanResults.lower()=='true')
if cleanres:
    cr=runningres.copy()
    cr['Pratt']=cr['MRBEEPostSelec_MVMR_Est']*cr['MRBEE_UVMR_Est']
    cr=cr[['Gene','geneBP','Chromosome','GeneSelected','MRBEEPostSelec_MVMR_Est','MRBEEPostSelec_MVMR_SE','LocusR2','Pratt','CHRspecificGroupID']]
    cr=cr.rename(columns={'MRBEEPostSelec_MVMR_Est':'Est','MRBEEPostSelec_MVMR_SE':'SE'})
    cr.to_csv(fp1,sep='\t')
else:
    runningres.to_csv(fp1,sep='\t')

runningdiagnostics.to_csv(fp2,sep='\t')
# also save a copy of runningres to HORNET/res.csv so it can be read by plotres.r
copyfpout=os.path.abspath(os.getcwd())
runningres.to_csv(copyfpout+'/res.csv') # full (not super cleaned) data will still be written
print('Results are written to '+fp1)
print('Analysis diagnostics are written to '+fp2)
# save executing commands of plotres.r for later - don't want to cause an early error bc users don't have R installed

# delete iteratively saven files if user chose to iteratively save results
#if platform.system()!='Windows':
#    op=os.path.abspath(args.out+'_tempresults.txt')
#    if (args.iterativelySave.lower()=='true') | (args.iterativelySave.lower()=='yes'):
#        delcall=callDelete()
#        out=subprocess.call([delcall, op],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
#        awd=os.path.abspath(writableDir)
#        out=subprocess.call([delcall, awd+'/tempOut.ld'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
#        out=subprocess.call([delcall, awd+'/tempOut.log'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
#        out=subprocess.call([delcall, awd+"/myExtract.txt"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
#        out=subprocess.call([delcall, awd+"/outcomeplinkout.clumped"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
#        out=subprocess.call([delcall, awd+"/outcomeplinkout.log"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
#        out=subprocess.call([delcall, awd+"/outcomePsout.txt"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

### network construction (done in Python)
# first find top args.networksInTopKLoci loci
graphNetworks=(args.graphNetworks.lower()=='yes') | (args.graphNetworks.lower()=='true')
hasSaved=False
if (runningres.shape[1]>1) & (graphNetworks):
    vals=numpy.sort(runningres['LocusR2'].unique())[::-1]
    for _ in range(0,len(vals)):
        if vals[_]<args.networkR2Thres:
            continue
        dcut=runningres[runningres['LocusR2']==vals[_]]
        genes=dcut['Gene'].values.tolist(); genes=[genes[_].split(',')[0] for _ in range(0,len(genes))]
        if dcut.shape[0]<1:
            continue
        chrom,grp=int(dcut['Chromosome'][0]),dcut['CHRspecificGroupID'][0]
        edges=edgeLists['Chromosome'+str(chrom)][grp]
        if sum(sum(edges))==0:
            continue
        G=addEdges(edges)
        fl=flatten_list([args.networkDiseaseLabel,genes])
        G=nx.relabel_nodes(G,mapping=mapNodeNames(edges,fl)) # add new labels
        # color_map=colorCoreGenes(edges) # add colors for core and peripheral genes
        fig=matplotlib.pyplot.figure(figsize=(12,8))
        warnings.filterwarnings("ignore")
        nx.draw(G, with_labels=True,ax=fig.add_subplot(),alpha=0.5,font_size=10,node_color='skyblue') # all same color - less likely to produce an error
        leadgene=dcut.loc[abs(dcut['MRJonesEst']).idxmax()]['Gene']
        leadgene=leadgene.split('.')[0]
        fig.savefig(os.path.abspath(args.networkGraphsOut)+'/'+args.graphNetworksFilePrefix+'_'+leadgene+'_graph.png',dpi=500)
        hasSaved=True

print(' ')
if graphNetworks & hasSaved:
    print('Network graphs are saved in the '+os.path.abspath(graphoutdir)+' folder')

# finally, if user has R installed, execute commands in HORNET/plotres.r
# from subprocess import Popen, PIPE
# proc = Popen(["which", "R"],stdout=PIPE,stderr=PIPE)
# exit_code = proc.wait()
# if exit_code == 0:
#     cmd=['Rscript','plotres.r']
#     subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

print(' ')
print('HORNET ended '+time.ctime())

