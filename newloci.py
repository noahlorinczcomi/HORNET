
#!/usr/bin/env python
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas
from scipy.stats import norm
import argparse
import subprocess
import numpy
import platform
def callPlink():
    sys=platform.system()
    if sys=='Linux':
        call_='./plinkdir/linux/plink'
    elif sys=='Darwin':
        call_='./plinkdir/mac/plink'
    elif sys=='Windows':
        call_='plinkdir/windows/plink.exe'
    return call_

def callDelete():
    sys=platform.system()
    if (sys=='Linux') | (sys=='Darwin'):
        call_='rm'
    elif sys=='Windows':
        call_='del'
    return call_

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
parser=argparse.ArgumentParser(prog='',description='This program will search your results.txt file to find novel loci; \n i.e., loci that were not detected at a standard significance threshold in the original disease GWAS but were detected by HORNET. \n Negative confounding between genes is one explanation for why this may happen.')

### data
parser.add_argument('-m','--message',action='store',type=str,default='yes',help='(Optional) should the warning message be printed? Default is "yes".')
parser.add_argument('-r','--results',action='store',type=str,default='results.txt',help='(Required) filepath to your results file (with file extension).')
parser.add_argument('-pg','--phenoGWAS',action='store',type=str,help='(Required) path to your phenotype GWAS file of summary statistics (with file extension).')
parser.add_argument('-ld','--LDRef',action='store',type=str,default=os.path.abspath('data/ldref/1kg3TRANS'),help='(Optional) path to the LD reference panel file (without extension). the default is"data/ldref/1kg3TRANS".')
parser.add_argument('-rsid','--phenoRSID',action='store',type=str,default='SNP',help='(Required) column name in phenotype GWAS data indicating rsID.')
parser.add_argument('-pz','--phenoZ',action='store',type=str,default='Z',help='(Required if no P-value column in phenotype GWAS data) column name in phenotype GWAS indicating Z-statistics for associations between SNPs and phenotype.')
parser.add_argument('-pp','--phenoP',action='store',type=str,default='missing',help='(Required if no Z-statistic column in phenotype GWAS data) column name indicating P-value for association in the phenotype GWAS data.')
parser.add_argument('-psig','--phenoPSignifThres',action='store',type=float,default=5e-8,help='(Optional) P-value threshold for classifying an SNP as "significantly" associated with the outcome/phenotype. The default is 5e-8.')
parser.add_argument('-pr2','--phenoR2Thres',action='store',type=float,default=0.1,help='(Optional) SNPs in LD beyond this threshold with a lead outcome/phenotype SNP will be considered part of the same locus as the lead outcome SNP. The defaul is 0.1.')
parser.add_argument('-pkb','--phenoKbThres',action='store',type=int,default=500,help='(Optional) Size of window in which PLINK should clump SNPs to identify idependent outcome signals. The default is 500, which means 500Kb.')
parser.add_argument('-hpv','--hornetSignifPThres',action='store',type=float,default=5e-5,help='(Optional) Genes with a P-value from HORNET less than this threshold will be considered detected/statistically significant. The default is 5e-5.')
parser.add_argument('-hr','--hornetRSquaredThres',action='store',type=float,default=0.10,help='(Optional) Genes in a locus that explains at least this much outcome/phenotype variance will be considered detected/statistically significant. The default is 0.1.')
parser.add_argument('-hpr','--hornetPrattThres',action='store',type=float,default=0.05,help='(Optional) Genes with a Pratt index value (gene-specific outcome variance attributable) from HORNET greater than this threshold will be considered detected/statistically significant. The default is 0.05.')
parser.add_argument('-n','--novelKbWindow',action='store',type=float,default=1,help='(Optional) Genes detected in HORNET must be at least this many Kb away from any known outcome signals to be considered novel. The default is 500, which means 500Kb.')
parser.add_argument('-o','--out',action='store',type=str,default='newgenes',help='(Optional) Name of file to which results should be written (without file extension; it will be .csv). The default is "newgenes".')
parser.add_argument('-p','--print',action='store',type=str,default='yes',help='(Optional) Should results also be printed to the console? The default is "yes".')
args=parser.parse_args()
if args.message.lower() in ['yes','true']:
    print('\n\n NOTE! Use of this is only informative if you used HORNET to analyze loci that \n were "less significant" than what you are calling a "significant" outcome locus here\n\n you can silence this message by adding the "-m no" flag.\n\n')
### steps
# 1) make sure PLINK can read the outcome GWAS
# 2) perform clumping
# 3) load results
# 4) determine if any results are within +- user-defined window
# 5) if not, novel! print/save those

# if no P-value column, will need to add it
if args.phenoP=='missing':
    spp=args.pz.split('@')
    datain=pandas.read_csv(os.path.abspath(args.phenoGWAS),sep=None,engine='python')
    if len(spp)>1:
        betacol=spp[0]
        secol=spp[1]
        # need to save out data with a Z-score or P-value column
        datain['Z']=datain[betacol]/datain[secol]
        datain['P']=2*norm.cdf(-abs(datain['Z']))
    else:
        datain['P']=2*norm.cdf(-abs(datain[args.phenoZ]))
    datain[[args.phenoRSID,'P']].to_csv(os.path.abspath(os.getcwd())+'/a.txt') # will be removed later    
else:
    cmd=[callPlink(),"--bfile", os.path.abspath(args.LDRef), "--clump", os.path.abspath(args.phenoGWAS), "--clump-kb", str(args.phenoKbThres), "--clump-p1", str(args.phenoPSignifThres),"--clump-p2", str(args.phenoPSignifThres), "--clump-r2", str(args.phenoR2Thres), "--out", os.path.abspath(os.getcwd())+'/a','--clump-snp-field',args.phenoRSID]
out=subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) # executing the command in the terminal
oc=pandas.read_fwf('a.clumped', sep='\t') # load clumps; will have columns SNP, CHR, BP, etc.
oc=oc.dropna() # for some reason some NAs get appended to the end
res=pandas.read_csv(os.path.abspath(args.results),sep=None,engine='python')
# did user give me complete results or cleaned results?
if 'Est' in res.columns.tolist(): # if this is in there, they gave me clean results
    res['Q']=(res['Est']/res['SE'])**2
else:
    res['Q']=(res['MRBEEPostSelec_MVMR_Est']/res['MRBEEPostSelec_MVMR_SE'])**2
    res['Pratt']=res['MRBEEPostSelec_MVMR_Est']*res['MRBEEPostSelec_MVMR_SE']
res=res[['Gene','Chromosome','geneBP','Q','RsquaredMRJones','Pratt']]
res=res[res['Q']>(norm.ppf(args.hornetSignifPThres)**2)]
if res.shape[0]==0:
    raise ValueError('\n\n no significant genes given your definition of "significant"')
res=res[res['RsquaredMRJones']>args.hornetRSquaredThres]
if res.shape[0]==0:
    raise ValueError('\n\n no significant genes given your definition of "significant"')
res=res[res['Pratt']>args.hornetPrattThres]
if res.shape[0]==0:
    raise ValueError('\n\n no significant genes given your definition of "significant"')
chrs=res['Chromosome'].unique()
novelGenes=[]
novelchrs=[]
closestSNPs=[]
closestBPs=[]
for _ in range(0,len(chrs)): # for each chromosome in the results
    novelGenes=flatten_list(novelGenes)
    novelchrs=flatten_list(novelchrs)
    closestSNPs=flatten_list(closestSNPs)
    closestBPs=flatten_list(closestBPs)
    # 
    reschrdf=res[res['Chromosome']==chrs[_]]
    clumpchrdf=oc[oc['CHR']==chrs[_]]
    if clumpchrdf.shape[0]==0:
        novelGenes.append(reschrdf['Gene'].unique().tolist())
        novelchrs.append((reschrdf['Chromosome'].values.tolist()))
        closestSNPs.append([float('nan')]*reschrdf.shape[0])
        closestBPs.append([float('nan')]*reschrdf.shape[0])
        continue
    for j in range(0,reschrdf.shape[0]):
        # is the significant gene in HORNET within args.novelKbWindow of any known signals?
        alldist=abs((reschrdf['geneBP'].values)[j]-clumpchrdf['BP'].values)
        # this code is actually recording all genes meeting this criteria, even if they are in the same locus as each other
        if all(alldist>(args.novelKbWindow*1000)):
            novelGenes.append((reschrdf['Gene'].values)[j])
            novelchrs.append((reschrdf['Chromosome'].values)[j])
            closeix=numpy.argmin(alldist)
            closestSNPs.append((clumpchrdf['SNP'].values)[closeix])
            closestBPs.append((clumpchrdf['BP'].values)[closeix])
novelGenes=numpy.array(flatten_list(novelGenes))
closestSNPs=numpy.array(flatten_list(closestSNPs))
closestBPs=numpy.array(flatten_list(closestBPs))
# be prepared to tell the user which outcome signals were the closest to these genes
pdf=pandas.DataFrame({'NovelGene':novelGenes,'Chromosome':novelchrs,'ClosestOutcomeSNP':closestSNPs,'ClosestOutcomeSNPBP':closestBPs})
# add other info like pratt index etc
res=res.rename(columns={'Gene':'NovelGene','RsquaredMRJones':'LocusR2','Q':'CausalChisq'})
pdf=pandas.merge(pdf,res[['NovelGene','geneBP','Pratt','CausalChisq','LocusR2']],how='left',left_on='NovelGene',right_on='NovelGene')
pdf['NovelGene']=pdf['NovelGene'].apply(lambda x: x.split('.')[0])
# load gene name key and attach actual gene names
gk=pandas.read_csv(os.path.abspath('data/biomart_gene_legend.csv.gz'),sep=None,engine='python')
pdf=pandas.merge(pdf,gk.rename(columns={'EnsemblID':'NovelGene'}),how='left',left_on='NovelGene',right_on='NovelGene')
pdf=pdf.rename(columns={'GeneName':'NovelGeneName'})
pdf=pdf[['NovelGene','NovelGeneName','Chromosome','geneBP','CausalChisq','Pratt','LocusR2','ClosestOutcomeSNP','ClosestOutcomeSNPBP']]
doprint=(args.print.lower()=='yes') | (args.print.lower()=='true')
if doprint:
    print(pdf)
# save out
fpout=os.path.abspath(args.out)+'.csv'
pdf.to_csv(fpout)
print('\n Novel genes are saved to {}'.format(fpout))
# remove stuff I wrote out

cmd=[callDelete(),os.path.abspath(os.getcwd())+'/a.log']
out=subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) # executing the command in the terminal
cmd=[callDelete(),os.path.abspath(os.getcwd())+'/a.clumped']
out=subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) # executing the command in the terminal
if args.phenoP=='missing':
    cmd=[callDelete(),os.path.abspath(os.getcwd())+'/a.txt']
    out=subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) # executing the command in the terminal
