import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas
import argparse
from scipy.stats import norm
import numpy
import os

parser=argparse.ArgumentParser(prog='Summarise HORNET Results',description='This program subsets the full results to only significant genes, where you defined "significant"')

### data
parser.add_argument('-r','--results',action='store',type=str,help='(Required) Full filepath (w/ extension) to results outputted by hornet.py')
parser.add_argument('-pr','--prattThreshold',action='store',type=float,default=0.1,help='(Optional) Lower threshold for Pratt index values. The default is 0.1.')
parser.add_argument('-v','--localVarianceExplainedThreshold',action='store',type=float,default=0.5,help='(Optional) The locus in which a prioritized gene sits must have explained at least this much local genetic variance in the disease phenotype. The default is 0.5.')
parser.add_argument('-pv','--pValueThreshold',action='store',type=float,default=5e-5,help='(Optional) Only genes with P<this threshold will be considered as priorities. The default is 5e-5.')
parser.add_argument('-tk','--topKGenes',action='store',type=int,default=20,help='(Optional) Only the top --topKGenes will be printed (sorted by Pratt index). The default is 20.')
parser.add_argument('--print',action='store',type=str,default='yes',help='(Optional) Should the results be printed in addition to being saved? The default is yes.')
parser.add_argument('-o','--out',action='store',type=str,default='',help='(Optional) Filename (without extension, with possible prefix directory) to which results should be written. The default is "", which means no results are saved.')

args=parser.parse_args()

# load data
data=pandas.read_csv(args.results,engine='python',sep=None)
### varios subsetting
data=data[['Gene','CHRspecificGroupID','geneBP','Chromosome','MRBEEPostSelec_MVMR_Est','MRBEEPostSelec_MVMR_SE','MRBEE_UVMR_Est','RsquaredMRJones']].dropna()
data['Pratt']=data['MRBEEPostSelec_MVMR_Est']*data['MRBEE_UVMR_Est']
data=data[(data['Pratt']>0) & (data['Pratt']<=1)]
sub1=data['Pratt']>args.prattThreshold
sub2=data['RsquaredMRJones']>args.localVarianceExplainedThreshold
z=data['MRBEEPostSelec_MVMR_Est']/data['MRBEEPostSelec_MVMR_SE']
data['CausalPvalue']=2*norm.cdf(-abs(z),0,1)
sub3=abs(z)>norm.ppf(1-args.pValueThreshold)
subs=numpy.column_stack((sub1,sub2,sub3))
keepers=subs.sum(axis=1)==3
if sum(keepers)==0:
    raise ValueError('\n\n No genes met the criteria you specified. \n Consider changing the parameters, which you can check using: python sumres.py -h')
cdata=data[keepers].copy()
# top gene in each locus (largest effect size)
cdata['cs']=cdata['Chromosome'].astype(str)+cdata['CHRspecificGroupID']
topg=cdata.groupby('cs').agg({'Pratt': 'max'})
topg=topg.sort_values('Pratt',ascending=False)
topg=topg.iloc[:args.topKGenes,:]
pgenes=pandas.merge(topg,cdata,how='left',left_on='Pratt',right_on='Pratt') # note I am merging on Pratt - should work, unlikely to get 6-9 digits of the same by chance
pgenes['GeneMarkername']=pgenes['Chromosome'].astype(str)+':'+pgenes['geneBP'].astype(str)
pgenes=pgenes[['Gene','GeneMarkername','Pratt','RsquaredMRJones','CausalPvalue']]
printbool=(args.print.lower()=='yes') | (args.print.lower()=='true')
if printbool:
    print(pgenes)

# save
if args.out!='':
    pgenes.drop_duplicates().to_csv(args.out+'.csv')
