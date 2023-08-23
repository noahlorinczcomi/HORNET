
#!/usr/bin/env python
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas
from scipy.stats import norm
import argparse
import subprocess
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

parser=argparse.ArgumentParser(prog='',description='This program will search your results.txt file to find novel loci; \n i.e., loci that were not detected at a standard significance threshold in the original disease GWAS but were detected by HORNET. \n Negative confounding between genes is one explanation for why this may happen.')

### data
parser.add_argument('-m','--message',action='store',type=str,default='yes',help='(Optional) path to your results file.')
parser.add_argument('-r','--results',action='store',type=str,default='results.txt',help='(Required) path to your results file.')
parser.add_argument('-pg','--phenoGWAS',action='store',type=str,help='(Required) path to your results file.')
parser.add_argument('-ld','--LDRef',action='store',type=str,default=os.path.abspath('data/ldref/1kg3TRANS'),help='(Optional) path to your results file.')
parser.add_argument('-rsid','--phenoRSID',action='store',type=str,default='SNP',help='(Required) path to your results file.')
parser.add_argument('-pz','--phenoZ',action='store',type=str,default='Z',help='(Required if no P-value column in phenotype GWAS) path to your results file.')
parser.add_argument('-pp','--phenoP',action='store',type=str,default='missing',help='(Required) path to your results file.')
parser.add_argument('-psig','--phenoPSignifThres',action='store',type=float,default=5e-8,help='(Optional) path to your results file.')
parser.add_argument('-pr2','--phenoR2Thres',action='store',type=float,default=0.1,help='(Optional) path to your results file.')
parser.add_argument('-pkb','--phenoKbThres',action='store',type=int,default=500,help='(Optional) path to your results file.')
parser.add_argument('-n','--novelKbWindow',action='store',type=float,default=1,help='(Optional) path to your results file.')
parser.add_argument('-o','--out',action='store',type=str,default='newloci.csv',help='(Optional) path to your results file.')
args=parser.parse_args()
if args.message.lower() in ['yes','true']:
    print('\n\n NOTE! This is only meaningful if you used HORNET to analyze loci that \n were "less significant" than what you are calling a "significant" outcome locus here\n\n you can silence this message by adding the "-m no" flag.\n\n')
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
print(type(oc))
print(oc.shape)
oc=oc.dropna().values # for some reason some NAs get appended to the end
print('no errors')
#print(oc.head())
res=pandas.read_csv(os.path.abspath(args.results),sep=None,engine='python')
reschrs=res['Chromosome'].unique()
for _ in range(0,)
