#!/usr/bin/env python
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas
import numpy
import gzip
import argparse
import os

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

parser=argparse.ArgumentParser(prog='This tool can help guide you in choosing which tissues are good choices for performing genome-wide eQTL-MVMR with MR-Jones',
                    description='This program relies on pre-calculated heritability scores, which are proportional to SNP heritability, from cis-SNPs and GTEx v8 summary data to create an ordered list of tissues in which your candidate genes have the strongest eQTL signals.')

parser.add_argument('--candidateGenes',action='store',type=str,help='(Required) This can either be a comma-separated list of gene IDs or the path to a file (with extension) that contains one gene ID per row and no header. NOTE, all gene IDs must follow the Ensembl gene ID format (e.g., "ENSG..."), but they do not need the additional ".<xx>" identifier. For example, "ENSG00000227232" is sufficient, you do not need to use "ENSG00000227232.5"')
parser.add_argument('--nTopTissues',action='store',type=int,default=10,help='(Required) The number of top tissues to display. The default is 10. Note, if you put "yes" for --saveAsFile, all results will be written out, not just these that are printed.')
parser.add_argument('--saveAsFile',action='store',type=str,default='no',help='(Optional) If you would like the results to be saved to a file, put "true" or "yes", otherwise ignore this flag or put "false" or "no"')
parser.add_argument('--outFilename',action='store',type=str,default='tissue_results.txt',help='(Optional) If you put "true" or "yes" for the --saveAsFile flag, you can specify the filename (without extension) of the output file here. The default is "tissue_results.txt"')
parser.add_argument('--printResults',action='store',type=str,default='yes',help='(Optional) Please give this flag "no" if you do not want the results to be printed to the console. Note that if you give this a "no", you will want to make sure the results are saved in a file, otherwise you will never see the results. The default is to print results to the console.')


args=parser.parse_args()
printres=(args.printResults.lower()=='yes') | (args.printResults.lower()=='true')
# load hscores.txt.gz data
file=gzip.open('data/gtexv8_tissues_h2scores.txt.gz', 'rb')
hscores=pandas.read_table(file,sep=None,engine='python')
hscores['Gene']=[(hscores['Gene'][_].split('.'))[0] for _ in range(0, hscores.shape[0])]
isfile=os.path.exists(args.candidateGenes) # candidateGenes='ENSG00000130203.9,ENSG00000087237.10,ENSG00000169174.10': APOE,CETP,PCSK9
dosave=(args.saveAsFile.lower()=='true') | (args.saveAsFile.lower()=='yes')
if isfile: 
    cg=pandas.read_csv(args.candidateGenes,engine='python',header=None).values.tolist()
    cg=flatten_list(cg)
else: 
    cg=args.candidateGenes.split(',') 

cg=[(cg[_].split('.'))[0].upper() for _ in range(0,len(cg))] # if ensembl IDs.versionCTRL were given in candidateGenes, drop version CTRL

# load key file of <Gene name>, <Ensembl ID> pairs
file=gzip.open('data/biomart_gene_legend.csv.gz', 'rb')
genekey=pandas.read_table(file,sep=None,engine='python')[['EnsemblID','GeneName']]
# try to determine if user gave me gene names or Ensembl IDs (don't have to consider Ensembl version control from here on out)
match_ens=genekey[genekey['EnsemblID'].isin(cg)].shape[0]
match_name=genekey[genekey['GeneName'].isin(cg)].shape[0]
if (match_ens==0) & (match_name==0):
    raise ValueError('\n\n\n None of these genes are found in the data. Are you sure you specified their IDs correctly? You can use either Ensembl IDs or gene names. If they are correct, then none of these genes have eQTLs in any of these 49 tissues in GTEx v8 \n')

# if they are names, use names from here on out
if match_name>match_ens: # if the user gave gene names instead of Ensembl IDs
    hscores=pandas.merge(hscores,genekey,how='inner',left_on='Gene',right_on='EnsemblID') # merge hscores with gene names key
    hscores=hscores.rename(columns={'Gene': 'EnsemblID2'}) # rename 'Gene' column, which is Ensembl ID, in hscores to something else and make new 'Gene' column gene name 
    hscores=hscores.rename(columns={'GeneName': 'Gene'})

# find genes that user gave that may not be in the data
foundgenes=hscores[hscores['Gene'].isin(cg)]['Gene'].unique()
mask=[cg[_] in foundgenes for _ in range(0,len(cg))]; mask=numpy.array(mask)
if (sum(mask==False)>0) & printres:
       notfound=numpy.array(cg)[mask==False].tolist()
       mess=', '.join(notfound)
       print('\n Warning! The following genes were not found in the GTEx v8 data: '+mess)

# collecting results
t5=pandas.DataFrame()
for _ in range(0,len(cg)):
    toadddf=hscores[hscores['Gene']==cg[_]].sort_values(by='h2score',ascending=False)
    toadd=toadddf.iloc[:5,:]
    if toadd.shape[0]==0:
        continue
    t5=t5.append(toadd.copy(),ignore_index=True)

# summarise the results in some way
tissues=t5['Tissue'].unique()
t5.groupby(['Tissue']).agg({'h2score': 'max','Tissue': 'count'})

t5g=t5.groupby(['Tissue']).agg(
    Maxh2score=pandas.NamedAgg(column='h2score',aggfunc=max),
    nSignifSNPs=pandas.NamedAgg(column='nSignifSNPs',aggfunc=sum),
    TissueCount=pandas.NamedAgg(column='h2score',aggfunc='count'),
    Genes=pandas.NamedAgg(column='Gene',aggfunc=lambda x: ','.join([str(_) for _ in x]))
)
t5g=t5g.sort_values(by=['TissueCount','nSignifSNPs','Maxh2score'],ascending=False)

if dosave: # save a new file with all results for each gene
    t5.to_csv(args.outFilename+'.txt')
    print('full results written to'+args.outFilename+'.txt')

if printres:
    print('\n The top 10 tissues for these genes are the following: \n')
    print(t5g.iloc[:args.nTopTissues,:])
    print('\n')
