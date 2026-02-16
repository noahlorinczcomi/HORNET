# standard imports
import pandas as pd
import numpy as np
import gzip
import os

# custom imports
from cli import parse_args
from utils import flatten_list


def main():
    args = parse_args()
    return None


if __name__ == "__main__":
    args = parse_args()
    main(args)

printres = (args.printResults.lower() == "yes") | (args.printResults.lower() == "true")
# load hscores.txt.gz data
file = gzip.open("data/gtexv8_tissues_h2scores.txt.gz", "rb")
hscores = pd.read_table(file, sep=None, engine="python")
hscores["Gene"] = [
    (hscores["Gene"][_].split("."))[0] for _ in range(0, hscores.shape[0])
]
isfile = os.path.exists(
    args.candidateGenes
)  # candidateGenes='ENSG00000130203.9,ENSG00000087237.10,ENSG00000169174.10': APOE,CETP,PCSK9
dosave = (args.saveAsFile.lower() == "true") | (args.saveAsFile.lower() == "yes")
if isfile:
    cg = pd.read_csv(args.candidateGenes, engine="python", header=None).values.tolist()
    cg = flatten_list(cg)
else:
    cg = args.candidateGenes.split(",")

cg = [
    (cg[_].split("."))[0].upper() for _ in range(0, len(cg))
]  # if ensembl IDs.versionCTRL were given in candidateGenes, drop version CTRL

# load key file of <Gene name>, <Ensembl ID> pairs
file = gzip.open("data/biomart_gene_legend.csv.gz", "rb")
genekey = pd.read_table(file, sep=None, engine="python")[["EnsemblID", "GeneName"]]
# try to determine if user gave me gene names or Ensembl IDs (don't have to consider Ensembl version control from here on out)
match_ens = genekey[genekey["EnsemblID"].isin(cg)].shape[0]
match_name = genekey[genekey["GeneName"].isin(cg)].shape[0]
if (match_ens == 0) & (match_name == 0):
    raise ValueError(
        "\n\n\n None of these genes are found in the data. Are you sure you specified their IDs correctly? You can use either Ensembl IDs or gene names. If they are correct, then none of these genes have eQTLs in any of these 49 tissues in GTEx v8 \n"
    )

# if they are names, use names from here on out
if match_name > match_ens:  # if the user gave gene names instead of Ensembl IDs
    hscores = pd.merge(
        hscores, genekey, how="inner", left_on="Gene", right_on="EnsemblID"
    )  # merge hscores with gene names key
    hscores = hscores.rename(
        columns={"Gene": "EnsemblID2"}
    )  # rename 'Gene' column, which is Ensembl ID, in hscores to something else and make new 'Gene' column gene name
    hscores = hscores.rename(columns={"GeneName": "Gene"})

# find genes that user gave that may not be in the data
foundgenes = hscores[hscores["Gene"].isin(cg)]["Gene"].unique()
mask = [cg[_] in foundgenes for _ in range(0, len(cg))]
mask = np.array(mask)
if (sum(mask == False) > 0) & printres:
    notfound = np.array(cg)[mask == False].tolist()
    mess = ", ".join(notfound)
    print("\n Warning! The following genes were not found in the GTEx v8 data: " + mess)

# collecting results
t5 = pd.DataFrame()
for _ in range(0, len(cg)):
    toadddf = hscores[hscores["Gene"] == cg[_]].sort_values(
        by="h2score", ascending=False
    )
    toadd = toadddf.iloc[:5, :]
    if toadd.shape[0] == 0:
        continue
    t5 = t5.append(toadd.copy(), ignore_index=True)

# summarise the results in some way
tissues = t5["Tissue"].unique()
t5.groupby(["Tissue"]).agg({"h2score": "max", "Tissue": "count"})

t5g = t5.groupby(["Tissue"]).agg(
    Maxh2score=pd.NamedAgg(column="h2score", aggfunc=max),
    nSignifSNPs=pd.NamedAgg(column="nSignifSNPs", aggfunc=sum),
    TissueCount=pd.NamedAgg(column="h2score", aggfunc="count"),
    Genes=pd.NamedAgg(column="Gene", aggfunc=lambda x: ",".join([str(_) for _ in x])),
)
t5g = t5g.sort_values(by=["TissueCount", "nSignifSNPs", "Maxh2score"], ascending=False)

if dosave:  # save a new file with all results for each gene
    t5.to_csv(args.outFilename + ".txt")
    print("full results written to" + args.outFilename + ".txt")

if printres:
    print("\n The top 10 tissues for these genes are the following: \n")
    print(t5g.iloc[: args.nTopTissues, :])
    print("\n")
