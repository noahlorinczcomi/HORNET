# standard imports
import os
import pandas as pd
import subprocess
import numpy as np
import platform
from scipy.stats import norm

# custom imports
from cli import parse_args
import utils


def main(args):
    """
    steps
    1) make sure PLINK can read the outcome GWAS
    2) perform clumping
    3) load results
    4) determine if any results are within +- user-defined window
    5) if not, novel! print/save those
    """
    args = parse_args()

    # if no P-value column, will need to add it
    datain = pd.read_csv(os.path.abspath(args.phenoGWAS), sep=None, engine="python")
    if args.phenoP == "missing":
        spp = args.phenoZ.split("@")
        if len(spp) > 1:
            betacol = spp[0]
            secol = spp[1]
            # need to save out data with a Z-score or P-value column
            datain["Z"] = datain[betacol] / datain[secol]
            datain["P"] = 2 * norm.cdf(-abs(datain["Z"]))
        else:
            datain["P"] = 2 * norm.cdf(-abs(datain[args.phenoZ]))
        datain = datain.rename(columns={args.phenoRSID: "SNP"})
    else:
        datain = datain.rename(columns={args.phenoRSID: "SNP", args.phenoP: "P"})
    datain[["SNP", "P"]].to_csv("a.txt", sep="\t", index=False)
    cmd = [
        utils.callPlink(),
        "--bfile",
        os.path.abspath(args.LDRef),
        "--clump",
        os.path.abspath("a.txt"),
        "--clump-kb",
        str(args.phenoKbThres),
        "--clump-p1",
        str(args.phenoPSignifThres),
        "--clump-p2",
        str(args.phenoPSignifThres),
        "--clump-r2",
        str(args.phenoR2Thres),
        "--out",
        "a",
    ]
    _ = subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    print(*cmd)
    fpin = os.path.abspath("a.clumped")
    print(fpin)
    print(os.path.isfile(fpin))
    oc = pd.read_fwf(
        fpin, sep="\t"
    )  # load clumps; will have columns SNP, CHR, BP, etc.
    oc = oc.dropna()  # for some reason some NAs get appended to the end
    res = pd.read_csv(os.path.abspath(args.results), sep=None, engine="python")
    # did user give me complete results or cleaned results?
    if (
        "LocusR2" in res.columns.tolist()
    ):  # if this is in there, they gave me clean results
        res["Q"] = (res["Est"] / res["SE"]) ** 2
        res = res.rename(columns={"LocusR2": "RsquaredMRJones"})
    else:
        res["Q"] = (res["MRBEEPostSelec_MVMR_Est"] / res["MRBEEPostSelec_MVMR_SE"]) ** 2
        res["Pratt"] = res["MRBEEPostSelec_MVMR_Est"] * res["MRBEEPostSelec_MVMR_SE"]

    res = res[["Gene", "Chromosome", "geneBP", "Q", "RsquaredMRJones", "Pratt"]]
    res = res[res["Q"] > (norm.ppf(args.hornetSignifPThres) ** 2)]
    if res.shape[0] == 0:
        raise ValueError(
            '\n\n no significant genes given your definition of "significant"'
        )
    res = res[res["RsquaredMRJones"] > args.hornetRSquaredThres]
    if res.shape[0] == 0:
        raise ValueError(
            '\n\n no significant genes given your definition of "significant"'
        )
    res = res[res["Pratt"] > args.hornetPrattThres]
    if res.shape[0] == 0:
        raise ValueError(
            '\n\n no significant genes given your definition of "significant"'
        )
    chrs = res["Chromosome"].unique()
    novelGenes = []
    novelchrs = []
    closestSNPs = []
    closestBPs = []
    closestPs = []
    for _ in range(0, len(chrs)):  # for each chromosome in the results
        novelGenes = utils.flatten_list(novelGenes)
        novelchrs = utils.flatten_list(novelchrs)
        closestSNPs = utils.flatten_list(closestSNPs)
        closestBPs = utils.flatten_list(closestBPs)
        closestPs = utils.flatten_list(closestPs)
        #
        reschrdf = res[res["Chromosome"] == chrs[_]]
        clumpchrdf = oc[oc["CHR"] == chrs[_]]
        if clumpchrdf.shape[0] == 0:
            novelGenes.append(reschrdf["Gene"].unique().tolist())
            novelchrs.append(reschrdf["Chromosome"].values.tolist())
            closestSNPs.append([float("nan")] * reschrdf.shape[0])
            closestBPs.append([float("nan")] * reschrdf.shape[0])
            closestPs.append([float("nan")] * reschrdf.shape[0])
            continue
        for j in range(0, reschrdf.shape[0]):
            # is the significant gene in HORNET within args.novelKbWindow of any known signals?
            alldist = abs((reschrdf["geneBP"].values)[j] - clumpchrdf["BP"].values)
            # this code is actually recording all genes meeting this criteria, even if they are in the same locus as each other
            if all(alldist > (args.novelKbWindow * 1000)):
                novelGenes.append((reschrdf["Gene"].values)[j])
                novelchrs.append((reschrdf["Chromosome"].values)[0])
                closeix = np.argmin(alldist)
                closestSNPs.append((clumpchrdf["SNP"].values)[closeix])
                closestBPs.append((clumpchrdf["BP"].values)[closeix])
                closestPs.append((clumpchrdf["P"].values)[closeix])
    if len(novelGenes) == 0:
        raise ValueError(
            '\n\n no significant genes given your definition of "significant"'
        )
    novelGenes = np.array(utils.flatten_list(novelGenes))
    novelchrs = np.array(utils.flatten_list(novelchrs))
    closestSNPs = np.array(utils.flatten_list(closestSNPs))
    closestBPs = np.array(utils.flatten_list(closestBPs))
    closestPs = np.array(utils.flatten_list(closestPs))

    pdf = pd.DataFrame(
        {
            "NovelGene": novelGenes,
            "Chromosome": novelchrs,
            "ClosestPhenoSNP": closestSNPs,
            "ClosestPhenoBP": closestBPs,
            "ClosestPhenoPvalue": closestPs,
        }
    )
    # add other info like pratt index etc
    res = res.rename(
        columns={"Gene": "NovelGene", "RsquaredMRJones": "LocusR2", "Q": "CausalChisq"}
    )
    pdf = pd.merge(
        pdf,
        res[["NovelGene", "geneBP", "Pratt", "CausalChisq", "LocusR2"]],
        how="left",
        left_on="NovelGene",
        right_on="NovelGene",
    )
    pdf["NovelGene"] = pdf["NovelGene"].apply(lambda x: x.split(".")[0])
    # load gene name key and attach actual gene names
    gk = pd.read_csv(
        os.path.abspath("data/biomart_gene_legend.csv.gz"), sep=None, engine="python"
    )
    pdf = pd.merge(
        pdf,
        gk.rename(columns={"EnsemblID": "NovelGene"}),
        how="left",
        left_on="NovelGene",
        right_on="NovelGene",
    )
    pdf = pdf.rename(columns={"GeneName": "NovelGeneName"})
    pdf = pdf[
        [
            "NovelGene",
            "NovelGeneName",
            "Chromosome",
            "geneBP",
            "CausalChisq",
            "Pratt",
            "LocusR2",
            "ClosestPhenoSNP",
            "ClosestPhenoBP",
            "ClosestPhenoPvalue",
        ]
    ]
    doprint = (args.print.lower() == "yes") | (args.print.lower() == "true")
    if doprint:
        print(pdf)
    # save out
    fpout = os.path.abspath(args.out) + ".csv"
    pdf.to_csv(fpout, index=False)
    print("\n Novel genes are saved to {}".format(fpout))
    # remove stuff I wrote out if on linux/mac
    if platform.system() != "Windows":
        cmd = [utils.callDelete(), os.path.abspath(os.getcwd()) + "/a.log"]
        _ = subprocess.call(
            cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
        )  # executing the command in the terminal
        cmd = [utils.callDelete(), os.path.abspath(os.getcwd()) + "/a.clumped"]
        _ = subprocess.call(
            cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
        )  # executing the command in the terminal
        if args.phenoP == "missing":
            cmd = [utils.callDelete(), os.path.abspath(os.getcwd()) + "/a.txt"]
            _ = subprocess.call(
                cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )  # executing the command in the terminal

    if args.message.lower() in ["yes", "true"]:
        print(
            '\n\n NOTE! Use of this is only informative if you used HORNET to analyze loci that \n could have been "less significant" than what you are calling a "significant" outcome locus here\n\n you can silence this message by adding the "-m no" flag.\n\n'
        )

    return None


if __name__ == "__main__":
    args = parse_args()
    main(args)
