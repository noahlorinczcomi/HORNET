# standard imports
import pandas as pd
import numpy as np

# custom imports
from MNet import *  # this is a source file in the HORNET directory


def geneGroupFinder(geneGroups, gene, isGtex=False):
    keynames = list(geneGroups)
    # if from GTEx, gene IDs will have a "." for some reason ...
    found = 0
    for s in range(0, len(keynames)):
        fl = geneGroups[keynames[s]]
        if isGtex:
            fl = [fl[i].split(".")[0] for i in range(0, len(fl))]
            gene = gene.split(".")[0]
        if gene in fl:
            found = 1
            break
    if found == 1:
        out = keynames[s]
        return out  # returns keyname of group the gene is in
    else:
        return "cannot find"


# function to load outcome data and determine the chromosome
def loadOutcomeGWASData(fp, effectAllele, z, rsid, ldRefDir):
    if "@" in z:  # if no Z-stat column in phenotype data
        sep = z.split("@", 1)  # user must have passed 'BETA@SE' for Z-stat column name
        uc2 = [
            rsid,
            sep[0],
            sep[1],
            effectAllele,
        ]  # only pick out these columns from pheno GWAS data
        dataPheno = pd.read_table(
            fp, sep=None, engine="python", usecols=uc2
        )  # load pheno GWAS data
        dataPheno["phenoZ"] = (
            dataPheno[sep[0]] / dataPheno[sep[1]]
        )  # create my own Z-stat column
        dataPheno = dataPheno.rename(
            columns={uc2[0]: "phenoSNP", uc2[3]: "phenoEffectAllele"}
        )  # rename to names I'll now later
        dataPheno = dataPheno[
            ["phenoSNP", "phenoZ", "phenoEffectAllele"]
        ]  # drop unneccessary columns
    else:
        uc2 = [rsid, z, effectAllele]
        dataPheno = pd.read_table(fp, sep=None, engine="python", usecols=uc2)
        dataPheno = dataPheno.rename(
            columns={uc2[0]: "phenoSNP", uc2[1]: "phenoZ", uc2[2]: "phenoEffectAllele"}
        )
    # filter to only LD ref dir SNPs
    onekgSNPs = pd.read_table(
        ldRefDir + ".bim", sep=None, engine="python", header=None
    )  # has column names CHR, SNP (SNP is for rsID)
    onekgSNPs = onekgSNPs[
        [1]
    ]  # the SNP column name is actually an index for .bim files
    merged = pd.merge(dataPheno, onekgSNPs, left_on="phenoSNP", right_on=1)
    merged = merged[["phenoSNP", "phenoEffectAllele", "phenoZ"]]
    return (
        merged  # this is dataPeno and 1kg merged; NOT exposure and pheno/outcome merged
    )


# function to load exposure data and merge with outcome data (separated from outcome data bc many exposure files will be loaded)
# with specifying sample size column name
def loadExposureGWASData(
    fp,
    effectAllele="",
    z="",
    rsid="",
    snpBP="",
    geneLabel="",
    geneBPLabel="",
    nLabel="",
    isRawGzippedGTEx=False,
    mapwd="data/maps/1kgPhase3maps",
    nr=1e100,
):  # outcome data already filtered to 1kg SNPs
    dataLoaded = False
    # If this is raw GTEx data, I already know what to expect for column names
    if isRawGzippedGTEx == True:  # I take care of all column names here
        from io import BytesIO  # Converting bytes to bytes input file
        import gzip
        import pyarrow  # Fast reading of parquets

        file = gzip.open(fp, "rb")
        # speed of loading full parquet data seems to be same as loading some subset
        # pf=pyarrow.parquet.ParquetFile(file);first_ten_rows=next(pf.iter_batches(batch_size=10));df=pa.Table.from_batches([first_ten_rows]).to_pd()
        df = file.read()
        df = pd.read_parquet(BytesIO(df))
        dataLoaded = True
        dataGene = attachRsidToGTEx(df, mapwd=mapwd)
        dataGene["geneZ"] = dataGene["slope"] / dataGene["slope_se"]
        dataGene = dataGene.rename(
            columns={
                "rsID": "geneSNP",
                "phenotype_id": "Gene",
                "ALT": "geneEffectAllele",
            }
        )  # already has snpBP column from attachRsidToGTEx()
        if is_number(nLabel):
            dataGene["geneN"] = nLabel
        else:
            dataGene = dataGene.rename(columns={nLabel: "geneN"})
        dataGene["tempGeneBP"] = (
            dataGene["snpBP"] - dataGene["tss_distance"]
        )  # snpBP is in hg38 since it is from GTEx
        aggd = dataGene.groupby("Gene").agg({"tempGeneBP": "mean"})
        aggd = pd.DataFrame(
            {"geneBP": aggd.tempGeneBP.values, "Gene": list(aggd.index)}
        )
        dataGene = pd.merge(dataGene, aggd, left_on="Gene", right_on="Gene")
        dataGene = dataGene[
            ["geneSNP", "geneZ", "Gene", "geneEffectAllele", "geneBP", "snpBP", "geneN"]
        ]  # drop unneccessary columns
    else:  # I need to take care of the column names now
        if "@" in z:
            sep = z.split(
                "@", 1
            )  # user must have passed 'BETA@SE' for Z-stat column name
            uc1 = [
                rsid,
                sep[0],
                sep[1],
                geneLabel,
                effectAllele,
                geneBPLabel,
                snpBP,
            ]  # only pick out these columns from exposure GWAS data
            if dataLoaded == False:
                dataGene = pd.read_table(fp, sep=None, engine="python", nrows=nr)
                # load exposure GWAS data
            dataGene["geneZ"] = (
                dataGene[sep[0]] / dataGene[sep[1]]
            )  # create my own Z-stat column
            dataGene = dataGene.rename(
                columns={
                    uc1[0]: "geneSNP",
                    uc1[3]: "Gene",
                    uc1[4]: "geneEffectAllele",
                    uc1[5]: "geneBP",
                    uc1[6]: "snpBP",
                }
            )
            if is_number(nLabel):
                dataGene["geneN"] = nLabel
            else:
                dataGene = dataGene.rename(columns={nLabel: "geneN"})
            dataGene = dataGene[
                [
                    "geneSNP",
                    "geneZ",
                    "Gene",
                    "geneEffectAllele",
                    "geneBP",
                    "snpBP",
                    "geneN",
                ]
            ]  # drop unneccessary columns
        else:
            uc1 = [rsid, z, geneLabel, effectAllele, geneBPLabel, snpBP]
            if dataLoaded == False:
                dataGene = pd.read_table(fp, sep=None, engine="python")
            dataGene = dataGene.rename(
                columns={
                    uc1[0]: "geneSNP",
                    uc1[1]: "geneZ",
                    uc1[2]: "Gene",
                    uc1[3]: "geneEffectAllele",
                    uc1[4]: "geneBP",
                    uc1[5]: "snpBP",
                }
            )
            if is_number(nLabel):
                dataGene["geneN"] = nLabel
            else:
                dataGene = dataGene.rename(columns={nLabel: "geneN"})
            dataGene = dataGene[
                [
                    "geneSNP",
                    "geneZ",
                    "Gene",
                    "geneEffectAllele",
                    "geneBP",
                    "snpBP",
                    "geneN",
                ]
            ]  # drop unneccessary columns
    # if data is from MetaBrain, change SNP column to only contain rsIDs
    dataGene["geneSNP"] = dataGene["geneSNP"].apply(lambda x: rsidFinder(x))
    return dataGene


# merge exposure and outcome data
def mergeExposureAndOutcomeGWAS(dataGene, dataPheno):
    merged = pd.merge(
        left=dataGene, right=dataPheno, left_on="geneSNP", right_on="phenoSNP"
    )
    # harmonise exposure and outcome effect sizes by effect alleles
    mask = (
        merged.geneEffectAllele.str.upper() != merged.phenoEffectAllele.str.upper()
    )  # find where [uppercase] alleles do not match
    merged.loc[mask, "geneZ"] = (
        -1 * merged.loc[mask, "geneZ"]
    )  # change sign of Z-stat for gene expression if effect alleles do not match
    # implicitly makes the outcome effect allele the effect allele
    return merged


# function to attach rsIDs to GTEx data that is in hg38 and SNPs are named as eg chr9_10536_C_T_b38
def attachRsidToGTEx(dataGTEx, mapwd="data/maps/1kgPhase3maps"):
    # the mapfile is in a pre-determined location that will be set when the user clones our github repo:
    # alt allele in GTEx is effect allele
    # separate information in rsidcolumnname of dataGTEx
    splitted = dataGTEx["variant_id"].str.split("_", n=10)
    splitted = np.array(splitted)
    marker_ = []
    chr_ = []
    bp_ = []
    alt_ = []
    ref_ = []
    for _ in range(0, splitted.shape[0]):
        chr_.append(int(splitted[_][0].replace("chr", "")))
        bp_.append(int(splitted[_][1]))
        ref_.append(splitted[_][2])
        alt_.append(splitted[_][3])
        marker_.append(str(chr_[_]) + "_" + str(bp_[_]))
    dataGTEx["markername"] = marker_
    dataGTEx["CHR"] = chr_
    dataGTEx["REF"] = ref_  # GTEx and merged_ will now have these columns, too
    dataGTEx["ALT"] = alt_  # effect allele according to documentation
    dataGTEx["snpBP"] = bp_
    # load mapfile corresponding to this chromosome (should all be one CHR at a time from dataGene so using any values in chr_ will work)
    # names of mapfiles do not change and are not set by the user
    mapfile = pd.read_csv(
        mapwd + "/" + "1kgPhase3hg19Chrom" + str(chr_[0]) + "MapToHg38.txt.gz",
        engine="python",
        sep=None,
    )  # uses chromosome-specific mapfiles. chromsome determined using data
    mapfile["markername"] = (
        mapfile["CHROM_b38"].astype(str) + "_" + mapfile["POS_b38"].astype(str)
    )  # make out of hg38 values bc that's what's in GTEx
    mapfile = mapfile[["rsID", "markername"]]  # drop unnecessary columns
    # merge
    merged_ = pd.merge(dataGTEx, mapfile, left_on="markername", right_on="markername")
    return merged_


### !!! NOTE !!! this function is pretty slow ... loading phenotype data always seems to take a long time
def attachRsidToPhenoData(
    fpdata,
    chrom,
    bp,
    effectAllele,
    z,
    rsid,
    ldRefDir,
    build="NA",
    mapwd="data/maps/fullMaps/1kgPhase3MapAllSNPs.txt",
    nSNPsToTest=1000,
):
    if (is_number(build) == False) & (build.lower() != "na"):
        raise ValueError(
            "please enter either 37 (hg19) or 38 (hg38) as the build parameter"
        )
    data = pd.read_csv(fpdata, sep=None, engine="python")
    mapdf = pd.read_csv(mapwd, sep=None, engine="python")
    # all mapfiles should have rsID  CHROM_b37  CHROM_b38   POS_b37   POS_b38 (I know this because I will be the one to create them)
    if build.lower() == "na":  # if the user does not know the build
        build = determineBuild(data, chrom, bp, mapdf, nSNPsToTest=nSNPsToTest)
    mapdf = mapdf.rename(
        columns={"CHROM_b" + str(build): "usechr", "POS_b" + str(build): "usebp"}
    )  # identifying columns in map corresponding to build
    mapdf["usemarker__"] = (
        mapdf["usechr"].astype(str) + ":" + mapdf["usebp"].astype(str)
    )  # making markername column to merge map on
    data["usemarker__"] = (
        data[chrom].astype(str) + ":" + data[bp].astype(str)
    )  # making markername column to merge data on
    data = data.drop_duplicates(
        subset=["usemarker__"]
    )  # cannot have duplicates when merging
    dataPheno = pd.merge(
        data,
        mapdf[["rsID", "usemarker__"]],
        left_on="usemarker__",
        right_on="usemarker__",
    )  # merge data and map
    nMapped = dataPheno.shape[0]  # number of mapped SNPs
    dataPheno = dataPheno.drop(
        columns=["usemarker__"]
    )  # dropping markername column I made
    del data  # to save memory
    # note that any SNPs in data that are not in the map are dropped - you should therefore choose a map consistent with the LD ref panel since all
    # SNPs will be subsetted to those in the LD ref panel anyway
    #
    # instead of saving, I should now give it to the user in a cleaned way as I would if they used loadPhenoData()
    if "@" in z:  # if no Z-stat column in phenotype data
        sep = z.split(
            "@", 1
        )  # user must have passed (eg) 'BETA@SE' for Z-stat column name
        uc2 = [
            rsid,
            sep[0],
            sep[1],
            effectAllele,
        ]  # only pick out these columns from pheno GWAS data (rsid should be rsid in the map. if i made the map then it is 'rsID')
        dataPheno["phenoZ"] = (
            dataPheno[sep[0]] / dataPheno[sep[1]]
        )  # create my own Z-stat column
        dataPheno = dataPheno.rename(
            columns={uc2[0]: "phenoSNP", uc2[3]: "phenoEffectAllele"}
        )  # rename to names I'll now later
        dataPheno = dataPheno.rename(columns={"rsID": "phenoSNP"})[
            ["phenoSNP", "phenoZ", "phenoEffectAllele"]
        ]  # drop unneccessary columns. mapdf always attacheds 'rsID'
        dataPheno["phenoEffectAllele"] = dataPheno[
            "phenoEffectAllele"
        ].str.upper()  # make effect allele uppercase
    else:
        uc2 = ["rsID", z, effectAllele]
        dataPheno = dataPheno.rename(
            columns={uc2[0]: "phenoSNP", uc2[1]: "phenoZ", uc1[2]: "phenoEffectAllele"}
        )
        dataPheno["phenoEffectAllele"] = dataPheno[
            "phenoEffectAllele"
        ].str.upper()  # make effect allele uppercase
    # filter to only LD ref dir SNPs
    onekgSNPs = pd.read_table(
        ldRefDir + ".bim", sep=None, engine="python", header=None
    )  # has column names CHR, SNP (SNP is for rsID)
    onekgSNPs = onekgSNPs[
        [1]
    ]  # the SNP column name is actually an index for .bim files
    merged = pd.merge(dataPheno, onekgSNPs, left_on="phenoSNP", right_on=1)
    merged = merged.drop(columns=1)
    return merged


def determineBuild(data, chrom, bp, mapdf, nSNPsToTest=1000):
    # mapdf always has columns rsID  CHROM_b37  CHROM_b38   POS_b37   POS_b38
    mapdf["map__marker37"] = (
        mapdf["CHROM_b37"].astype(str) + ":" + mapdf["POS_b37"].astype(str)
    )  # funny name so nobody else has it picked already
    mapdf["map__marker38"] = (
        mapdf["CHROM_b38"].astype(str) + ":" + mapdf["POS_b38"].astype(str)
    )
    mask = data[chrom].values == mapdf["CHROM_b37"].values[0]
    datacut = data[
        mask
    ]  # data will be for all CHRs since it is outcome data, but map data will be CHR-specific. This is fine, just subset outcome data this this CHR
    datacut = datacut.iloc[
        : (nSNPsToTest * 100), :
    ]  # assuming only 1/100th of outcome GWAS SNPs will be in map
    datacut["data__marker"] = datacut[chrom].astype(str) + ":" + datacut[bp].astype(str)
    m1 = pd.merge(
        datacut, mapdf, left_on="data__marker", right_on="map__marker37"
    )  # for build 37
    m2 = pd.merge(
        datacut, mapdf, left_on="data__marker", right_on="map__marker38"
    )  # for build 38
    if (m1.shape[0] / 2) > m2.shape[0]:
        out = 37
    elif m1.shape[0] < (m2.shape[0] / 2):
        out = 38
    else:
        out = "cannot reliably determine build - please increase nSNPsToTest until this message disappears"
    return out


def rsidFinder(st):
    # x should be a string
    sp = st.split(":")
    if len(sp) == 1:
        return st
    w = ["rs" in sp[_].lower() for _ in range(0, len(sp))]
    if True in w:
        out = sp[w.index(True)]
    else:
        out = float("nan")
    return out
