# standard imports
import warnings
import os
import tqdm
from time import time

# custom imports
import constants
from cli import parse_args

from utils import (

)
from MNet import (

)


def main(args):
    """ """
    dontHide = (args.hideVersion.lower() == "yes") | (
        args.hideVersion.lower() == "true"
    )
    if dontHide:
        f = open("hornet.txt", "r")
        con = f.read()
        print(con)
        f.close()

    print(" ")
    print("HORNET started " + time.ctime())

    ### first checking if all of the files/direcetories they gave actually exist
    fileChecker(os.path.abspath(args.eQTLGWAS), "eQTL GWAS directory")
    fileChecker(os.path.abspath(args.phenoGWAS), "phenotype GWAS filepath")
    fileChecker(
        os.path.abspath(args.LDRef + ".bed"),
        "LD reference panel filepath (w/o extension)",
    )
    fileChecker(os.path.abspath(args.writableDir), "temporing working directory")
    if args.whereSaveRawData != "":
        fileChecker(os.path.abspath(args.whereSaveRawData), "saved raw data directory")
    od = (args.out.split("/"))[:-1]
    fileChecker("/".join(od), "output/results directory")
    # if no files in eQTL GWAS directory, there's a problem
    if (os.listdir(args.eQTLGWAS)) == 0:
        raise ValueError(
            "It looks like there are no files in the eQTL GWAS directory ({})".format(
                args.eQTLGWAS
            )
        )
    # output directory for graphs
    # graphoutdir=args.networkGraphsOut+'/' if args.networkGraphsOut!='plots' else os.getcwd()+'/plots/'
    if os.path.isdir(os.path.abspath(args.networkGraphsOut)) == False:
        raise ValueError(
            "\n It looks like the directory you want to save network graph plots to does not exist({})".format(
                args.networkGraphsOut
            )
        )
    # check other flags
    if (args.geneScreener in ["gscreen", "mrjones"]) == False:
        raise ValueError(
            'the "--geneScreener" flag accepts either "gscreen" or "mrjones", but you gave it {}'.format(
                args.geneScreener.lower()
            )
        )

    ##
    candidateGenesIsFile = os.path.isfile(args.candidateGenes)
    if candidateGenesIsFile == False:
        candidateGenes = (
            args.candidateGenes.split(",") if len(args.candidateGenes) > 0 else []
        )
    else:
        candidateGenes = (
            pandas.read_table(args.candidateGenes, header=None)
            .iloc[:, 0]
            .values.tolist()
        )

    # drop Ensembl version type if it is present (if it is not, this code will have no effect)
    if len(candidateGenes) > 0:
        candidateGenes = [
            candidateGenes[i].split(".")[0] for i in range(0, len(candidateGenes))
        ]

    snplabs = args.snpLabels.split(",")
    zlabs = args.zLabels.split(",")
    ealabs = args.effectAlleles.split(",")
    startwd = (
        os.getcwd()
    )  # starting working directory - will be the absolute path of the HORNET software
    impute = (args.imputeMissingEQTLs.lower() == "true") | (
        args.imputeMissingEQTLs.lower() == "yes"
    )

    ### exposure parameters
    dirGene = os.path.abspath(args.eQTLGWAS)
    isRawGTEx = (args.isRawGTEx.lower() == "true") | (args.isRawGTEx.lower() == "yes")
    effectAlleleGene = ealabs[0]
    zGene = zlabs[0]
    rsidGene = snplabs[0]
    geneLabelGene = args.eQTLGeneLabel
    geneBPLabelGene = args.eQTLGeneBP
    snpBPGene = args.eQTLSNPBPLabel
    nLabel = args.eQTLGWASn
    # if GTEx, you need to choose these carefully or else no analysis can be done
    ldUpperLimit = args.ldMax
    ldOtherLociOtherPt = args.POtherLoci
    ldOtherLociR2 = (args.LDOtherLoci) ** 2
    ldOtherLociWindow = args.otherLociMbWindow
    q0Correls = 1 - norm.ppf(args.PvalueThresholdBiasCorrection)
    nMinCorrels = args.nSNPsBiasCorrection
    jointChiGenesP = args.MVMRIVPThreshold
    assumedMissingMean = 0
    opAlpha = 0
    verbose = True
    nMinIVs = args.minMVMRIVs  # default is 30
    hessMinScale = args.hessMinScale  # default is 2.5
    silence = (args.silence.lower() == "true") | (args.silence.lower() == "yes")
    UniMRIVPThreshold = args.UVMRIVPThreshold
    adjustSEsForInflation = (args.adjustSEsForInflation.lower() == "true") | (
        args.adjustSEsForInflation.lower() == "yes"
    )
    analysisInPhenotypeLoci = (args.analysisInPhenotypeLoci.lower() == "true") | (
        args.analysisInPhenotypeLoci.lower() == "yes"
    )
    outcomeLociMbWindow = args.onlyWithinMbOfPhenoLocus
    q0geneGroups = norm.ppf(1 - args.geneGroupPvalueThreshold)
    assumeNoSampleOverlap = (args.assumeNoSampleOverlap.lower() == "true") | (
        args.assumeNoSampleOverlap.lower() == "yes"
    )
    shrinkBiasCorrection = True  # I don't think we want this adjusted; keep it simple
    saveData = (args.saveRawData.lower() == "true") | (
        args.saveRawData.lower() == "yes"
    )
    numIndexGenesToFormGroup = args.numIndexGenesToFormGroup

    ### phenotype parameters
    fpPheno = os.path.abspath(args.phenoGWAS)
    writableDir = os.path.abspath(args.writableDir)
    ldRefDir = os.path.abspath(args.LDRef)
    # mapwd is fixed - users must use one of the maps I created
    mapwd = os.path.abspath("data/maps/1kgPhase3maps")
    rsidPheno = snplabs[1]
    effectAllelePheno = ealabs[1]
    zPheno = zlabs[1]
    nLabel = args.phenotypeGWASn
    outcomeClumpingKBWindow = args.phenoLociMbWindow
    outcomeClumpingPthreshold = args.phenoLociPvalue
    outcomeClumpingR2 = args.phenoLociR2
    ### loading phenotype data
    print("Loading phenotype GWAS data")
    dataPheno = loadOutcomeGWASData(
        fpPheno, effectAllelePheno, zPheno, rsidPheno, ldRefDir
    )
    dataPheno = findOutcomeSignals(
        dataPheno,
        ldRefDir,
        writableDir,
        outcomeClumpingKBWindow,
        outcomeClumpingPthreshold,
        outcomeClumpingR2,
    )
    dataPheno[dataPheno["isOutcomeClump"] == True].to_csv(
        os.path.abspath("results/outcomeloci.csv")
    )

    ### loading key/dictionary/lookup/reference data
    bim = pandas.read_csv(
        ldRefDir + ".bim", sep="\t", names=["chr", "rsid", "x", "bp", "a1", "a2"]
    )  # for figuring out which chromosome is being used later
    genekey = pandas.read_csv(
        os.getcwd() + "/data/biomart_gene_legend.csv.gz", engine="python"
    )  # to be used only if candidateGenes are specified

    ### running
    # os.chdir(dirGene)
    runningres = pandas.DataFrame()  # to be filled in
    runningdiagnostics = pandas.DataFrame()  # to be filled in
    infls = list()
    edgeLists = {}
    usedCandidateGenes = []
    for _ in range(0, len(os.listdir(dirGene))):
        t0 = time.perf_counter()  # start timer
        # make sure I haven't already used all candidate genes
        mask = [
            candidateGenes[i] in usedCandidateGenes
            for i in range(0, len(candidateGenes))
        ]
        if all(mask) & (
            len(candidateGenes) > 0
        ):  # if all candidate genes have already been used
            print("All available candidate genes have been tested")
            break  # don't need to consider additional chromosomes
        fpGene = os.listdir(dirGene)[_]  # define filepath for exposure data
        dataGene = loadExposureGWASData(
            dirGene + "/" + fpGene,
            effectAlleleGene,
            zGene,
            rsidGene,
            snpBPGene,
            geneLabelGene,
            geneBPLabelGene,
            nLabel=nLabel,
            isRawGzippedGTEx=isRawGTEx,
            mapwd=mapwd,
        )  # load genes
        merged = mergeExposureAndOutcomeGWAS(
            dataGene, dataPheno
        )  # merge with outcome/phenotype data
        merged = pandas.merge(
            merged, bim[["rsid", "chr", "a1"]], left_on="geneSNP", right_on="rsid"
        )
        # harmonise effect sizes to the 'effect' allele in the bim file, which is 'a1' (effect sizes currently harmonized to phenoEffectAllele)
        mask = merged.phenoEffectAllele.str.upper() != merged.a1.str.upper()
        merged.loc[mask, "geneZ"] = (
            -1 * merged.loc[mask, "geneZ"]
        )  # harmonizing geneZ to bim effect allele
        merged.loc[mask, "phenoZ"] = (
            -1 * merged.loc[mask, "phenoZ"]
        )  # harmonizing phenoZ to bim effect allele
        chromosome = merged["chr"].unique()
        chromosome = int(chromosome)
        del dataGene
        # delete data we no longer need
        # check merged data for candidate gene(s) if user provided any (not much speed advantage checking merged vs dataGene bc merging is very fast)
        merged["Gene"] = merged["Gene"].apply(lambda x: x.split(".")[0])
        if len(candidateGenes) > 0:
            ch = merged.copy()
            ch = ch[
                ch["Gene"].isin(candidateGenes)
            ]  # candidateGenes already got splitted by '.'
            usedCandidateGenes.append(
                ch["Gene"].unique().tolist()
            )  # keeping track of which candidate genes I've used so I know when to stop
            usedCandidateGenes = flatten_list(usedCandidateGenes)
            if ch.shape[0] == 0:
                print(
                    "skipping chromsome "
                    + str(chromosome)
                    + " because there are no candidate genes on this chromosome"
                )
                del ch
                continue

        print("Starting chromosome " + str(chromosome))
        if (analysisInPhenotypeLoci == True) & (
            len(candidateGenes) == 0
        ):  # form gene groups
            geneGroups, ggKeys, lens, usedGenes = defGeneGroupsByOutcome(
                q0geneGroups,
                merged,
                KbWindow=outcomeClumpingKBWindow,
                closestK=numIndexGenesToFormGroup,
            )
        else:
            geneGroups, ggKeys, lens, usedGenes = defGeneGroups(q0geneGroups, merged)
        # [print(geneGroupFinder(geneGroups,candidateGenes[i],isGtex=True if isRawGTEx else False)) for i in range(0,len(candidateGenes))]
        thingsMonitored, outerDict, edgeD = MVMRworkhorse(
            merged,
            geneGroups,
            ggKeys,
            writableDir=writableDir,
            ldRefDir=ldRefDir,
            isGtex=isRawGTEx,
            analysisOnlyInOutcomeLoci=analysisInPhenotypeLoci,
            outcomeLociMbWindow=outcomeLociMbWindow,
            ldUpperLimit=ldUpperLimit,
            ldOtherLociOtherPt=ldOtherLociOtherPt,
            ldOtherLociR2=ldOtherLociR2,
            ldOtherLociWindow=ldOtherLociWindow,
            q0Correls=q0Correls,
            nMinCorrels=nMinCorrels,
            jointChiGenesP=jointChiGenesP,
            opAlpha=opAlpha,
            nMinIVs=nMinIVs,
            hessMinScale=hessMinScale,
            silence=silence,
            UniMRIVPThreshold=UniMRIVPThreshold,
            candidateGenes=candidateGenes,
            assumeNoSampleOverlap=assumeNoSampleOverlap,
            shrinkBiasCorrection=shrinkBiasCorrection,
            networkR2Thres=args.networkR2Thres,
            impute=impute,
            saveData=saveData,
            geneSelector=args.geneScreener.lower(),
        )
        edgeLists["Chromosome" + str(chromosome)] = edgeD  # save edges
        if len(outerDict) == 0:
            print(
                "  Skipping chromosome "
                + str(chromosome)
                + " because of insufficient phenotype and/or eQTL signals"
            )
            continue
        res = organizeMetaResults(outerDict)
        res["CHR"] = chromosome
        res["infl"] = 1
        res = prepRes(res)  # store results
        diagnostics = organizeThingsMonitored(thingsMonitored)  # store diagnostics
        diagnostics["Chromosome"] = chromosome
        runningdiagnostics = pandas.concat([runningdiagnostics, diagnostics.copy()])
        if adjustSEsForInflation:  # should SEs be adjusted for inflation?
            infl, ccs, ccps, ug, ms, bps = calcInflationFactor(
                merged, ldRefDir, writableDir, mbWindow=1, ldUpperLimit=ldUpperLimit
            )
            infls.append(infl)
            res["nullInflation"] = numpy.max((1, infl))
            res = adjustInflation(res, numpy.max((1, infl)))
        runningres = pandas.concat([runningres, res.copy()])
        # if user wanted to save raw data, do that now inside of the chosen directory
        if saveData:
            fpout = os.path.abspath(args.whereSaveRawData) + "/"
            print(
                'raw data are saved in the "bybx.txt", "UU.txt", "LD.txt", "snps.txt", "genes.txt", and "geneBPs.txt" files in the {} directory'.format(
                    fpout
                )
            )
            ln = list(outerDict)[0]
            numpy.savetxt(
                fpout + "bybx.txt",
                numpy.column_stack((outerDict[ln]["by"], outerDict[ln]["bX"])),
            )
            numpy.savetxt(fpout + "UU.txt", outerDict[ln]["UU"])
            numpy.savetxt(fpout + "LD.txt", outerDict[ln]["regLD"])
            numpy.savetxt(
                fpout + "effectAlleles.txt",
                numpy.array(outerDict[ln]["data"]["EffectAlleles"]),
                fmt="%s",
            )
            numpy.savetxt(
                fpout + "genes.txt",
                numpy.array(outerDict[ln]["data"]["Gene"]),
                fmt="%s",
            )
            numpy.savetxt(
                fpout + "snps.txt", numpy.array(outerDict[ln]["data"]["IVs"]), fmt="%s"
            )
            numpy.savetxt(
                fpout + "geneBPs.txt",
                numpy.array(outerDict[ln]["data"]["GeneBP"]),
                fmt="%s",
            )

        # LDLeadGeneSNPs(res,thischr,writableDir,ldRefDir) # write out results to be read in by an R program that will make figures etc.
        t1 = time.perf_counter() - t0
        print(
            "  Chromosome "
            + str(chromosome)
            + " complete "
            + str(round(t1 / 60, 1))
            + " minutes later"
        )
        # iteratively save so I don't lose any progress between genes
        if (args.iterativelySave.lower() == "true") | (
            args.iterativelySave.lower() == "yes"
        ):
            runningres.to_csv(args.out + "_tempresults.txt", sep="\t")
            runningdiagnostics.to_csv(args.out + "_diagnostics.txt", sep="\t")

    # if no chromosomes could be analyzed, there is nothing to save - tell the user and exit
    if runningres.shape[0] == 0:
        raise ValueError(
            "\n\n No genes could be analyzed in these data.\n Consider adjusting your criteria for forming \ngene groups, IV sets, etc., or the data itself\nas the cause of this result\n\n"
        )

    fp1 = os.path.abspath(args.out) + "_results.txt"
    fp2 = os.path.abspath(args.out) + "_diagnostics.txt"
    # did the user want us to clean up the results?
    cleanres = (args.cleanResults.lower() == "yes") | (
        args.cleanResults.lower() == "true"
    )
    if cleanres:
        cr = runningres.copy()
        cr["Pratt"] = cr["MRBEEPostSelec_MVMR_Est"] * cr["MRBEE_UVMR_Est"]
        cr = cr[
            [
                "Gene",
                "geneBP",
                "Chromosome",
                "GeneSelected",
                "MRBEEPostSelec_MVMR_Est",
                "MRBEEPostSelec_MVMR_SE",
                "LocusR2",
                "Pratt",
                "CHRspecificGroupID",
            ]
        ]
        cr = cr.rename(
            columns={"MRBEEPostSelec_MVMR_Est": "Est", "MRBEEPostSelec_MVMR_SE": "SE"}
        )
        cr.to_csv(fp1, sep="\t")
    else:
        runningres.to_csv(fp1, sep="\t")

    runningdiagnostics.to_csv(fp2, sep="\t")
    # also save a copy of runningres to HORNET/res.csv so it can be read by plotres.r
    copyfpout = os.path.abspath(os.getcwd())
    runningres.to_csv(
        copyfpout + "/res.csv"
    )  # full (not super cleaned) data will still be written
    print("Results are written to " + fp1)
    print("Analysis diagnostics are written to " + fp2)
    # save executing commands of plotres.r for later - don't want to cause an early error bc users don't have R installed

    ### network construction (done in Python)
    # first find top args.networksInTopKLoci loci
    graphNetworks = (args.graphNetworks.lower() == "yes") | (
        args.graphNetworks.lower() == "true"
    )
    hasSaved = False
    if (runningres.shape[1] > 1) & (graphNetworks):
        vals = numpy.sort(runningres["LocusR2"].unique())[::-1]
        for _ in range(0, len(vals)):
            if vals[_] < args.networkR2Thres:
                continue
            dcut = runningres[runningres["LocusR2"] == vals[_]]
            genes = dcut["Gene"].values.tolist()
            genes = [genes[_].split(",")[0] for _ in range(0, len(genes))]
            if dcut.shape[0] < 1:
                continue
            chrom, grp = int(dcut["Chromosome"][0]), dcut["CHRspecificGroupID"][0]
            edges = edgeLists["Chromosome" + str(chrom)][grp]
            if sum(sum(edges)) == 0:
                continue
            G = addEdges(edges)
            fl = flatten_list([args.networkDiseaseLabel, genes])
            G = nx.relabel_nodes(G, mapping=mapNodeNames(edges, fl))  # add new labels
            # color_map=colorCoreGenes(edges) # add colors for core and peripheral genes
            fig = matplotlib.pyplot.figure(figsize=(12, 8))
            warnings.filterwarnings("ignore")
            nx.draw(
                G,
                with_labels=True,
                ax=fig.add_subplot(),
                alpha=0.5,
                font_size=10,
                node_color="skyblue",
            )  # all same color - less likely to produce an error
            leadgene = dcut.loc[abs(dcut["MRJonesEst"]).idxmax()]["Gene"]
            leadgene = leadgene.split(".")[0]
            fig.savefig(
                os.path.abspath(args.networkGraphsOut)
                + "/"
                + args.graphNetworksFilePrefix
                + "_"
                + leadgene
                + "_graph.png",
                dpi=500,
            )
            hasSaved = True

    print(" ")
    if graphNetworks & hasSaved:
        print(
            "Network graphs are saved in the "
            + os.path.abspath(graphoutdir)
            + " folder"
        )

    print(" ")
    print("HORNET ended " + time.ctime())
    return None


if __name__ == "__main__":
    args = parse_args()
    main(args)
