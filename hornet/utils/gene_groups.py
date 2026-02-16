# standard imports
import numpy as np

# custom imports
from MNet import *  # this is a source file in the HORNET directory


# identify gene groups for multivariable MR
def defGeneGroups(q0geneGroups, merged, allowSingletons=False):
    # merged=merged.sort_values(by='geneBP',axis=0) # sort by BP position of genes - not sure why this matters but it does
    genes = merged["Gene"].unique().tolist()
    # genes=list(genes)
    mergedSig = merged[abs(merged["geneZ"]) > q0geneGroups]
    usedGenes = []
    lens = []
    geneGroups = {}
    for ii in range(0, len(genes)):
        if (genes[ii] in usedGenes) == True:
            continue
        else:
            qts = mergedSig[mergedSig["Gene"] == genes[ii]][
                "geneSNP"
            ].values  # all SNPs that are cis-eQTLs for this gene (genes[ii])
        if len(qts) == 0:  # if no cis-eQTLs, skip gene
            continue
        qtsGenes = mergedSig[mergedSig["geneSNP"].isin(qts)][
            "Gene"
        ].unique()  # all genes for which these SNPs (qts) are cis-eQTLs
        if (len(qtsGenes) == 1) & (
            allowSingletons == False
        ):  # I will not form a group of 1 ... why not?
            continue
        toAdd = list(qtsGenes)  # genes sharing one of these eQTLs
        boo = [
            (toAdd[x] in usedGenes) == False for x in range(0, len(toAdd))
        ]  # boolean: True=gene has not been grouped; False=gene has been grouped already
        if (
            sum(boo) == 0
        ):  # If all genes sharing this eQTL are already in another group, next
            continue
        else:
            toAdd = (
                np.array(toAdd)[boo]
            ).tolist()  # else, only use the ungrouped genes to form a group with the current gene
        if (allowSingletons == False) & (sum(boo) == 1):
            continue
        usedGenes.append(toAdd)
        usedGenes = flatten_list(usedGenes)
        geneGroups[genes[ii]] = (
            toAdd  # geneGroups is a dictionary with names that are the index genes in a group
        )
        lens.append(len(toAdd))
    ggKeys = [list(geneGroups)[x] for x in range(0, len(geneGroups))]
    return geneGroups, ggKeys, lens, usedGenes


def defGeneGroupsByOutcome(q0geneGroups, merged, KbWindow=100, closestK=2):
    # assumes the "isOutcomeClump" column is in merged
    # KbWindow should be set to the same, or smaller, value as the value for outcomeClumpingKBWindow
    clumpBPs = merged[merged["isOutcomeClump"] == True]["snpBP"].unique()
    usedGenes = []
    lens = []
    geneGroups = {}
    mergedcut = merged.copy()[
        abs(merged["geneZ"]) > q0geneGroups
    ]  # subsetting only to significant eQTLs
    for i in range(0, len(clumpBPs)):
        cutdf = mergedcut.copy()[
            abs(mergedcut["snpBP"] - clumpBPs[i]) < (KbWindow * 1e3)
        ]  # subsetting to only those within KbWindow either side of outcome clump
        cutdf["absDist"] = abs(
            cutdf["snpBP"].values - clumpBPs[i]
        )  # calculating absolute distance from eQTL to outcome clump
        closestKgenes = (cutdf["Gene"].unique())[
            :closestK
        ]  # it puts them in order of appearance in data, which is now sorted ascending by distance from outcome clump
        # now find genes sharing eQTLs with any of these closestK genes closest to the clump
        theseeqtls = cutdf[cutdf["Gene"].isin(closestKgenes)][
            "geneSNP"
        ]  # find eQTLs for closest closestKgenes genes
        genestogroup = cutdf[cutdf["geneSNP"].isin(theseeqtls)][
            "Gene"
        ].unique()  # this is the group that is formed
        boo = [
            (genestogroup[x] in usedGenes) == False for x in range(0, len(genestogroup))
        ]  # True if not grouped, False otherwise
        if sum(boo) < 1:  # this will allow singletons, but will skip empty ones
            continue
        genestogroup = genestogroup[boo]  # only group those not previously grouped
        geneGroups["group" + str(i)] = (
            genestogroup  # key is a non-informative number - else there can be issues for some reason
        )
        usedGenes.append(
            genestogroup.tolist()
        )  # add newly grouped genes to running list of those already used
        usedGenes = flatten_list(usedGenes)
    ggKeys = list(geneGroups)
    lens = [len(geneGroups[ggKeys[x]]) for x in range(0, len(ggKeys))]
    return geneGroups, ggKeys, lens, usedGenes
