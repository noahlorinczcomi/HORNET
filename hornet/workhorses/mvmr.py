# perform analysis
def MVMRworkhorse(
    merged,
    geneGroups,
    ggKeys,
    writableDir,
    ldRefDir,
    isGtex=False,
    analysisOnlyInOutcomeLoci=True,
    outcomeLociMbWindow=1,
    ldUpperLimit=0.5,
    ldOtherLociOtherPt=0.0001,
    ldOtherLociR2=0.1,
    ldOtherLociWindow=1,
    q0Correls=3.290527,
    nMinCorrels=50,
    jointChiGenesP=5e-8,
    opAlpha="dynamic",
    nMinIVs=50,
    hessMinScale=5,
    silence=False,
    UniMRIVPThreshold=1e-8,
    candidateGenes="",
    assumeNoSampleOverlap=True,
    shrinkBiasCorrection=True,
    networkR2Thres=0.25,
    impute=True,
    saveData=False,
    geneSelector="mrjones",
):
    # ldOtherLociOtherPt: only SNPs with P<this threshold will be considered when removing SNPs in IV set that are in LD with SNPs outside of the LD set
    # making sure the user hasn't passed anything crazy for ldOtherLociWindow
    if (ldOtherLociWindow > 10) & (silence == False):
        raise ValueError(
            "are you sure you gave `ldOtherLociWindow` the right value? This is {} megabases! If so, set the `-silence` flag to `true` to silence this message next time".format(
                str(ldOtherLociWindow * 1e6)
            )
        )
    # merged: full data set of merged gene expression and phenotype GWAS data - already went through QC (e.g., harmonizing alleles)
    # geneGroups: dictionary with list elements of genes that belong together in groups
    # ggKeys: list(geneGroups) - ie names of keys of geneGroups dictionary
    # writableDir: a directory that the user can write in - usually the current WD from which this software is executed
    # ldRefDir: filepath to LD reference panel (without extension - must have .bed, .fam, bim accompanying files)
    # analysisOnlyInOutcomeLoci: boolean, if true, will only perform MR within/near outcome clumps
    # outcomeLociMbWindow: float; if analysisOnlyInOutcomeLoci==True, will only perform MR for groups with any genes within outcomeLociMbWindow Mb of outcome an clump
    # ldUpperLimit: for 2 SNPs with LD>this threshold, one will be randomly dropped
    # q0Correls: Only SNPs with P<this quantile will be considered when estimating correlations between measurement errors
    # nMinCorrels: If <this number of nonsig and nonmissing SNPs are avaiable to calcualte correls b/w measurement errors, set as missing
    # jointChiGenesP: only SNPs with P<this threshold in a p-degree of freedom joint test for p genes will be included as IVs
    # assumedMissingMean: assume that SNPs with missing association estimates (Z-stats) have this value as their mean
    # opAlpha: alpha to use in MR-Jones (note no grid searching is currently allowed so as to keep the program fast), or 'dynamic' (chosen from correlations; <1)
    # verbose: True/False - should progress be printed
    # nMinIVs: minimum number of IVs to use in MR; else group will be skipped
    # hessMinScale: diag(bx.T@LDinv@bx)/diag(m*Suu) must have be greater than hessMinScale
    # silence: a warning if it seems like the user set the ldOtherLociWindow to be extremely large
    # UniMRIVPThreshold: P-value threshold for univariable MR using the perGeneSMRBias() function
    ###
    ### want to record positions of outcome clumps up front
    # merged['Gene']=merged['Gene'].str.split(pat='.',n=0,expand=True).iloc[:,0] # dropping the '.<x>' suffixes from gene IDs (applies to GTEx, for example)
    outcomeclumpbps = merged[merged["isOutcomeClump"] == True][
        "snpBP"
    ].unique()  # workhorse() receives CHR-specific data, so just need BP position
    ### first, if user just wants to test a set of candidate genes, make sure this CHR has some of them, else next
    if len(candidateGenes) > 0:  # redundant, but just in case
        candidateGenes = [
            candidateGenes[i].split(".")[0] for i in range(0, len(candidateGenes))
        ]
        analysisOnlyInOutcomeLoci = False  # perform analysis in candidateGenes loci
        nMinIVs = 10
        hessMinScale = 0  # I want to let the user get what they ask for, instead of restrciting them
    progs = list(numpy.linspace(0, len(geneGroups), 10))
    progs = [int(progs[x]) for x in range(0, len(progs))]  # for progress printing
    thingsMonitored = {}
    outerDict = {}
    edgeDict = {}
    for ogene in range(0, len(geneGroups)):  # 0, len(geneGroups)
        thingsToMonitor = {}
        ############################################################### start top 1/3rd (preparing data) (avg ~7 seconds per gene group)
        # 1) subsetting full data to only selected genes
        theseGenes = geneGroups[ggKeys[ogene]]  # identify which genes are in this group
        if isinstance(theseGenes, str):
            theseGenes = [theseGenes]  # if only one gene considered
        # if user only wants to test candidate genes, restrict to only those groups containing at least one of the candidate genes
        if len(candidateGenes) > 0:
            tt = numpy.array(
                [theseGenes[x] in candidateGenes for x in range(0, len(theseGenes))]
            )
            if sum(tt) == 0:
                continue

        thingsToMonitor["startingNumGenes"] = len(
            theseGenes
        )  # number of starting genes
        genedf = merged[
            merged["Gene"].isin(theseGenes)
        ]  # filter full data set to only these genes
        genedf = genedf.dropna(
            subset=["phenoZ"]
        )  # drop rows with missing values for outcome (I will not want to impute these)
        genedf = genedf.drop_duplicates(
            subset=["geneSNP", "Gene"], keep="first"
        )  # cannot keep SNPs tested >once for a gene
        sn = list(
            genedf["geneSNP"].unique()
        )  # list of unique SNPs present in data for these genes (will be number of rows in wide data matrix)
        ### if the user wants to only look in significant outcome loci (see arguments analysisOnlyInOutcomeLoci and outcomeLociMbWindow)
        if analysisOnlyInOutcomeLoci == True:
            sbps = genedf["snpBP"].unique()
            tt = [
                any(abs(sbps - outcomeclumpbps[x]) < (outcomeLociMbWindow * 1e6))
                for x in range(0, len(outcomeclumpbps))
            ]
            if (
                sum(tt) == 0
            ):  # if no genedf SNPs (candidate IVs) are within outcomeLociMbWindow of any outcome clumps, move to the next group
                thingsToMonitor["locusSkippedBCNoOutcomeSignal"] = ggKeys[
                    ogene
                ]  # defining the gene locus using the index gene
                continue

        ### may still want to impose a restriction on which SNPs can be included - some SNPs have >75% missing for genes - do we really want to keep these?
        gnKeep = genedf["Gene"].value_counts()
        gnKeep = gnKeep[
            gnKeep > numpy.min([100, len(sn) * 0.1])
        ]  # restricting genes to only those with >100 SNPs or at least 10% observed
        genedf = genedf[genedf["Gene"].isin(list(gnKeep.index))]
        # restricting and noting which genes got removed
        newGenes = genedf[
            "Gene"
        ].unique()  # list of theseGenes cut to those with a enough observed data to continue
        lostGenes = [
            theseGenes[x]
            for x in range(0, len(theseGenes))
            if (theseGenes[x] in newGenes) == False
        ]  # which genes got removed?
        thingsToMonitor["droppedGenes"] = lostGenes
        if len(newGenes) == 1:  # if all but one gene just got removed
            reshaped = genedf.copy()
            reshaped = reshaped.rename(columns={"geneZ": newGenes[0] + "_Z"})
        else:
            reshaped = reshapeCorrelDf(
                genedf, newGenes, "Gene", joinType="outer"
            )  # reshaping data to be wide format (will have some NaN's if joinType='outer')

        ### identify lead SNP for each gene (its ok if some of these will not be included in MR, this is for later infererence and plotting)
        genedf["absZ"] = numpy.abs(genedf["geneZ"])
        leadGeneSNPDf = genedf.loc[genedf.groupby("Gene")["absZ"].idxmax()][
            ["geneSNP", "Gene", "geneZ"]
        ]
        genedf = genedf.drop(columns=["absZ"], axis=1)
        cn = [
            newGenes[x] + "_Z" for x in range(0, len(newGenes))
        ]  # column names for Z-statistics
        thingsToMonitor["genesDroppedBCMissing"] = len(theseGenes) - len(
            cn
        )  # how many genes got dropped because of missing?
        thingsToMonitor["startingNSNPs"] = reshaped.shape[0]  # starting number of SNPs
        geneLocs = genedf[
            ["Gene", "geneBP"]
        ].drop_duplicates()  # record locations of all genes in this group now before any get dropped
        ############################################################### end top 1/3rd (preparing data) (avg ~7 seconds per gene group)
        ############################################################### start middle 1/3rd (imputation, measurement error correlations, IV selection, LD estimation) (avg ~7 seconds per gene group)
        ### Imputation of missing in actual gene exposure data matrix
        # Yihe's reduced rank matrix completion (imputation with 0s or GAM predictions can both be problematic)
        X = reshaped[cn].values
        # actually a lot of missing if the number of genes in a group is very large
        if len(cn) == 1:  # if only one gene, still want this to be a matrix
            X = X.reshape((X.shape[0], 1))
        # cn and columns of X are in the same order. How to find geneBPs that are in this same order:
        cnstripped = [cn[x].split("_", 1)[0] for x in range(0, len(cn))]
        gg = pandas.merge(
            pandas.DataFrame({"hi": cnstripped}),
            geneLocs,
            left_on="hi",
            right_on="Gene",
        )["geneBP"].values
        mask = numpy.isnan(X)
        Del = mask + 1 - 1
        newX = X.copy()
        if (X.shape[1] > 1) & (impute == True):
            L, Factor, Gamma = ffimputeBIC(
                X,
                Kmax=15,
                dg=gg,
                dim_proj=0,
                W=makeI(X.shape[0]),
                iter_max=15,
                eps=0.001,
                tuning=0,
            )
            newX[mask] = L[mask]
        elif (X.shape[1] > 1) & (impute == False):
            newX[mask] = 0

        # Omega=numpy.isnan(X); Xc=X.copy(); Xc[Omega]=assumedMissingMean; d0=numpy.linalg.svd(Xc,compute_uv=False)
        # newX,out2,out3=fullRankMatrixCompletion(X.copy(),missImp=0,lambda_=min(d0)/2,eps=0.01*X.shape[1]) # lambda_ chosen somewhat arbitrarily
        mm = numpy.isnan(X).sum(axis=0) / X.shape[0]
        thingsToMonitor["minPropMissingZs"] = numpy.min(mm)
        thingsToMonitor["maxPropMissingZs"] = numpy.max(mm)
        thingsToMonitor["minNumNonMissing"] = X.shape[0] - numpy.max(
            numpy.isnan(X).sum(axis=0)
        )
        # probably need to do some more work on the method in fullRankMatrixCompletion(), but I think generally the idea is appropriate.
        ### calculation of R=DSD where S=SigmaUU
        ccDf = reshaped[cn]
        # `reshaped` with only columns with names in `cn`. Must reset index
        copycat = (ccDf.copy()).reset_index()
        ccDf = ccDf.assign(
            myInd=list(copycat.index)
        )  # add column to end of ccDf indicating row number (indices get all messed up)
        cond = numpy.array(abs(ccDf) > q0Correls)
        # for each pair, find all SNPs where cond==False for both, which may ==True for others
        p = ccDf.shape[1] - 1
        cc = makeI(p)
        nn = numpy.zeros((p, p))
        for i in range(0, p):
            for j in range(0, p):
                if i > j or i == j:
                    continue
                else:
                    condi = flatten_list(cond[:, i])
                    condj = flatten_list(
                        cond[:, j]
                    )  # finding logical for column i and j
                    keepijs = [
                        x
                        for x in range(0, len(condi))
                        if (condi[x] == False and condj[x] == False)
                    ]  # want where both nonsignificant
                    subCCDF = ccDf.copy()
                    subCCDF = subCCDF.iloc[
                        keepijs, [i, j, p]
                    ].dropna()  # and remove missing; p+1 is to keep the row number column I made
                    mask = numpy.zeros(len(condi), dtype="bool")
                    mask[subCCDF["myInd"].values] = True
                    subCCDF = subCCDF.drop("myInd", axis=1)
                    nn[i, j] = subCCDF.shape[0]
                    if (
                        nn[i, j] < nMinCorrels
                    ):  # I previously chose 50 somehwat arbitrarily
                        # print('problem with '+str(i)+' and '+str(j)) # probably want to comment this out
                        cc[i, j] = float(
                            "nan"
                        )  # make missing - to be filled in with MaxDet imputation
                        continue
                    # correlation using only nonsig
                    # recall indexing by rows and columns must happen sequentially not simultaneously
                    cc[i, j] = numpy.corrcoef(subCCDF.values, rowvar=False)[
                        0, 1
                    ]  # no need to premultiply by LD matrix bc since nonsig assume not affected by LD
                    # correlation using all with resampling
                    # cc[i,j]=rs,vstotal=subSamCalcR(subCCDF,k=subCCDF.shape[0])

        cc = cc + cc.T
        numpy.fill_diagonal(cc, 1)
        del ccDf
        Omega = numpy.isnan(cc)
        thingsToMonitor["numMissingMECorrelations"] = (
            sum(Omega.reshape((Omega.shape[0] * Omega.shape[1],))) / 2
        )  # number of missing values in correlation matrix
        ### imputation of missing in R=DSD where S=SigmaUU
        # using MaxDet to fill in missing & drop any genes that you need to (where MaxDet sol not in {-1,-1} (if no missing this code can still run and nothing will happen)
        ccComplete, ccDropped, ccNDropped = MaxDet(cc)
        mask = numpy.ones(
            (cc.shape[0],), dtype="bool"
        )  # creating a mask to drop genes/exposures (see below)
        mask[ccDropped] = (
            False  # if gene was dropped bc of a bad MaxDet solution, do not consider it below
        )
        Omega = Omega[mask, :][
            :, mask
        ]  # find indices of missing for non-bad MaxDet solutions (ie, find which were successfully imputed)
        cnMask = numpy.ones(len(cn), dtype="bool")
        cnMask[ccDropped] = False
        newX = newX[:, mask]
        cn = numpy.array(cn)[cnMask].tolist()
        thingsToMonitor["numBadMaxDetSolutions"] = (
            cc.shape[0] - ccComplete.shape[0]
        )  # number of genes dropped bc of bad MaxDet solutions
        # do not drop any SNPs - this may change in the future, but for now do not
        ccComplete, alpha, deco = posDefifyCorrMat(
            ccComplete, epsilon=1e-5
        )  # making ccComplete positive definite via shrinkage if it is not
        thingsToMonitor["alphaToMakeMECorrelMatPosDef"] = (
            alpha  # factor that off-diagonal elements of ccComplete needed to multiplied by to make it posdef
        )
        # remove genes with norm(newX[:,j],2)^2<5m
        toKeep = numpy.diag(newX.T @ newX) > (
            hessMinScale * newX.shape[0] / 2
        )  # technically >hessMinScale*newX.shape[0]*diag(Suu), but diag(Suu)=(1)
        if sum(toKeep) == 0:
            continue
        newX = newX[:, toKeep]
        ccComplete = ccComplete[toKeep, :][:, toKeep]
        droppedGenes = numpy.array(cn)[toKeep == False].tolist()
        cn = numpy.array(cn)[toKeep].tolist()
        reshaped = reshaped.drop(labels=droppedGenes, axis=1)
        # perform p-degree of freedom chi-square test to test SNP assoc with at least one gene at P<tau. drop SNPs not meeting this criteria
        # Note I need SigmaUU, but not the LD matrix, for this so it happens here and not earlier
        tau = jointChiGenesP
        # no need to pre-multiply by sparsePrecisionSqrt bc I am considering one SNP at a time
        chi_ = numpy.diag(
            newX @ numpy.linalg.inv(ccComplete.copy()) @ newX.T
        )  # the inverse matrix is the inverse of SigmaUU
        chiPs = [
            1 - stats.chi2.cdf(chi_[x], newX.shape[1]) for x in range(0, newX.shape[0])
        ]
        mask = numpy.array(chiPs) < tau
        mask = mask.reshape((mask.shape[0],))
        reshaped = reshaped[mask]
        newX = newX[mask, :]
        if newX.shape[0] < nMinIVs:
            continue  # also need to go to next if not many SNPs here

        # 3) retrieving LD matrix for selected genes
        # `reshaped` contains the SNPs I want to use in MR (`ccDf` was only used for calculating the correlation matrix)
        # generating SNP LD matrix
        reshaped = reshaped.reset_index()
        reshaped = reshaped.drop_duplicates(
            subset=["geneSNP"], keep=False
        )  # drop duplicate SNPs (cannot appear in file given to PLINK)
        mask = numpy.zeros((newX.shape[0]), dtype="bool")
        mask[reshaped.index.tolist()] = True
        newX = newX[mask, :]  # in case any duplicates got removed
        outDir1 = (
            writableDir + "/myExtract.txt"
        )  # (for saving SNPs for LD calculations) write out to a directory that must be writable (since the user is there)
        outDir2 = writableDir + "/tempOut"  # (for saving LD matrix)
        reshaped["index"] = list(
            range(0, reshaped.shape[0])
        )  # this is the current ordering of rows of newX
        reshaped = reshaped.sort_values(
            by="snpBP", ascending=True
        )  # need to order by BP so I know order of SNPs in ldmat from PLINK
        newX = newX[reshaped["index"].values]
        mstart = newX.shape[
            0
        ]  # Need to make sure the data I already created (newX/bX) matches this order
        minBP = min(reshaped["snpBP"].values)
        maxBP = max(reshaped["snpBP"].values)
        thingsToMonitor["BPSpanOfIVSet"] = maxBP / 1e6 - minBP / 1e6

        if ldOtherLociWindow > 0:
            # defaults: ldOtherLociOtherPt=5e-5,ldOtherLociR2=0.1
            mask1 = merged[
                (merged["snpBP"] < minBP)
                & (merged["snpBP"] > (minBP - ldOtherLociWindow * 1e6))
            ]
            mask2 = merged[
                (merged["snpBP"] > maxBP)
                & (merged["snpBP"] < (maxBP + ldOtherLociWindow * 1e6))
            ]
            mask1 = mask1[abs(mask1["geneZ"]) > norm.ppf(1 - ldOtherLociOtherPt)]
            mask2 = mask2[abs(mask2["geneZ"]) > norm.ppf(1 - ldOtherLociOtherPt)]
            mask1 = (
                mask1[["geneSNP", "snpBP"]]
                .drop_duplicates(keep="first")
                .sort_values(by="snpBP", ascending=True)
            )
            mask2 = (
                mask2[["geneSNP", "snpBP"]]
                .drop_duplicates(keep="first")
                .sort_values(by="snpBP", ascending=True)
            )
            # write out big list in this order: mask1, reshaped, mask2
            big_list = pandas.concat(
                [mask1["geneSNP"], reshaped["geneSNP"], mask2["geneSNP"]]
            )
            big_list.to_csv(
                outDir1, header=False, index=False, sep=" "
            )  # writing out list of SNPs to give PLINK to make LD calculations
            cmd = [
                callPlink(),
                "--r",
                "square",
                "--bfile",
                ldRefDir,
                "--extract",
                outDir1,
                "--out",
                outDir2,
            ]  # the command for PLINK
            out = subprocess.call(
                cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )  # executing the command in the terminal
            ldmat = numpy.genfromtxt(
                outDir2 + ".ld", dtype="f"
            )  # reading LD matrix back in
            # now need to know if any SNPs in reshaped are in LD>ldOtherLociR2(=0.1 be default) with SNPs outside of reshaped
            maskreshaped = numpy.ones((ldmat.shape[0],), dtype="bool")
            maskreshaped[: mask1.shape[0]] = False
            maskreshaped[(mask1.shape[0] + reshaped.shape[0]) :] = False
            inds1 = ldmat[mask1.shape[0] : (mask1.shape[0] + reshaped.shape[0]), :][
                :, : mask1.shape[0]
            ]
            inds1 = numpy.nan_to_num(inds1, nan=0)
            inds2 = ldmat[mask1.shape[0] : (mask1.shape[0] + reshaped.shape[0]), :][
                :, (mask1.shape[0] + reshaped.shape[0]) :
            ]
            inds2 = numpy.nan_to_num(inds2, nan=0)
            ldmask1 = (((inds1**2) > ldOtherLociR2).sum(axis=1)) > 0
            ldmask2 = (((inds2**2) > ldOtherLociR2).sum(axis=1)) > 0
            ldmask = (ldmask1 + ldmask2) == 0
            reshaped = reshaped[numpy.array(ldmask).astype(bool)]
            newX = newX[ldmask]
            ldmat = ldmat[maskreshaped, :][
                :, maskreshaped
            ]  # subset to only SNPs in reshaped
            ldmat = ldmat[ldmask, :][
                :, ldmask
            ]  # subset again to only SNPs in reshaped uncorrelated with SNPs in surrounding window
            # some LD correlations are totally missing for a SNP (rare, but possible). remove these. (possible why? observed AF=0 or 1 so no variation?)
            mask = (
                (numpy.isnan(ldmat).sum(axis=0)) == ldmat.shape[0]
            ) == False  # create mask: True if non-missing LD; False if all missing LD (I've never seen partial missing)
            reshaped = reshaped[mask]
            newX = newX[mask, :]
            ldmat = ldmat[mask, :][:, mask]  # drop SNPs with missing LD
            del mask1, mask2, ldmask, inds1, inds2, ldmask1, ldmask2
        else:
            # personal check (should be all True's): [sum(newX[:,xx]==reshaped[cn[xx]].values)==newX.shape[0] for xx in range(0,newX.shape[1])]
            reshaped["geneSNP"].to_csv(
                outDir1, header=False, index=False, sep=" "
            )  # writing out list of SNPs to give PLINK to make LD calculations
            cmd = [
                callPlink(),
                "--r",
                "square",
                "--bfile",
                ldRefDir,
                "--extract",
                outDir1,
                "--out",
                outDir2,
            ]  # the command for PLINK
            out = subprocess.call(
                cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )  # executing the command in the terminal
            ldmat = numpy.genfromtxt(
                outDir2 + ".ld", dtype="f"
            )  # reading LD matrix back in
            # ldmat is ordered by BP position, NOT by order of SNPs in myExtract.txt. That is why I sorted reshaped by snpBP
            # some LD correlations are totally missing for a SNP (rare, but possible). remove these. (possible why? observed AF=0 or 1 so no variation?)
            mask = (
                (numpy.isnan(ldmat).sum(axis=0)) == ldmat.shape[0]
            ) == False  # create mask: True if non-missing LD; False if all missing LD (I've never seen partial missing)
            reshaped = reshaped[mask]
            newX = newX[mask, :]
            ldmat = ldmat[mask, :][:, mask]  # drop SNPs with missing LD

        thingsToMonitor["snpsDroppedBCCHPCorrectionByLD"] = mstart - newX.shape[0]
        ldc = abs(ldmat.copy())
        ss = (ldc > ldUpperLimit).sum(axis=0).sum(axis=0)
        if (ldmat.shape[0] < nMinIVs) or (
            ss == ldmat.shape[0] ** 2
        ):  # if not enough SNPs or all are in near-perfect LD with each other
            continue
        thingsToMonitor["snpsDroppedBCMissingLD"] = sum(
            mask == False
        )  # save what happened here
        snps = reshaped["geneSNP"].values  # find SNPs currently considered
        someDropped, snpsKept = LDShrinkUpper(ldmat, ldUpperLimit, snps)
        ld0 = (
            someDropped.copy()
        )  # find those not in ~perfect LD (>tau) with any others (makes ldmat pos def)
        mask = reshaped["geneSNP"].isin(snpsKept)
        reshaped = reshaped[mask]  # restrict to only SNPs from line above
        newX = newX[mask]
        thingsToMonitor["snpsDroppedBCHighLD"] = (
            ldmat.shape[0] - ld0.shape[0]
        )  # number of SNPs dropped because of high LD correlations
        ######################### if the number of IVs is very small, we can't really trust any results
        if ld0.shape[0] < nMinIVs:
            continue
        #########################
        # R=ld0.copy();kTaper=int(ld0.shape[0]/4);kBand=ld0.shape[0]+1;lambda_=0.01;epsilon=1e-4;toCorr=True
        ld0, alpha, deco = regularizeLD(
            ld0.copy(),
            int(ld0.shape[0] / 4),
            ld0.shape[0] + 1,
            0.1,
            epsilon=1e-4,
            toCorr=True,
            justMCP=False,
        )  # int rounds down
        w, v = deco[0], deco[1]
        sparsePrecision = v @ numpy.diag(1 / w) @ v.T
        sparsePrecisionSqrt = v @ numpy.diag(1 / w**0.5) @ v.T
        thingsToMonitor["alphaToMakeLDMatPosDef"] = alpha
        ############################################################### end middle 1/3rd (correlations and imputation)

        ############################################################### start final 1/3rd (MR and saving data) (avg ~2 seconds per gene group)
        # organizing data
        bx = numpy.array(newX)
        by = reshaped["phenoZ"].values
        bx0 = bx.copy()
        by0 = by.copy()
        p = 1 if len(bx.shape) == 1 else bx.shape[1]
        q = 1 if len(by.shape) == 1 else by.shape[1]
        m = bx.shape[0]
        bx = bx.reshape((m, p))
        by = by.reshape((m, q))
        bx0 = bx0.reshape((m, p))
        by0 = by0.reshape((m, q))
        thingsToMonitor["finalNumIVs"] = bx.shape[0]
        thingsToMonitor["outcomeEstH2InLocus"] = (by.T @ by) - by.shape[0]
        ##########################################
        # remove LD structure at the outset - seems to cause problems in standard error estimation - I currently have not found the problem or another solution
        bx = sparsePrecisionSqrt @ bx
        by = sparsePrecisionSqrt @ by
        # if isGtex: # if GTEx, all effect estimates are already joint
        #     bx=bx0.copy()
        #     by=by0.copy()
        #     ld0=numpy.eye(m)
        #     sparsePrecision=ld0.copy()
        #     sparsePrecisionSqrt=ld0.copy()
        ld0og = ld0.copy()
        ld0 = numpy.eye(bx.shape[0])
        sparsePrecision = ld0.copy()
        sparsePrecisionSqrt = ld0.copy()
        # this has implications for how nclumps are estimated - that's why I keep a copy of the original LD matrix
        ##########################################
        ########################################## (yes I am doing this for a second time - works well in both cases)
        # remove genes with norm(newX[:,j],2)^2<5m
        toKeep = numpy.diag(bx.T @ bx) > (
            hessMinScale * bx.shape[0] / 2
        )  # can change factor if you want; default is hessMinScale=5
        if sum(toKeep) < 2:  # if less than 2 genes available
            continue
        bx = bx[:, toKeep]
        bx0 = bx0[:, toKeep]
        ccComplete = ccComplete[toKeep, :][:, toKeep]
        droppedGenes = numpy.array(cn)[toKeep == False].tolist()
        cn = numpy.array(cn)[toKeep].tolist()
        reshaped = reshaped.drop(labels=droppedGenes, axis=1)
        p = sum(toKeep)
        ##########################################
        # ells=ldscore(ld0); ells=numpy.column_stack(([1]*ells.shape[0],ells)) # only approximate LD scores (only considers LD in ld0)
        # (numpy.linalg.inv(ells.T@ells)@ells.T@bx).shape
        # variance-covariance matrices of measurement errors
        VV = numpy.array(1).reshape((1, 1))
        UU = (
            ccComplete.copy()
        )  # maybe not memory friendly to make a copy but both object names work well in context
        # now to calculate UV
        if assumeNoSampleOverlap:  # if user wants to assume there is no sample overlap
            UV = numpy.zeros((p, 1))
        else:
            bybx = numpy.hstack(
                (by, bx)
            )  # already no missing because I removed it from newX (bx)
            mask = (
                abs(bybx > q0Correls).sum(axis=1) == 0
            )  # True values in mask are the SNPs to use in correls calculations
            mask = numpy.where(mask)[
                0
            ]  # for some reason returns a tuple so I add the [0]
            if mask.shape[0] < nMinCorrels:
                uv0 = numpy.zeros((p, q)) * float("nan")
                a_ = numpy.hstack((VV, uv0.T))
                b_ = numpy.hstack((uv0, UU))
                c_ = numpy.vstack((a_, b_))
                md = MaxDet(c_, doToBadSolutions="drop")[0]
                mask = abs(md) > 1
                md[mask] = 0
                UV = md.copy()
                UV = UV[:, 0]
                UV = UV[1:]
                UV = UV.reshape((p, 1))
            else:
                UV = numpy.corrcoef(
                    sparsePrecisionSqrt[mask, :][:, mask] @ bybx[mask, :], rowvar=False
                )[
                    :, 0
                ]  # Corr(by,bx) will be in first column
                UV = numpy.delete(UV, 0)
                UV = UV.reshape((p, 1))  # drop the 1 for Corr(by,by)

        # add intercept-relevant terms for MRJones
        bx_ = numpy.column_stack(([1] * m, bx))
        bx0_ = numpy.column_stack(([1] * m, bx0))
        UU_ = numpy.column_stack(([0] * p, UU))
        UU_ = numpy.row_stack(([0] * (p + 1), UU_))
        UV_ = numpy.row_stack((0, UV))
        ### shrinking of UU and UV to make sure they don't introduce big problems in causal estimation
        if shrinkBiasCorrection:
            UU = UU * 0.5
            UU_ = UU_ * 0.5
            UV = UV * 0.5
            UV_ = UV_ * 0.5
        if (opAlpha == "dynamic") & (bx.shape[0] > 10):
            c_c = numpy.corrcoef(bx, rowvar=False)
            numpy.fill_diagonal(c_c.reshape((p, p)), 0)
            opAlphaVal = float(numpy.max(abs(c_c)))
        elif is_number(opAlpha):
            opAlphaVal = opAlpha
        else:
            opAlphaVal = 0
        ## first defining hyperparameter lambda1, lambda2/tau grids
        # does the user want to use MR-Jones or G-Screen to perform gene selection?
        lamvec = numpy.linspace(0.01, 3 / 4, 20)
        if geneSelector == "mrjones":
            t1 = numpy.linspace(0.01, 1 / 2, 10)
            t2 = numpy.linspace(1, 3, 10)
            t3 = numpy.linspace(5, 10, 10)
            t4 = numpy.linspace(20, 40, 20)
            tauvec = numpy.concatenate((t1, t2, t3, t4))
            out = MRJones(
                by0,
                bx0_,
                ld0og,
                UU_,
                UV_,
                lamvec=lamvec,
                tauvec=tauvec,
                rho_theta=3 / 2,
            )
            deltas = out["gamma"]
            finalEsts = out["theta"]
            thingsToMonitor["propNonzeroDeltas"] = sum(deltas != 0) / m
        else:  # G-screen
            finalEsts, w0, bic = gscreen(
                bx,
                by,
                UU,
                gamma=mad(by.squeeze()) * 2,
                lamvec=lamvec,
                rho=3 / 2,
                eps=1e-2,
                max_iter=15,
            )
            thingsToMonitor["propNonzeroDeltas"] = 1 - sum(w0) / m
            finalEsts = finalEsts.reshape(
                bx.shape[1],
            )
            finalEsts = numpy.concatenate(
                [[0], finalEsts]
            )  # needs an intercept, which I forced to 0 in the gscreen() function
        n = reshaped["geneN"].values
        n = numpy.nanmedian(n)
        # variance explained
        if all(finalEsts == 0):
            r2mr = 0
        else:
            em = numpy.column_stack((by0, bx0_ @ finalEsts))
            if numpy.all(em[:, 0] == em[0, 0]) | numpy.all(em[:, 1] == em[1, 1]):
                r2mr = 0
            else:
                (numpy.sum(em[:, 1] ** 2) - m * numpy.sum(finalEsts**2)) / (
                    numpy.sum(by0**2) - m
                )
                r2mr = numpy.corrcoef(em, rowvar=False)[0, 1] ** 2
        finalEsts = finalEsts[1:]
        paraDf = pandas.DataFrame.from_dict(
            {"gg": numpy.array(cn), "GeneSelected": finalEsts != 0}
        )
        # add lead SNPs for each gene to use in later plotting
        leadGeneSNPDf["Gene"] = leadGeneSNPDf["Gene"] + "_Z"  # for merging
        paraDf = pandas.merge(
            paraDf, leadGeneSNPDf, left_on="gg", right_on="Gene"
        )  # merge
        # conditional F-statistics (MRJones recommended for this; lamvec determined automatically from within the function)
        if bx.shape[1] < 2:
            condFs = float("nan")
        else:
            condFs = conditionalFMRBEE(
                bx,
                UU,
                sparsePrecisionSqrt,
                tauLowerQuant=0.75,
                tauUpperQuant=0.95,
                opAlphaVal=opAlphaVal,
                method="mrbee",
            )  # 'notfree': MRBEE; 'free': MR-Jones
        paraDf["conditionalF"] = condFs
        # adding gene BP to results
        gl = []
        cn_ = [
            cn[x].split("_", 1)[0] for x in range(0, len(cn))
        ]  # cn_ drops '_Z' from cn (all cn==paraDf['gg'])
        for _ in range(0, len(cn_)):
            mask = cn_[_] == geneLocs["Gene"].values
            ta = geneLocs[mask]["geneBP"].values
            gl.append(float(ta))
        # finally, perform analysis using MVMR and uni MR with MR-Egger (this is equivalent to SMR and the Porcu et al methods)
        ivwout = IVW(
            bx, by, sparsePrecision, newX, uniIVPThreshold=1e-3
        )  # only IVs with P<uniIVPThreshold will be considered in univariable MR
        ivwout["gg"] = (
            cn  # cn is the order of genes in paraDf (check: [cn[x]==paraDf['gg'].values[x] for x in range(0,len(cn))])
        )
        paraDf = pandas.merge(paraDf, ivwout, left_on="gg", right_on="gg")
        # add BP position for lead SNP to paraDf
        paraDf = pandas.merge(
            paraDf,
            genedf[["geneSNP", "snpBP"]].drop_duplicates(keep="first"),
            left_on="geneSNP",
            right_on="geneSNP",
        )  # add snpBP to paraDf for lead SNP
        paraDf["Gene"] = [
            paraDf["Gene"][_].split("_", 1)[0] for _ in range(0, paraDf.shape[0])
        ]  # for merging
        paraDf = pandas.merge(
            paraDf,
            genedf[["Gene", "geneBP"]].drop_duplicates(keep="first"),
            left_on="Gene",
            right_on="Gene",
        )  # add gene BP to paraDf for each gene
        paraDf = paraDf.sort_values("geneBP")  # order by gene BP
        # approximate bias in using SMR
        smrbias = perGeneSMRBias(
            bx,
            by,
            ld0,
            sparsePrecision,
            newX,
            UU,
            UV,
            cn_,
            IVPThresh=UniMRIVPThreshold,
            r2Thresh=0.25,
            minNIVs=2,
        )
        paraDf = pandas.merge(paraDf, smrbias, left_on="Gene", right_on="gene")
        # # MRBEE with ridge
        # mrbeeRidge=MRBEERidge(bx,by,UU,UV,VV,sparsePrecision,nLambda=20)
        # ridgeDf=pandas.DataFrame({'Gene': cn_, 'MRBEERidgeEst': mrbeeRidge[0].squeeze(), 'MRBEERidgeSE': (numpy.diag(mrbeeRidge[1])**0.5).squeeze()})
        # paraDf=pandas.merge(paraDf,ridgeDf,left_on='Gene',right_on='Gene')
        ### MRBEE Post-Selection
        mask = finalEsts != 0  # True where gene should be included, False where not
        if sum(mask) == 0:
            toa = numpy.array([0] * len(mask)) * float("nan")  # vector of NaN's
            toAdd = pandas.DataFrame(
                {"Gene": cn_, "MRBEEPosSelEst": toa, "MRBEEPosSelSE": toa}
            )
            paraDf = pandas.merge(paraDf, toAdd, left_on="Gene", right_on="Gene")
        else:
            bxkeep = bx[:, mask]
            UUkeep = UU[mask, :][:, mask]
            UVkeep = UV[mask]
            p1, p2, p3, p4, p5 = imrbee(
                bxkeep,
                by,
                UUkeep,
                UVkeep,
                VV,
                sparsePrecision,
                ld0,
                0.05 / bxkeep.shape[0] ** 0.5,
                boot=False,
            )
            p1 = p1[1:]
            p2 = p2[1:, :][:, 1:]  # remove intercept-relevant terms
            toAdd = pandas.DataFrame(
                {
                    "Gene": numpy.array(cn_)[mask].tolist(),
                    "MRBEEPosSelEst": p1.squeeze(),
                    "MRBEEPosSelSE": numpy.diag(p2) ** 0.5,
                }
            )
            paraDf = pandas.merge(
                paraDf, toAdd, how="outer", left_on="Gene", right_on="Gene"
            )
        # single SNP
        singlesnpdf = singleSNP(bx, by, cn_)
        paraDf = pandas.merge(
            paraDf, singlesnpdf, how="outer", left_on="Gene", right_on="Gene"
        )
        # estimate gene expression heritabilities and numbers of SNPs with nonzero associations with gene expression (from theory)
        # basically performing clumping without window size restriction within the SNP set now: for each gene, number of P(Z)<0.05 and independent with other SNPs
        nclumps = []
        r2 = 0.01  # like r2 in PLINK
        ld0 = (
            ld0og.copy()
        )  # since all MR calculations are done, resetting LD matrix to be actual LD matrix (ie not identity matrix bc I made a transformation)
        for p_ in range(0, bx.shape[1]):
            mask = numpy.abs(bx0[:, p_]) > norm.ppf(
                1 - 1e-5
            )  # (1e-5)th quantile; only considered LD between P(Z)<0.05 SNPs
            ld0mask = ld0[mask, :][:, mask] ** 2  # squaring all values for convenience
            ninds = []
            for s_ in range(0, ld0mask.shape[0]):
                m1 = numpy.zeros((ld0mask.shape[0],), dtype="bool")
                m1[s_] = True  # creating mask to index this specific SNP
                ldm1 = ld0mask[m1, :][:, m1 == False].squeeze()
                ss = numpy.sum(
                    ldm1 > (r2)
                )  # SNPs with squared LD < r2 will be considered independent
                if ss == 0:  # if no other SNPs in LD with this one
                    ninds.append(1)  # put 1 (I will sum up later)
            if len(ninds) == 0:
                nclumps.append(
                    1
                )  # if can't find any for this gene, just put 1 - assume there is at least 1
            else:
                nclumps.append(sum(ninds))  # else, put total number I counted
        # nclumps is an estimate of the number of independent causal SNPs in this locus
        nclumps = numpy.array(nclumps)
        eps = (
            numpy.array(nclumps) / bx.shape[0]
        )  # an estimate of the proportion of causal SNPs in this locus
        medns = reshaped["geneN"].values
        N = numpy.diag(1 / medns)
        medns = numpy.nanmedian(medns)
        cisEstNaive = (
            ((numpy.diag(bx.T @ bx) - 1 / medns) / bx.shape[0] - 1) * nclumps / medns
        )  # medns is the median sample sizes for each gene
        toAdd = pandas.DataFrame(
            {"Gene": cn_, "nClumps": nclumps, "h2CisEstNaive": cisEstNaive}
        )
        paraDf = pandas.merge(paraDf, toAdd, left_on="Gene", right_on="Gene")
        paraDf["LocusR2"] = r2mr
        # pretty sure h2 cannot be reliably estimated from the data because bx=bhat*nk/sigk where I do not know sigk and there is evidence that it is not 1.
        # d1['data']=reshaped # probably don't want to save the data ... will take up a lot of memory
        d1 = {}
        d1["data"] = {
            "Gene": cn_,
            "GeneBP": gl,
            "IVs": reshaped["geneSNP"].values.tolist(),
            "EffectAlleles": reshaped["phenoEffectAllele"].values.tolist(),
        }
        if (
            saveData
        ):  # currently not used after this point. data gets outputted by MVMRWorkhorse() but is not saved after that
            d1["bX"] = bx
            d1["by"] = by
            d1["UU"] = UU
            d1["UV"] = UV
            d1["regLD"] = ld0
        paraDf.index = paraDf["Gene"]
        paraDf = paraDf.reindex(cn_)
        d1["meta"] = paraDf  # put whole data frame
        # outerDict[ggKeys[ogene]]=d1
        outerDict["group" + str(ogene)] = d1
        thingsMonitored[ggKeys[ogene]] = (
            thingsToMonitor  # the key is the key for the group in geneGroups
        )
        ### MN networks
        try:
            # back-tabbed the if: else: below twice
            if (r2mr > networkR2Thres) & (bx.shape[0] > 25) & (bx.shape[1] > 1):
                XX = numpy.column_stack((by, bx))
                XXSE = numpy.ones(XX.shape)
                Rnoise = numpy.eye(bx.shape[1] + 1)
                Rnoise[1:, 1:] = UU
                Rnoise = posDefifyCorrMat(Rnoise)[0]
                lv1 = lamvec = numpy.linspace(1, 5, 10) / 100
                lv2 = lamvec = numpy.linspace(7, 12, 10) / 100
                lv3 = lamvec = numpy.linspace(15, 20, 10) / 100
                lv = numpy.concatenate((lv1, lv2, lv3))
                net = entropy_mcp_spearman_sampling(
                    XX, XXSE, Rnoise / 5, lamvec=lv
                )  # make Rnoise small
                c1 = sum(sum(numpy.isnan(net["Theta"]))) > 0
                c2 = sum(sum(numpy.isnan(net["Theta_Alt"]))) > 0
                if (c1 == True) & (c2 == False):
                    edges = (net["Theta_Alt"] != 0).astype(int)
                elif (c1 == False) & (c2 == True):
                    edges = (net["Theta"] != 0).astype(int)
                elif (c1 == False) & (c2 == False):
                    edges = (net["Theta"] != 0).astype(int)
                else:  # could not conpute network
                    edges = numpy.zeros(XX.shape)
                numpy.fill_diagonal(edges, 0)
            else:
                edges = numpy.zeros((bx.shape[1] + 1, bx.shape[1] + 1))
        except:
            edges = numpy.zeros((bx.shape[1] + 1, bx.shape[1] + 1))
        edgeDict["group" + str(ogene)] = edges
        ############################################################### final 1/3rd (MR and saving data)
    return thingsMonitored, outerDict, edgeDict
