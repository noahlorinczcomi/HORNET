# standard imports
import argparse

# custom imports


def parse_args():
    parser = argparse.ArgumentParser(
        prog="HORNET: Genome-wide robust multivariable MR with gene expression",
        description="This program performs univariable and multivariable MR with MR-Jones to identify putative causal genes and their regulatory networks.",
    )

    # ---- GWAS and eQTL DATA ----
    parser.add_argument(
        "--eQTLGWAS",
        action="store",
        type=str,
        help="(Required) path to directory/folder containing multiple files with .gz extensions that are the eQTL GWAS summary statistics for gene expression, where each file corresponds to each chromosome. There should be nothing else in this directory/folder other than these files. If you only want to perform the analysis on certain chromosomes, create a directory/folder containg the eQTL GWAS summary statistics for only those chromosomes.",
    )
    parser.add_argument(
        "--phenoGWAS",
        action="store",
        type=str,
        help="(Required) filepath with .gz extension to GWAS for outcome phenotype",
    )
    parser.add_argument(
        "--LDRef",
        action="store",
        type=str,
        default=os.getcwd() + "/data/ldref/1kg3TRANS",
        help="(Required) filepath without extension to the LD reference panel. The .bed, .bim, .fam files must each be present.",
    )
    parser.add_argument(
        "--isRawGTEx",
        action="store",
        type=str,
        default="no",
        help='(Required) If the eQTL association estimates are directly from the "eQTL Tissue-Specific All SNP Gene Associations" data sets from the GTEx Consortium, put True; else, put False.',
    )
    parser.add_argument(
        "--writableDir",
        action="store",
        type=str,
        default="tempfiles",
        help="(Optional) A directory/folder in which temporary files (eg, LD matrices) can be written and removed.",
    )

    # ---- GWAS and eQTL METADATA ----
    parser.add_argument(
        "--snpLabels",
        action="store",
        type=str,
        help="(Required,Required) comma-separated list of rsID column names for gene expression (first) and outcome phenotype (second) GWAS, respectively. If either of the data do not contain an rsID column, please put NA in the list.",
    )
    parser.add_argument(
        "--eQTLSNPBPLabel",
        action="store",
        type=str,
        help="(Required) column name of SNP base pair position values in the eQTL GWAS data sets",
    )
    parser.add_argument(
        "--zLabels",
        action="store",
        type=str,
        help="(Required,Required) comma-separated list of Z-statistic column names for gene expression (first) and outcome trait (second) GWAS, respectively. If Z-statistics are not present in one or either data set, please input BETA@SE where BETA is the name of the effect size column and SE is the name of the standard error column,",
    )
    parser.add_argument(
        "--effectAlleles",
        action="store",
        type=str,
        help="(Required, Required) column names in the eQTL GWAS (first) and phenotype GWAS (second) representing the effect alleles, respectively",
    )
    parser.add_argument(
        "--eQTLGeneLabel",
        action="store",
        type=str,
        help="(Required) the name of the column in the eQTL GWAS that indicates the gene name/label",
    )
    parser.add_argument(
        "--eQTLGeneBP",
        action="store",
        type=str,
        help="(Required) the name of the column in the eQTL GWAS that indicates the base pair location of the centers of genes or their transcription start sites",
    )
    parser.add_argument(
        "--eQTLGWASn",
        action="store",
        type=int,
        default=500,
        help="(Optional) the sample size of the eQTL GWAS",
    )
    parser.add_argument(
        "--phenotypeGWASn",
        action="store",
        type=int,
        default=30000,
        help="(Optional) the sample size of the outcome phenotype GWAS",
    )

    # --- ANALYSIS PARAMETERS ----
    parser.add_argument(
        "--candidateGenes",
        action="store",
        type=str,
        default="",
        help="(Optional) If you only want to perform the analysis for a list of candidate genes, please specify that comma-separated list here. If your list of candidate genes is large, you can place the directory to that file here. Each row of this file should contain gene IDs that are present in the eQTL GWAS data.",
    )
    parser.add_argument(
        "--analysisInPhenotypeLoci",
        action="store",
        type=str,
        default="true",
        help="(Optional) One of True or False. Should the the analyses only be performed in loci that are significantly associated with the outcome phenotype? See -phenoLociPvalue and -phenoLociMbWindow for how these loci are defined. The default is True.",
    )
    parser.add_argument(
        "--phenoLociPvalue",
        action="store",
        type=float,
        default=5e-5,
        help="(Optional) P-value threshold for considering a locus to be significantly associated with the outcome phenotype. The default is 5e-5.",
    )
    parser.add_argument(
        "--phenoLociMbWindow",
        action="store",
        type=float,
        default=0.5,
        help="(Optional) Mb window size for considering two phenotype loci as distinct. The default is 0.5.",
    )
    parser.add_argument(
        "--phenoLociR2",
        action="store",
        type=float,
        default=0.1,
        help="(Optional) SNPs associated with the outcome at significance level -phenoLociPvalue and in squared LD beyond this thresold will be clumped together using PLINKs clumping procedure. The default is 0.1.",
    )

    # ---- FORMING GENE GROUPS ----
    parser.add_argument(
        "--geneGroupPvalueThreshold",
        action="store",
        type=float,
        default=5e-8,
        help="(Optional) only genes sharing cis-eQTLs with P<this threshold will be grouped together for multivariable MR. The default is 5e-8, and we recommend a threshold at or very close to genome-wide significance (5e-8) because this ensures that gene groups do not become exceedingly large. If they become exceedingly large, causal estimation may be unstable and there is a greater dependence on reliable missing data imputation.",
    )
    parser.add_argument(
        "--onlyWithinMbOfPhenoLocus",
        action="store",
        type=float,
        default=1,
        help="(Optional) Only form gene groups within this many Mb of an outcome locus if -analysisInPhenotypeLoci is True. The default is 1 (1Mb).",
    )
    parser.add_argument(
        "--numIndexGenesToFormGroup",
        action="store",
        type=int,
        default=4,
        help="(Optional) If -analysisInPhenotypeLoci is True, the -numIndexGenesToFormGroup closest to each lead SNP indexing a phenotype locus will be used to form a group of genes. Note, these genes must also be within -onlyWithinMbOfPhenoLocus of the lead phenotype SNP for the locus. The default is 4.",
    )

    # ---- FORMING MR IV SETS ----
    parser.add_argument(
        "--ldMax",
        action="store",
        type=float,
        default=0.5,
        help="(Optional) two cis-eQTLs with an LD correlation greater than this threshold will be considered as the same cis-eQTL and the one with the largest eQTL GWAS P-value will be excluded",
    )
    parser.add_argument(
        "--LDOtherLoci",
        action="store",
        type=float,
        default=0.2,
        help="(Optional) IVs that are in absolute LD beyond this threshold with SNPs in other loci may be excluded (see -POtherLoci). The default is 0.2.",
    )
    parser.add_argument(
        "--POtherLoci",
        action="store",
        type=float,
        default=0.05,
        help="(Optional) IVs that are in absolute LD beyond -LDOtherLoci with SNPs in other loci that have P-values less than this threshold may be excluded (see -OtherLociWindow). The default is 0.05.",
    )
    parser.add_argument(
        "--otherLociMbWindow",
        action="store",
        type=float,
        default=0.05,
        help="(Optional) IVs that meet the conditions of -LDOtherLoci and -POtherLoci and are within this Mb window will be excluded. The default is 1, corresponding to 1Mb.",
    )
    parser.add_argument(
        "--MVMRIVPThreshold",
        action="store",
        type=float,
        default=5e-8,
        help="(Optional) IVs must have P<this threshold in a joint test for multiple genes to be considered as IVs in multivariable MR. The default is 5e-8.",
    )
    parser.add_argument(
        "--UVMRIVPThreshold",
        action="store",
        type=float,
        default=5e-5,
        help="(Optional) IVs must have P<this threshold in the eQTL GWAS to be considered as IVs in univariable MR. The default is 5e-5.",
    )

    # ---- MR CAUSAL ESTIMATION ----
    parser.add_argument(
        "--imputeMissingEQTLs",
        action="store",
        type=str,
        default="yes",
        help='(Optional) Should missing values in the MVMR IV set (which is a union set of gene-specific IV sets) be imputed using the matrix completion method described in Lorincz-Comi & Yang et al (2023)? If you put "no" or "False" here, missing values will still be imputed, but with all 0s. The default is "yes".',
    )
    parser.add_argument(
        "--nSNPsBiasCorrection",
        action="store",
        type=int,
        default=50,
        help="(Optional) The minimum number of non-significant SNPs to be used when calculating bias-correction terms used by MR-Jones, following the methods of Lorincz-Comi & Yang et al. (2023). The default is 50.",
    )
    parser.add_argument(
        "--PvalueThresholdBiasCorrection",
        action="store",
        type=float,
        default=0.01,
        help="(Optional) Only SNPs with P<this threshold will be used to calculatee the bias-correction terms in MR-Jones, following the methods of Lorincz-Comi & Yang et al. (2023). The default is 0.01.",
    )
    parser.add_argument(
        "--assumeNoSampleOverlap",
        action="store",
        type=str,
        default="yes",
        help='(Optional) Should it be assumed that the eQTL and disease GWAS do not share any of the same participants? If this assumption should not be made, bias-correction for sample overlap will be made following the methods of Lorincz-Comi & Yang et al. (2023). The default is "yes".',
    )
    parser.add_argument(
        "--minMVMRIVs",
        action="store",
        type=int,
        default=30,
        help="(Optional) The minimum number of IVs that may be used in multivariable MR with MR-Jones. If less than this threshold are available for the locus, the locus will be skipped. The default is 30.",
    )
    parser.add_argument(
        "--hessMinScale",
        action="store",
        type=float,
        default=2.5,
        help="(Optional) The signal-to-noise ratio for associations between SNPs and gene expression. Larger values will restrict the analysis to only genes with the strongest signals. Values equal to or less than 1/2 may lead to a Hessian matrix with diagonal elements that are not the correct sign, which is introduced by the MRBEE/MR-Jones correction for weak instruments. Note that the bias-correction term in the Hessian matrix of MR-Jones is automatically multiplied by the factor 1/2. The default is 2.5, corresponding to an effective signal-to-noise ratio of 5.",
    )
    parser.add_argument(
        "--geneScreener",
        action="store",
        type=str,
        default="mrjones",
        help='(Optional) The tool to use to perform gene selection in the screening step of penalized eQTL-MVMR. The options are "gscreen" for the G-Screen method and "mrjones" for the MR-Jones method. The default is "mrjones"',
    )

    # --- POST-ESTIMATION ANALYSES AND OUTPUT OPTIONS ---
    parser.add_argument(
        "--adjustSEsForInflation",
        action="store",
        type=str,
        default="yes",
        help='(Optional) One of True or False. Should the standard errors of causal effect estimates be adjusted for inflation due to misspecified LD structure? The default is "yes".',
    )
    parser.add_argument(
        "--saveRawData",
        action="store",
        type=str,
        default="false",
        help='(Optional) One of "true" "yes", "false" or "no" . Should the raw data used in multivariable MR be saved? This includes association estimates with each gene in a group, their standard errors, correspondingly the same for the outcome phenotype, and the estimated LD matrix. NOTE!!! This really only works if you are exploring a single candidate gene/locus',
    )
    parser.add_argument(
        "--whereSaveRawData",
        action="store",
        type=str,
        default="",
        help='(Optional) If you gave --saveRawData "yes" or "true", this should be the directory/folder in which you want the raw data saved. The files will have the names "bybx.txt" (where -by- is in the first column),"UU.txt", and "LD.txt"',
    )
    parser.add_argument(
        "--out",
        action="store",
        type=str,
        default="results",
        help="(Optional) filename without extension of location in which results should be written, in tab-separated .txt format. The default is results. Note that a tab-separated file named diagnostics.txt will automatically be written. This file contains information about the performance of multivariable MR.",
    )
    parser.add_argument(
        "--iterativelySave",
        action="store",
        type=str,
        default="true",
        help="(Optional) Should results be saved iteratively as each chromosome is completed? This is helpful if you anticipate the analysis may take a relatively long time and you do not want to lose progress in case your access to the machine it is running on expires. The default is True.",
    )
    parser.add_argument(
        "--cleanResults",
        action="store",
        type=str,
        default="yes",
        help='(Optional) Should the minimum interested results be saved? Put "no" or "false" if you want the complete results, including causal estiamtes from other MR methods. The default is "yes", to onnly output the most important results',
    )

    # ---- NETWORK OPTIONS ----
    parser.add_argument(
        "--graphNetworks",
        action="store",
        type=str,
        default="yes",
        help='(Optional) Should graph images of causal networks be stored? If "yes", they will go into the "plots" subfolder. The default is "yes".',
    )
    parser.add_argument(
        "--graphNetworksFilePrefix",
        action="store",
        type=str,
        default="",
        help='(Optional) The default prefix for the filenames of graphed regulatory networks. Note that all network graph files will contain the lead gene in the filename automatically. What you put here will come before that. The default is "".',
    )
    parser.add_argument(
        "--networkR2Thres",
        action="store",
        type=float,
        default=0.25,
        help="(Optional) Only loci with at least this much genetic variance in the disease phenotype explained will be considered for graph analysis using MNet. The default is 0.25.",
    )
    parser.add_argument(
        "--networkDiseaseLabel",
        action="store",
        type=str,
        default="disease",
        help='(Optional) In graphs of causal networks, the disease node will have this label. This default is "disease"',
    )
    parser.add_argument(
        "--networkGraphsOut",
        action="store",
        type=str,
        default="plots",
        help='(Optional) The directory/folder in which you would like plots of causal network graphs to be saved. The default is the "plots/" subdirectory/subfolder of the "hornet/" subdirectory/subfolder.',
    )

    # --- OUTPUT PRINTING OPTIONS ---
    parser.add_argument(
        "--silence",
        action="store",
        type=str,
        default="no",
        help='(Optional) Should warnings about the size of the CHP window outside of the target locus be ignored? Put "true" or "yes". The default is "no".',
    )
    parser.add_argument(
        "--hideVersion",
        action="store",
        type=str,
        default="yes",
        help='(Optional) Should the version number of HORNET be printed at the beginning? If not, put "false" or "no". The default is "yes".',
    )

    args = parser.parse_args()

    return args
