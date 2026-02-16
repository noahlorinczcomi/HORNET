import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        prog="",
        description="This program will search your results.txt file to find novel loci; \n i.e., loci that were not detected at a standard significance threshold in the original disease GWAS but were detected by HORNET. \n Negative confounding between genes is one explanation for why this may happen.",
    )

    # --- DATA ---
    parser.add_argument(
        "-pg",
        "--phenoGWAS",
        action="store",
        type=str,
        help="(Required) path to your phenotype GWAS file of summary statistics (with file extension).",
    )
    parser.add_argument(
        "-ld",
        "--LDRef",
        action="store",
        type=str,
        default=os.path.abspath("data/ldref/1kg3TRANS"),
        help='(Optional) path to the LD reference panel file (without extension). the default is"data/ldref/1kg3TRANS".',
    )
    parser.add_argument(
        "-r",
        "--results",
        action="store",
        type=str,
        default="results.txt",
        help="(Required) filepath to your results file (with file extension).",
    )

    # --- DATA METADATA ----
    parser.add_argument(
        "-rsid",
        "--phenoRSID",
        action="store",
        type=str,
        default="SNP",
        help="(Required) column name in phenotype GWAS data indicating rsID.",
    )
    parser.add_argument(
        "-pz",
        "--phenoZ",
        action="store",
        type=str,
        default="Z",
        help="(Required if no P-value column in phenotype GWAS data) column name in phenotype GWAS indicating Z-statistics for associations between SNPs and phenotype.",
    )
    parser.add_argument(
        "-pp",
        "--phenoP",
        action="store",
        type=str,
        default="missing",
        help="(Required if no Z-statistic column in phenotype GWAS data) column name indicating P-value for association in the phenotype GWAS data.",
    )

    # --- HYPERPARAMETERS ---
    parser.add_argument(
        "-psig",
        "--phenoPSignifThres",
        action="store",
        type=float,
        default=5e-8,
        help='(Optional) P-value threshold for classifying an SNP as "significantly" associated with the outcome/phenotype. The default is 5e-8.',
    )
    parser.add_argument(
        "-pr2",
        "--phenoR2Thres",
        action="store",
        type=float,
        default=0.1,
        help="(Optional) SNPs in LD beyond this threshold with a lead outcome/phenotype SNP will be considered part of the same locus as the lead outcome SNP. The default is 0.1.",
    )
    parser.add_argument(
        "-pkb",
        "--phenoKbThres",
        action="store",
        type=int,
        default=500,
        help="(Optional) Size of window in which PLINK should clump SNPs to identify idependent outcome signals. The default is 500, which means 500Kb.",
    )
    parser.add_argument(
        "-hpv",
        "--hornetSignifPThres",
        action="store",
        type=float,
        default=5e-5,
        help="(Optional) Genes with a P-value from HORNET less than this threshold will be considered detected/statistically significant. The default is 5e-5.",
    )
    parser.add_argument(
        "-hr",
        "--hornetRSquaredThres",
        action="store",
        type=float,
        default=0.10,
        help="(Optional) Genes in a locus that explains at least this much outcome/phenotype variance will be considered detected/statistically significant. The default is 0.1.",
    )
    parser.add_argument(
        "-hpr",
        "--hornetPrattThres",
        action="store",
        type=float,
        default=0.05,
        help="(Optional) Genes with a Pratt index value (gene-specific outcome variance attributable) from HORNET greater than this threshold will be considered detected/statistically significant. The default is 0.05.",
    )
    parser.add_argument(
        "-n",
        "--novelKbWindow",
        action="store",
        type=float,
        default=500,
        help="(Optional) Genes detected in HORNET must be at least this many Kb away from any known outcome signals to be considered novel. The default is 500, which means 500Kb.",
    )

    # --- OUTPUT ---
    parser.add_argument(
        "-m",
        "--message",
        action="store",
        type=str,
        default="yes",
        help='(Optional) should the warning message be printed? Default is "yes".',
    )
    parser.add_argument(
        "-p",
        "--print",
        action="store",
        type=str,
        default="yes",
        help='(Optional) Should results also be printed to the console? The default is "yes".',
    )
    parser.add_argument(
        "-o",
        "--out",
        action="store",
        type=str,
        default="newgenes",
        help='(Optional) Name of file to which results should be written (without file extension; it will be .csv). The default is "newgenes".',
    )

    args = parser.parse_args()

    return args
