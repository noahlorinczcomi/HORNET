import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        prog="This tool can help guide you in choosing which tissues are good choices for performing genome-wide eQTL-MVMR with MR-Jones",
        description="This program relies on pre-calculated heritability scores, which are proportional to SNP heritability, from cis-SNPs and GTEx v8 summary data to create an ordered list of tissues in which your candidate genes have the strongest eQTL signals.",
    )

    parser.add_argument(
        "--candidateGenes",
        action="store",
        type=str,
        help='(Required) This can either be a comma-separated list of gene IDs or the path to a file (with extension) that contains one gene ID per row and no header. NOTE, all gene IDs must follow the Ensembl gene ID format (e.g., "ENSG..."), but they do not need the additional ".<xx>" identifier. For example, "ENSG00000227232" is sufficient, you do not need to use "ENSG00000227232.5"',
    )
    parser.add_argument(
        "--nTopTissues",
        action="store",
        type=int,
        default=10,
        help='(Required) The number of top tissues to display. The default is 10. Note, if you put "yes" for --saveAsFile, all results will be written out, not just these that are printed.',
    )
    parser.add_argument(
        "--saveAsFile",
        action="store",
        type=str,
        default="no",
        help='(Optional) If you would like the results to be saved to a file, put "true" or "yes", otherwise ignore this flag or put "false" or "no"',
    )
    parser.add_argument(
        "--outFilename",
        action="store",
        type=str,
        default="tissue_results.txt",
        help='(Optional) If you put "true" or "yes" for the --saveAsFile flag, you can specify the filename (without extension) of the output file here. The default is "tissue_results.txt"',
    )
    parser.add_argument(
        "--printResults",
        action="store",
        type=str,
        default="yes",
        help='(Optional) Please give this flag "no" if you do not want the results to be printed to the console. Note that if you give this a "no", you will want to make sure the results are saved in a file, otherwise you will never see the results. The default is to print results to the console.',
    )

    args = parser.parse_args()
    return args
