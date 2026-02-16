import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        prog="Summarise HORNET Results",
        description='This program subsets the full results to only significant genes, where you defined "significant"',
    )

    # ---- DATA ----
    parser.add_argument(
        "-r",
        "--results",
        action="store",
        type=str,
        help="(Required) Full filepath (w/ extension) to results outputted by hornet.py",
    )

    # ---- HYPERPARAMETERS ----
    parser.add_argument(
        "-pr",
        "--prattThreshold",
        action="store",
        type=float,
        default=0.1,
        help="(Optional) Lower threshold for Pratt index values. The default is 0.1.",
    )
    parser.add_argument(
        "-v",
        "--localVarianceExplainedThreshold",
        action="store",
        type=float,
        default=0.5,
        help="(Optional) The locus in which a prioritized gene sits must have explained at least this much local genetic variance in the disease phenotype. The default is 0.5.",
    )
    parser.add_argument(
        "-pv",
        "--pValueThreshold",
        action="store",
        type=float,
        default=5e-5,
        help="(Optional) Only genes with P<this threshold will be considered as priorities. The default is 5e-5.",
    )
    parser.add_argument(
        "-tk",
        "--topKGenes",
        action="store",
        type=int,
        default=20,
        help="(Optional) Only the top --topKGenes will be printed (sorted by Pratt index). The default is 20.",
    )

    # ---- OUTPUT ----
    parser.add_argument(
        "--print",
        action="store",
        type=str,
        default="yes",
        help="(Optional) Should the results be printed in addition to being saved? The default is yes.",
    )
    parser.add_argument(
        "-o",
        "--out",
        action="store",
        type=str,
        default="",
        help='(Optional) Filename (without extension, with possible prefix directory) to which results should be written. The default is "", which means no results are saved.',
    )

    args = parser.parse_args()
    return args
