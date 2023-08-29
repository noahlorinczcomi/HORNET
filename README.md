This tutorial demonstrates how to download and use the HORNET software to perform genome-wide searches for genes whose expression may cause disease risk. 

HORNET performs multivariable Mendelian Randomization (MR) using the MR with JOint selectioN of exposurES and pleiotropy (MR-Jones) method. MR-Jones is an extension of the MR with Bias-corrected Estimating Equations (MRBEE) method to the high dimensional setting and adjusts for horizontal pleiotropy using penalized regression.

The HORNET, MR-Jones, and MRBEE papers can be found from the following references:
1. (**HORNET**) HORNET preprint
2. (**MR-Jones**) preprint
3. Lorincz-Comi, N., Yang, Y., Li, G., & Zhu, X. (2023). MRBEE: A novel bias-corrected multivariable Mendelian Randomization method. *bioRxiv*, 2023-01. DOI: https://doi.org/10.1101/2023.01.10.523480

# Downloading HORNET
All methods of downloading HORNET require that you are able to type 'python' at your command line and a Python (v3.7.0 or later) console opens. This can be easily achieved by downloading [Anaconda](https://www.anaconda.com/download/).

## Windows users
Download the `windowsconfig.bat` file from [here](https://github.com/noahlorinczcomi/HORNET/blob/main/windowsconfig.bat) (there is a download link on the page) and move it to some folder where you would like HORNET to be. Double-click on the file and HORNET will be installed. You can run HORNET from the command line by copy/pasting the commands in the tutorial below, or you can use the desktop GUI version.

## Linux/Mac users
If you are using Linux or a Mac, copy/paste the commands below into your terminal or Anaconda prompt.

```unix
git clone https://github.com/noahlorinczcomi/HORNET.git && cd HORNET
wget -O hornet_data.tar http://hal.case.edu/~njl96/hornet_data.tar.gz
tar -xvf hornet_data.tar && rm hornet_data.tar
mkdir tempfiles plots results
pip install -r requirements.txt
```

Lastly, depending if you're using Linux or Mac, execute on of the following:
```unix
chmod +x ./plinkdir/linux/plink
chmod +x ./plinkdir/mac/plink
```

The creation of `tempfiles/` is to help you. The program accepts a flag `--writableDir` which expects `tempfiles/` to be given to it. If `tempfiles/` does not exist, you must specify the `--writabledir` flag with any folder/directory that can be iteratively written into while HORNET runs. The creation of `plots/` and `results/` is also to help you. HORNET will try saving plots and results to `plots/` and `results/` folder by default. You can change these by giving something else to the `--networkGraphsOut` and `--out` flags.

<!---
Now, in the `HORNET` directory, you should see the following folders and files:
1. `data/`
    * `ldref/` # *LD reference panels*
        * `1kg_phase3_<EUR/AFR/EAS/SAS/HIS/TRANS>_only.<bed/bim/fam>`
    * `maps/` # *folder containing mapfiles for converting from hg19 to hg38 and vice versa*
        * `1kgPhase3maps/`
            * `1kgPhase3hg19Chrom<CHR>MapToHg38.txt.gz` # *CHR-specific hg19/hg38 maps of only 1kg Phase 3 SNPs*
        * `fullMaps/`
            * `1kgPhase3MapAllSNPs.txt` # *hg19/hg38 map using all 1kg Phase 3 SNPs*
        * `largeMaps/`
            * `hg19Chrom<CHR>MapToHg38.txt.gz` # *large hg19/hg38 map from GLGC (Graham et al. 2021) including 44M SNPs*
        * `testdata/` # *test/example data is included here for demonstrative purposes*
            * `GTEx/`
                * `GTEx_frontal_cortex_CHR1.txt.gz` # *all associations b/w SNPs and CHR 1 gene expression in frontal cortex (GTEx v8)*
                * `GTEx_frontal_cortex_CHR2.txt.gz` # *all associations b/w SNPs and CHR 2 gene expression in frontal cortex (GTEx v8)*
            * `eQTLGen/`
                * `eQTLGen_blood_CHR1.txt.gz` # *all associations b/w SNPs and CHR 1 gene expression in whole blood (eQTLGen)*
                * `eQTLGen_blood_CHR2.txt.gz` # *all associations b/w SNPs and CHR 2 gene expression in whole blood (eQTLGen)*
2. `tempfiles/` # *temporary files will be iteratively written here when running HORNET*
3. `plinkdir/`  # *the PLINK v1.9 software (Purcell etal) is here*
4. `LD_block_finder.r` # *an R program to find blocks in an LD matrix using the method in Lorincz-Comi et al. (***) *
5. `hornet.py`    # *executable HORNET program*
6. `functions.py` # *source file of functions that HORNET uses*
--->

# Tutorial
This tutorial will use the command line to perform analyses with HORNET. We still need access to python from our command prompt, which can be achieved if we are using the Anaconda prompt.

## Background on what we are going
We want to perform genome-wide searches for genes with evidence of causality with a phenotype when they are expressed in a particular tissue. (You can test multiple tissues separately using HORNET multiple times). This means that need eQTL GWAS summary statistics that contain estimates of association between SNPs and the expression of all genes within some Mb window of them (usually +-1Mb). Since the scale of these data can be enormous when combined across the entire measurable genome, one file for each chromosome should exist in a folder. Nothing else should exist in this folder.

For this tutorial, we have already downloaded GWAS summary data for Alzheimer's disease (AD) from [Jansen et al. (2019)](https://doi.org/10.1038/s41588-018-0311-9) and gene expression data in blood from the [eQTLGen Consortium](https://www.eqtlgen.org/) and in frontal cortical tissue from the [GTEx Consortium](https://gtexportal.org/home/) for chromsomes 1 and 2. *Note*, you can remove these data sets from the `HORNET/` directory to save space after completing this tutorial.

Gene expression GWAS data are stored in the `testdata/` directory under the names `eQTLGen/eQTLGen_blood_CHR<1,2>.txt.gz` and `GTEx/GTEx_frontal_cortex_CHR<1,2>.txt.gz`, where data sets from these two sources are placed into separate folders for reasons that will become apparent later.

AD GWAS data are stored directly in the `testdata/` directory under the name `AD_Jansen_etal.txt.gz`.

![hi](https://github.com/noahlorinczcomi/HORNET/blob/main/HORNET%20Computational%20Flowchart.svg)

## eQTL summary data NOT from GTEx
If the eQTL GWAS summary data is not from GTEx, HORNET must be told the following information:
* Filepath to the folder where the eQTL GWAS summary data are located (nothing else besides these data must be present in this folder)
* Filepath to the phenotype/disease GWAS summary data
* Names of columns representing rsID, Z-score OR effect size and standard error, and effect alleles in both eQTL and phenotype GWAS data sets
* Names of columns representing SNP base pair position, gene label, and gene base pair position in the eQTL GWAS data set

An example command (in UNIX) using the eQTLGen GWAS data and AD is this:

```unix
python hornet.py --eQTLGWAS data/testdata/eQTLGen \
 --phenoGWAS data/testdata/AD_Jansen_etal.txt.gz \
 --LDRef data/ldref/1kg3EUR \
 --isRawGTEx False \
 --snpLabels SNP,SNP \
 --eQTLSNPBPLabel SNPPos \
 --zLabels Zscore,Z \
 --effectAlleles AssessedAllele,A1 \
 --eQTLGeneLabel Gene \
 --eQTLGeneBP GenePos \
 --out results/AD_blood
```

where the second value in each comma-separated argument corresponds to the phenotype (AD) and the first to the eQTL GWAS data. See all files in the `data/ldref/` directory if you want to choose a different 1kg Phase 3 population as the LD reference panel.

**Be sure to replace '\\' with '^' if you're on a Windows machine.**

## eQTL summary data from GTEx
Certain aspects of the command given to HORNET change when the eQTL GWAS data is directly from GTEx. This is because HORNET knows exactly what to expect from eQTL GWAS data from GTEx. Here is the basic command to give HORNET in this case:

```unix
python hornet.py --eQTLGWAS data/testdata/GTEx \
 --phenoGWAS data/testdata/AD_Jansen_etal.txt.gz \
 --LDRef data/ldref/1kg3EUR \
 --isRawGTEx yes \
 --snpLabels gtex,SNP \
 --eQTLSNPBPLabel gtex \
 --zLabels gtex,Z \
 --effectAlleles gtex,A1 \
 --eQTLGeneLabel gtex \
 --eQTLGeneBP gtex \
 --out results/AD_frontal_cortex
```

Here, all arguments indicating column names that correspond to the eQTL GWAS data set were replaced with `gtex`, but the only key change is switching `-isRawGTEx` from 'no' to 'yes'. This tells HORNET all it needs to load the data. Raw GTEx summary data do not contain rsIDs, which HORNET will automatically add using mapfiles in the `data/maps/1kgPhase3maps/` directory. 

## Additional flags
There are additional arguments that HORNET accepts, which relate to how gene groups are formed, IVs are selected, LD is handled, and the causal estimates are made. These options can be inspected more closely by looking at the output of:

```unix
python hornet.py --help
```

# HORNET Results
## Causal estimates and network graphs
**Assume we ran the above command for GTEx data**. 

We gave the `--out` flag the value `results/AD_frontal_cortex`, which is the prefix of a file name in the `results/` folder. This will give us the following results:
1. `HORNET/results/AD_frontal_cortex_results.txt`
    * Causal estimates and model fit
    * Other filepaths can be specified by changing the argument given to the `--out` flag
2. `HORNET/results/AD_frontal_cortex_diagnostics.txt`
    * Summary information related to missingness and imputation, the initial size of the gene network, the size of the IV set after applying various QC, etc.
3. `HORNET/plots/<lead gene ID>_graph.png`
    * Plots* of gene regulatory networks that include the disease outcome, too
        - \* : Only networks for which the genetic variance in the disease that was explained by gene expression exceeded `--networkR2Thres` will be considered

### Summarised results
**Note:** The `AD_frontal_cortex_results.txt` file in this example will, by default, contain only necessary information. It is possible to save additional information, such as causal estimates from other multivariable and univariable MR methods, by specifying `--cleanResults no`. Whether you choose to do this, you can still view the top 20 genes (20 by default) using the `sumres.py` file like this:
```unix
python sumres.py --results results/AD_frontal_cortex_results.txt
```
View `python sumres.py -h` for additional paramaters you can set.

## Manhattan and volcano plots
If you're lucky enough to also have R installed and on your PATH (ie if you can type 'R' at your terminal and open an R console), you can also create Manhattan and volcano plots using the `plotres.r` file. Using our example above, you would type this:
```unix
Rscript plotres.r
```
These plots will be saved in the `HORNET/plots/` folder.
## Novel loci
It is possible that some genes detected using HORNET were not detected in the original disease/phenotype GWAS. We provide the `newloci.py` script to find these genes. In the `HORNET/` direectory, you can type `python newloci.py -h` to see how to use it, but here is the basic usage with the example from above:
```unix
python newloci.py --results results/AD_frontal_cortex_results.txt \
 --phenoGWAS data/testdata/AD_Jansen_etal.txt.gz \
 --phenoRSID SNP \
 --phenoP P \
 --novelKbWindow 1000
```
which will automatically print the top results to the console and will save the full results in `HORNET/newgenes.csv`.
 
 

