This tutorial demonstrates how to download and use the HORNET software to perform genome-wide searches for genes whose expression may cause disease risk. HORNET performs multivariable Mendelian Randomization (MR) using the MR with JOint selectioN of exposurES and pleiotropy (MR-Jones) method. MR-Jones is an extension of the MR with Bias-corrected Estimating Equations (MRBEE) method to the high dimensional setting and adjusts for horizontal pleiotropy using penalized regression. The HORNET, MR-Jones, and MRBEE papers can be found from the following references:
1. (**HORNET**) HORNET preprint
2. (**MR-Jones**) preprint
3. Lorincz-Comi, N., Yang, Y., Li, G., & Zhu, X. (2023). MRBEE: A novel bias-corrected multivariable Mendelian Randomization method. *bioRxiv*, 2023-01. DOI: https://doi.org/10.1101/2023.01.10.523480

# Downloading HORNET
This tutorial will use the Linux terminal to perform all operations and Python v3.7 must be executable from the command line. Note that later versions of Python may become incompatible with HORNET as some of its dependencies are updated and certain functionality may become deprecated. We therefore recommend using Python v3.7. 

Start by cloning the HORNET Github repository like this:

```unix
git clone https://github.com/noahlorinczcomi/HORNET.git
cd HORNET
```

Next, you must copy the `hornet_data` tarball from our public web page using the following command:
```unix
wget -O hornet_data.tar http://hal.case.edu/~njl96/hornet_data.tar.gz
tar -xvf hornet_data.tar
rm hornet_data.tar
```

Now, in the `HORNET` directory, you should see the following folders and files:
1. `data/`
    * `ldref/` # LD reference panels
        * `1kg_phase3_<EUR/AFR/EAS/SAS/HIS/TRANS>_only.<bed/bim/fam>`
    * `maps/` # mapfiles for converting from hg19 to hg38 and vice versa
        * `1kgPhase3maps/`
            * `1kgPhase3hg19Chrom<CHR>MapToHg38.txt.gz`
        * `fullMaps/`
            * `1kgPhase3MapAllSNPs.txt`
        * `largeMaps/`
            * `hg19Chrom<CHR>MapToHg38.txt.gz`
2. `tempfiles/` # temporary files will be iteratively written here when running HORNET
3. `testdata/`  # test/example data is included here for demonstrative purposes
    * `GTEx`
    * `eQTLGen`
4. `plinkdir/`  # the PLINK v1.9 software (Purcell etal) is here
5. `LD_block_finder.r` # an R program to find blocks in an LD matrix using the method in Lorincz-Comi et al. ***
6. `hornet.py`

The HORNET software requires a number of Python modules which may or may not be available to you already. 

# Downloading summary GWAS data
For this tutorial, we have already downloaded GWAS summary data for Alzheimer's disease (AD) from Jansen et al. (2019) and gene expression data in blood from the eQTLGen Consortium (--) and in frontal cortical tissue from the GTEx Consortium (-) for chromsomes 1 and 2. 

Gene expression GWAS data are stored in the `testdata/` directory under the names `eQTLGen/eQTLGen_blood_CHR1.txt.gz` and `GTEx/GTEx_frontal_cortex_CHR1.txt.gz`, where data sets from these two sources are placed into separate folders for reasons that will become apparent later.

AD GWAS data are stored in the `testdata/` directory under the name `AD_Jansen_etal.txt.gz`.

# Using HORNET
HORNET is a command-line tool. This means that the user defines a single command and executes it in the terminal. The command that the user gives HORNET is slightly different if they are using raw summary data from GTEx vs some other source, such as eQTLGen. First, we demonstrate how to use HORNET with eQTL GWAS data that is **not** from GTEx.

## eQTL summary data NOT from GTEx
If the eQTL GWAS summary data is not from GTEx, HORNET must be told the following information:
* Filepath to the folder where the eQTL GWAS summary data are located (nothing else besides these data must be present in this folder)
* Filepath to the phenotype/disease GWAS summary data
* Names of columns representing rsID, Z-score OR effect size and standard error, and effect alleles in both eQTL and phenotype GWAS data sets
* Names of columns representing SNP base pair position, gene label, and gene base pair position in the eQTL GWAS data set

An example command using the eQTLGen GWAS data and AD is this:

```unix
python hornet.py \ 
 --eQTLGWAS testdata/eQTLGen \
 --phenoGWAS testdata/AD_Jansen_etal.txt.gz \
 --LDRef data/ldref/1kg3 \
 --isRawGTEx False \
 --snpLabels SNP,SNP \
 --eQTLSNPBPLabel SNPPos \
 --zLabels Zscore,Z \
 --effectAlleles AssessedAllele,A1 \
 --eQTLGeneLabel Gene \
 --eQTLGeneBP GenePos \
```

where the second value in each comma-separated argument corresponds to the phenotype (AD) and the first to the eQTL GWAS data.

## eQTL summary data from GTEx
Certain aspects of the command given to HORNET change when the eQTL GWAS data is directly from GTEx. Here is the basic command to give HORNET in this case:

```unix
python hornet.py \ 
 --eQTLGWAS testdata/GTEx \
 --phenoGWAS testdata/AD_Jansen_etal.txt.gz \
 --LDRef data/ldref/1kg3EUR \
 --isRawGTEx True \
 --snpLabels gtex,SNP \
 --eQTLSNPBPLabel gtex \
 --zLabels gtex,Z \
 --effectAlleles gtex,A1 \
 --eQTLGeneLabel gtex \
 --eQTLGeneBP gtex \
```

Here, all arguments indicating column names that correspond to the eQTL GWAS data set were replaced with `gtex`. The key argument is switching `-isRawGTEx` from False to True. This tells HORNET all it needs to load the data. Raw GTEx summary data do not contain rsIDs, which HORNET will automatically add using mapfiles in the `data/maps/1kgPhase3maps/`. directory. 

There are additional arguments that HORNET accepts, which relate to how gene groups are formed, IVs are selected, LD is handled, and the causal estimates are made. These option can be inspected more closely by looking at the output of:

```unix
./hornet.py --help
```

