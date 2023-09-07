echo 'Starting Cortex & AD analysis'
python hornet.py ^
 --eQTLGWAS data/testdata/MetaBrain ^
 --phenoGWAS data/testdata/AD_Jansen_etal.txt.gz ^
 --LDRef data/ldref/1kg3EUR ^
 --isRawGTEx False ^
 --snpLabels SNP,SNP ^
 --eQTLSNPBPLabel SNPPos ^
 --zLabels MetaBeta@MetaSE,Z ^
 --effectAlleles SNPEffectAllele,A1 ^
 --eQTLGeneLabel Gene ^
 --eQTLGeneBP GenePos ^
 --MVMRIVPThreshold 5e-3 ^
 --networkR2Thres 0.1 ^
 --geneGroupPvalueThreshold 1e-3 ^
 --out results/AD_cortex

echo 'Starting Thyroid & AD analysis'
 python hornet.py ^
 --eQTLGWAS data/testdata/GTEx ^
 --phenoGWAS data/testdata/AD_Jansen_etal.txt.gz ^
 --LDRef data/ldref/1kg3EUR ^
 --isRawGTEx yes ^
 --snpLabels gtex,SNP ^
 --eQTLSNPBPLabel gtex ^
 --zLabels gtex,Z ^
 --effectAlleles gtex,A1 ^
 --eQTLGeneLabel gtex ^
 --eQTLGeneBP gtex ^
 --MVMRIVPThreshold 5e-3 ^
 --networkR2Thres 0.1 ^
 --geneGroupPvalueThreshold 1e-3 ^
 --out results/AD_thyroid
