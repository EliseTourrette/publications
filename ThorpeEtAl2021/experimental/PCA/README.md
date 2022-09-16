The PCA analysis is done on CDS (biallelic SNPs) of the core genome (with an additional setp of LD pruning for the preprocessing of the data), aligned on 26695 reference genome

REQUIREMENTS
R v4.1.3 (use the library vcfR v1.13.0)  
PLINK v1.90b6.21  

preprocessing steps done in R
> preprocess.r

then LD pruning step
```
plink --vcf coreCDS.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out coreCDS_LD
```

and finally the PCA
```
plink --vcf coreCDS.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract coreCDS_LD.prune.in --make-bed --pca --out coreCDS_LD_PCA
```

the PCA results are then plotted using R (library tidyverse needed)
> plot.r