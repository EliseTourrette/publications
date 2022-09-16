Use the software GRAPES to calculate the rate of non-adaptive evolution

in this pipeline
do not use only the vcf file for the divergent dataset 
as it will lead to an underestimation of the number of divergent substitution
hence use the fact that I have the core genome and the whole genome of the Hp ref sequence
to align the core genome over the whole genome
and then get the positions of the core genome (compared to the whole genome)
use these positions to check which CDS are actually in the core genome
and use these core CDS for the divergent dataset (calculate dN and dS)
and to calculate LN and LS

DATA
whole ref seq
core ref seq
aligned outgroup to the ref seq (whole seq)
ref gff annotation file

REQUIREMENTS
GRAPES v1.1.1
python v3.10.0
R v4.1.1 (ggplot2 v3.3.5)


the first step is to use mummer (nucmer) to align the core (reference) genome over the whole reference genome in order to get the coordinates of the sites that are part of the core genome. <br />
the perl script is an in-house perl script (Chao) to format mummer output into a fasta file
```
nucmer --mum -p 26695_core 26695.fasta 26695_core.fasta
perl /home/public/PublicSoftware/aln_extract.pl 26695_core.delta > 26695_core_aligned.aln
```

the second step is to align the outgroup to the reference sequence (already done for another analysis)

the third step is to get the coordinates of the core CDS that are without missing values
> preprocess_coordinateCDS.py

then create a VCF file with the outgroup added to it (use the outgroup sequence aligned on the whole reference sequence)
> preprocess_outgroupVCF.py
> preprocess_outgroupVCR.r

then run these commands to have a standard vcf file (correct header)
```
head -n 3 vcf/Hp_core.vcf > header.txt
cat header.txt coreOutgroup.vcf > outgroup.vcf
cat header.txt HpOutgroup.vcf > Hp_core.vcf
rm header.txt coreOutgroup.vcf HpOutgroup.vcf
mv *.vcf vcf_outgroup
```

the fourth step is to calculate the number of non-synonymous and synonymous sites
> preprocess_nbSites.py

additionally, need the popylmorphism and divergence data
> preprocess_nbSitesPolDivergence.r


Lastly run GRAPES in the terminal
 if there is divergence data
 ```
 for file in *; do
   echo $file
   grapes -in $file -out "$file"_out -model GammaZero,GammaExpo,ScaledBeta >> "$file"_basic.out &
done
 ```

if there is divergence data
 ```
 for file in *; do
   echo $file
   grapes -in $file -out "$file"_out -model GammaZero,GammaExpo,ScaledBeta -no_div_data >> "$file"_basic.out &
done
 ```
 
 and finally, analysis of the results
 > analysis.r
 