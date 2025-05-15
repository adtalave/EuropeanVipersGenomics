# VCF Filtering
Once we have successfully called our genotypes we can start filtering out not trustworthy SNPs.
## GATK Hard filtering
First, we followed the proposed filtering pipeline by GATK focused on six annotation parameters that describe each position, i.e., QD, FS, SOR, MQ, MQRankSum and ReadPosRankSum. 
We can generate a table with these parameters with GATK and subsequently produce density plots of each one in R. 

```bash
#!/bin/bash
gatk="path/to/gatk"
chromosomes="list/of/chromosomes.txt"

for i in $(cat $chromosomes)
do
(
$gatk VariantsToTable -V "$i"_dataset.vcf.gz -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum -O "$i"_.table
) & 
done
wait
cat *_.table > allCHR.table
```
Different thresholds based on the distribution of your parameters should be applied to different projects. Read more about GATK's Hard filtering [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).
To filter youre genotype VCF with the selected thresholds:
```bash
#!/bin/bash
gatk="path/to/gatk"
out="path/to/output/directory"
reference="path/to/your/reference.fasta"

$gatk SelectVariants -R $reference \
-V "$i"_dataset.vcf.gz \
-select "QD > 5.0 && FS < 50.0 && SOR < 4.0 && MQ > 55.0 && MQRankSum > -1.5 && MQRankSum < 1.5 && ReadPosRankSum > -3.0 && ReadPosRankSum < 3.0" \
-O $out/hardfiltered_dataset.vcf.gz
```
## Normalizing indels
```bash
#!/bin/bash
reference="path/to/your/reference.fasta"

bcftools norm -f $reference -o normalizedfiltered_dataset.vcf.gz hardfiltered_dataset.vcf.gz
```


## Removing 10 SNPs around indels
```bash
#!/bin/bash
reference="path/to/your/reference.fasta"
bcftools filter --SnpGap 10 -o 10SNPfiltered_dataset.vcf.gz normalizedfiltered_dataset.vcf.gz
```

## Removing SNPs in repetitive regions
```bash
#!/bin/bash
reference="path/to/your/reference.fasta"
masked="path/to/your/sofmasked/reference.repeats.bed" #you can get this file with RED for example

bedtools intersect -header -v -a 10SNPfiltered_dataset.vcf.gz -b "$masked" | bgzip -c > masked_dataset.vcf.gz
```
Additionally, you can apply other masks such as mappability, explained [here](http://lh3lh3.users.sourceforge.net/snpable.shtml). 
## Other filters
Adjust the minimum allele frequency to the number of samples of your VCF, taking into account whether you'd like to keep singletons, for example.
```bash
#!/bin/bash
vcftools --gzvcf masked_dataset.vcf.gz --min-alleles 2 --max-alleles 2 --remove-indels --maf 0.006 --minQ 30 --min-meanDP 4 --max-meanDP 150 --minDP 4 --maxDP 150 --max-missing 0.9 --recode --stdout | gzip -c > filtered_dataset.vcf.gz
```
## LD-pruning 
For population genomics analyses you can calculate how Linkage Disequilibrium decays and consequently prune your data to keep only unlinked SNPs.
To calculate LD decay within a "population":
```bash
PopLDdecay/bin/PopLDdecay -InVCF  1pop_autosomes.vcf.gz  -MaxDist 40 -OutStat 1pop.stat
```
Then you can thin your dataset accordingly:
```bash
vcftools --gzvcf filtered_dataset.vcf.gz --thin 500 --recode --stdout | gzip -c > unlinked_dataset.vcf.gz 
```
