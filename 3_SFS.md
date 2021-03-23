# Rice blast variant calling
# 3. Effect of SNP filtered sets on the Site Frequency Spectra


Program                    | Link
-------------------------- | ---------------------------------
*bcftools 1.11*            | (http://www.htslib.org/download/)
*calculate_proportions.sh* | (https://github.com/Burbano-Lab/rice-blast-variant-calling/blob/main/scripts/calculate_proportions.sh)


After the classification of the isolates in different genetic groups, we divided the SNPs segregating in each genetic group.

```bash
for group in group{1..4}.list; do
	bcftools view -S $group -a -f PASS ALL.snps.vcf.gz | bcftools view -m2 -M2 -g ^miss - | bgzip > ${group%.list}.snps.filtered.vcf.gz
done
```
We finally counted the within-genetic-group SNP type frequency, namely fixed, singleton or non-singleton (non-fixed SNPs present in more than one isolate).
```bash
for group in group{1..4}.snps.filtered.vcf.gz; do
	bash calculate_proportions.sh $group > ${group%.vcf.gz}.proportions.txt
done
```

