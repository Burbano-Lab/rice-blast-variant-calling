# Rice blast variant calling
# 1. Generation of the Gold standard variants dataset and set of filters

Program                              | Location
------------------------------------ | --------------------------------------------------
*AdapterRemoval v2*                  | (https://github.com/mikkelschubert/adapterremoval)
*BWA v.0.7.12*                       | (https://github.com/lh3/bwa)
*minimap2 v.2.11*                    | (https://github.com/lh3/minimap2)
*bcftools v.1.11*                    | (http://www.htslib.org/download/)
*GATK v.3.8*                         | (https://github.com/broadinstitute/gatk)
*vcftools v.0.1.16*                  | (https://github.com/vcftools/vcftools)
*get_genotypes_from_VCF.py*          | (https://github.com/Burbano-Lab/rice-blast-variant-calling/tree/main/scripts/)
*compare_short_vs_long_genotypes.py* | (https://github.com/Burbano-Lab/rice-blast-variant-calling/tree/main/scripts/)


## Alignment of short reads to reference genome

Raw .fastq sequences were trimmed with *AdapterRemoval2*.
```bash
AdapterRemoval --file1 $sample1.R1.fastq.gz --file2 $sample1.R2.fastq.gz --gzip --basename $sample.trimmed
```

We used the Pac-bio assembly Guy11 (https://doi.org/10.1016/j.molp.2017.08.008) as the reference genome and used *BWA* to create the index.
```bash
bwa index guy11.fa
```

*BWA mem* was used to map the trimmed reads to the Guy11 indexed reference genome.
```bash
bwa mem guy11.fa $sample1.trimmed.R1.fastq.gz $sample1.trimmed.R1.fastq.gz
```

## Mapping Nanopore alignments to reference genome

We used *minimap2* to align the Nanopore alignments into the same coordinate system.
```bash
minimap2 -x asm5 -a guy11.fa $sample1.nanopore_assembly.fa | samtools sort - > $sample1.long_reads.bam
```

## Calling genotypes in the Nanopore data
As the assemblies are a consensus among sequences, we used a simple approach to identify the haplotypes along the reference coordinate system with *bcftools*. We used *bcftools mpileup* and piped the produced vcf file into a custom python script *get_genotypes_from_VCF.py* to get the genotypes along the genome.

```bash
bcftools mpileup -I -f guy11.fa $sample1.long_reads.bam | python3 get_genotypes_from_VCF.py --long /dev/stdin | gzip > $sample1.long_read_genotypes.gz
```


## Calling variant positions per sample

```bash
gatk -T HaplotypeCaller -R guy11.fa -I $sample1.short_mapped_reads.bam --genotyping_mode DISCOVERY -ERC BP_RESOLUTION -stand_call_conf 30 -ploidy -o $sample1.short_reads.vcf.gz
```

We used the custom script *get_genotypes_from_VCF.py* to get the genotypes along the genome.

```bash
python3 get_genotypes_from_VCF.py --short $sample1.short_reads.vcf.gz | gzip > $sample1.short_read_genotypes.gz
```

## Comparing matching and unmatching positions between short and long reads individual callings

We have created the python script *compare_short_vs_long_genotypes.py* which compares the generated genotypes and outputs a codified comparison per site in which:

Code | Description
---- | -----------------------------
A:   | Illumina 0 & Nanopore 0
B:   | Illumina 1 & Nanopore 1
C:   | Illumina 0 & Nanopore 1
D:   | Illumina 1 & Nanopore 0
E:   | Illumina X & Nanopore NA/abs
F:   | Illumina NA/abs & Nanopore X
N:   | Illumina NA & Nanopore NA


```bash
python3 compare_short_vs_long_genotypes.py $sample1.short_read_genotypes.gz $sample1.long_read_genotypes.gz 2> $sample1.summary.txt | gzip > $sample1.short_vs_long_compared.gz
```

The generated information can be summarized by exatracting the union of the matching and unmatching positions from all samples. This will constitute our gold standard variants dataset (GSVD) and the non-GSVD.

```bash
# All matching positions
zcat $sample{1..9}.short_vs_long.compared.gz | awk '($3 == "A" || $3 == "B")' | sort | uniq | gzip  > matching_positions.gz
# Matching SNPs
zcat $sample{1..9}.short_vs_long.compared.gz | awk '$3 == "B"' | sort | uniq | gzip  > GSVD.gz
# All mismatching positions
zcat $sample{1..9}.short_vs_long.compared.gz | awk '($3 == "C" || $3 == "D")' | sort | uniq | gzip > nonGSVD.gz
```

## Jointly calling genotypes on the 9 samples

To have information from every single base along the reference genome we use the flag *-allbases* in the *GATK* program *GenotypeGVCFs*.

```bash
gatk -T GenotypeGVCFs -R guy11.fa -V samples1-to-9.list -ploidy 1 -o joint.vcf.gz -allbases
```

## Filtering GSVD and non-GSVD from the Joint call

We now can filter the GSVD and non-GSVD positions from the joint SNP call.

```bash
bcftools view -R GSVD.gz joint.vcf.gz | bgzip > joint_GSVD.vcf.gz
bcftools view -R nonGSVD.gz joint.vcf.gz | bgzip > joint_nonGSVD.vcf.gz
```


## Extracting features

To ascertain the features and the empirical distributions of the two datasets, different summary statistics from the filtered SNP datasets are extracted

```bash
# We capture all the features in the VCF
zcat joint_GSVD.vcf.gz | grep "##INFO | awk -F "ID=" '{print $2}' | cut -f1 -d ","  | sed 's/^/--get-INFO /g' | tr "\n" " " > info_features

vcftools --gzvcf joint.goldstd.vcf.gz $(cat info_features) --min-alleles 2 --max-alleles 2 --remove-indels --stdout | cut -f5- > features.joint_GSVD.tsv
vcftools --gzvcf joint.mismatch.vcf.gz $(cat info_features) --min-alleles 2 --max-alleles 2 --remove-indels --stdout | cut -f5- > features.joint_nonGSVD.tsv
```

The result files can be easily imported to produce figures, inspect their distributions and generate cutoffs to be used as hard filters.

## SNP datasets with all rice-infecting isolates

Once the different cutoffs per filtering scheme were established, we called de-novo SNPs on the joint set of samples between the 9 samples previously described and the 131 samples described in [Latorre *et al*, 2020](https://doi.org/10.1186/s12915-020-00818-z). We created individual genomic VCF files for each isolate as described before and created a joint SNP call

```bash
gatk -T GenotypeGVCFs -R guy11.fa -V allsamples.list -ploidy 1 -o all_samples_joint.vcf.gz
```

Finally, we used the different filtering criteria to create four different SNP sets

```bash
# Latorre et al, 2020 filter
gatk -T VariantFiltration -R guy11.fa -V all_samples_joint.vcf.gz \
--filterExpression "QD < 5.0" --filterName "QDFilter" \
--filterExpression "QUAL < 5000.0" --filterName "QualFilter" \
--filterExpression "ReadPosRankSum < -2.0" --filterName "ReadPosRankSum" \
--filterExpression "ReadPosRankSum > 2.0" --filterName "ReadPosRankSum" \
--filterExpression "MQRankSum < -2.0" --filterName "MQRankSum" \
--filterExpression "MQRankSum > 2.0" --filterName "MQRankSum" \
--filterExpression "BaseQRankSum < -2.0" --filterName "BaseQRankSum" \
--filterExpression "BaseQRankSum > 2.0" --filterName "BaseQRankSum" \
--filterExpression "MQ < 20.0" --filterName "MQ" \
-o latorre_et_al_filtered.snps.vcf.gz

# QD-based filter
gatk -T VariantFiltration -R guy11.fa -V all_samples_joint.vcf.gz \
--filterExpression "QD < 22.21431 || QD > 32.16569" --filterName "QD" \
-o QD_filtered.snps.vcf.gz

# Relaxed filter
gatk -T VariantFiltration -R guy11.fa -V all_samples_joint.vcf.gz \
--filterExpression "QD < 19.72884 || QD > 34.65116" --filterName "QD" \
--filterExpression "ReadPosRankSum < -2.234249 || ReadPosRankSum > 2.234249" --filterName "ReadPosRankSum" \
--filterExpression "FS > 2.6" --filterName "FS" \
--filterExpression "MQRankSum < -0.5 || MQRankSum > 4.0" --filterName "MQRankSum" \
-o relaxed_filtered.snps.vcf.gz

# Stringent filter
gatk -T VariantFiltration -R guy11.fa -V all_samples_joint.vcf.gz \
--filterExpression "QD < 25.71431 || QD > 28.66569" --filterName "QD" \
--filterExpression "ReadPosRankSum < -0.0971292 || ReadPosRankSum > 0.0971292" --filterName "ReadPosRankSum" \
--filterExpression "FS > 0.1" --filterName "FS" \
-o stringent_filtered.snps.vcf.gz

```

The VCF files can be downloaded at [/data/VCF_SNP_sets.tar.gz](/data/VCF_SNP_sets.tar.gz)
