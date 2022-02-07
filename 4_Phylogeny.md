# Rice blast variant calling
# 4. Effect of SNP filtered sets on the Phylogeny

Program                          | Link
-------------------------------- | -----------------------------------------
*bcftools* v.1.11                | (http://www.htslib.org/download/)
*PLINK* v.1.9                    | (https://www.cog-genomics.org/plink/1.9/)
*tped2nexus.py*                  | (https://github.com/Burbano-Lab/rice-blast-variant-calling/tree/main/scripts/)
*RAxML-NG*                       | (https://github.com/amkozlov/raxml-ng)
*root_to_tip_distances.py*       | (https://github.com/Burbano-Lab/rice-blast-variant-calling/tree/main/scripts/)
*BEAST* v.2.6.3                  | (https://github.com/CompEvol/beast2)
*compute_all_invariant_sites.py* | (https://github.com/Burbano-Lab/rice-blast-variant-calling/tree/main/scripts/)

First, we filtered out all individuals belonging to the recombinant genetic group. The analyses were carried out with full information positions. 
```bash
bcftools view -S ^recombinant_isolates.list -a -f PASS ALL.snps.vcf.gz | bcftools view -m2 -M2 -g ^miss - | bgzip > clonal.snps.vcf.gz
```
We transformed the output .vcf files into temporary transposed ped and fam files and converted them into a nexus format.
```bash
plink --aec --vcf-filter --geno 0 --recode transpose --vcf clonal.snps.vcf.gz --out clonal.snps

python3 tped2nexus.py clonal.snps > clonal.snps.nexus
```

In order to assess the presence of temporal signal, we first created ML-based phylogenies with RAxML-NG with 1,000 bootsrap replicates
```bash
raxml-ng --all --msa clonal.snps.fasta --model GTR+G --bs-trees 1000 --out clonal.snps
```
From the generated best trees we then calculated root-to-tip distances
```bash
python root_to_tip_distances.py clonal.snps.raxml.bestTree > clonal.snps.raxml.root-to-tip_distances.txt
```
We then correlate the phylogenetic distances with the collection dates
```r
#R#

dates <- read.table('clonal_collection_year.tsv', header = T)
RtoT_distances <- read.table('clonal.snps.raxml.root-to-tip_distances.txt' header = T)
dist_year <- merge(RtoT_distances, dates, by.x = 1, by.y = 1)

# Pearson's r correlation coefficient
cor.test(dist_year$Year, dist_year$Root_to_Tip)

# Scatter plot
plot(dist_year$Year, dist_year$Root_to_Tip, xlab = 'Collection Year', ylab = 'Root to Tip Distance')
abline(lm(dist_year$Root_to_Tip ~ dist_year$Year), col = 'red', lty = 2)

# Sampling by replacement and Permutation test
set.seed(12345)
cors <- c()
perms <- c()
for(i in 1:1000){
  sub <- dist_year[unique(sample(nrow(dist_year), replace = T)), ]
  ct <- cor(sub$Year, sub$Root_to_Tip, method = 'pear')
  cors <- c(cors, ct)
  reg <- lm(sub$Root_to_Tip ~ sub$Year)
  
  subperm <- cbind(dist_year$Root_to_Tip, dist_year[sample(1:nrow(dist_year), replace = F), 3])
  pt <- cor(subperm[,1], subperm[,2], method = 'pear')
  perms <- c(perms, pt)
}
boxplot(list(cors, perms), ylim = c(-0.45,0.45), names = c('Sampling by replacemen', 'Permutation'))

```

After assessing the presence of temporal signal, we can use the .nexus file as input for *beauti* *BEAST2* to produce a preliminary .xml file.
Finally, we calculated the amount of variant and invariant sites by calculating the intersection of the positions interrogated by every single individual .vcf file.
```bash
for sample in *.g.vcf.gz; do
	bcftools view -H $sample | cut -f1-2 | gzip > ${sample%.g.vcf.gz}.positions.gz
done

ls *.positions.gz > files_with_positions.list

# We parallelize the following script with 10 threads
python3 compute_all_invariant_sites.py files_with_positions.list clonal.snps.vcf.gz reference_genome.fasta 10 > invariant_sites.txt
```
A "constantSiteWeights=" tag was manually added to the xml file in order to incorporate this information in the phylogenetic reconstruction (see examples in the ./data/4_phtlogeny/ folder)

Finally, we submitted the .xml files into the CIPRES gateway (https://www.phylo.org/). Each dataset was independently run four times with a chain length of 10 million each.
