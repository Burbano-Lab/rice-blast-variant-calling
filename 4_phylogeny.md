# Rice blast variant calling
# 4. Effect of SNP filtered sets on the Phylogeny

Program                          | Link
-------------------------------- | -----------------------------------------
*bcftools* v.1.11                | (http://www.htslib.org/download/)
*PLINK* v.1.9                    | (https://www.cog-genomics.org/plink/1.9/)
*tped2nexus.py*                  | (./scripts/tped2nexus.py)
*BEAST* v.2.6.3                  | (https://github.com/CompEvol/beast2)
*compute_all_invariant_sites.py* | (./scripts/compute_all_invariant_sites.py)

First, we filtered out all individuals belonging to the recombinant genetic group. The analyses were carried out with full information positions. 
```bash
bcftools view -S ^recombinant_isolates.list -a -f PASS ALL.snps.vcf.gz | bcftools view -m2 -M2 -g ^miss - | bgzip > clonal.snps.vcf.gz
```
We transformed the output .vcf files into temporary transposed ped and fam files and converted them into a nexus format.
```bash
plink --aec --vcf-filter --geno 0 --recode transpose --vcf clonal.snps.vcf.gz --out clonal.snps

python3 tped2nexus.py clonal.snps > clonal.snps.nexus
```

Now, the .nexus file is suitable to be loaded into the *beauti* program from *BEAST2* to produce a preliminary .xml file.
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
