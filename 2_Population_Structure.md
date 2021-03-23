# Rice blast variant calling
# 2. Effect of SNP filtered datasets on population structure


Program              | Location
-------------------- | -----------------------------------------
*PLINK* v.1.9        | (https://www.cog-genomics.org/plink/1.9/)
*AdmixTools* v7.0    | (https://github.com/DReichLab/AdmixTools)
*vcftools* v.0.1.16  | (https://github.com/vcftools/vcftools)
*R v.4.0.3*          | (https://www.r-project.org/)


We used two main approaches to evaluate the impact of the different filters on the rice-infecting blast population structure analyses: i) f3-outgroup based method and, ii) PCA-based method

## f3-outgroup based method
We used *PLINK* to convert the .vcf files to .ped and .map format.

```bash
# We restricted these analyses for full information positions

plink --allow-extra-chr --geno 0 --recode --vcf-filter --vcf $input.vcf --out $output_prefix
```

We used *qp3Pop* from *AdmixTools* to compute all f3-outgroup tests in the following 3-taxa configuration:

*f3*(A, B; Outgroup)

; where A and B represent all possible pairwise comparisons and the Outgroup is the wheat-infecting isolate BTJP-4(12) in all instances.

An example of the input files is shown:

*isolates.list*
```bash
head -n 5 isolates.list

#BD0024	U	BD0024
#BR0026	U	BR0026
#CD0073	U	CD0073
#CD0203	U	CD0203
#CH0052	U	CH0052
```

*comparisons.txt*
```bash
head -n 5 comparisons.txt

#BD0024 BR0026 O-6-BTJP4-12-S99
#BD0024 CD0073 O-6-BTJP4-12-S99
#BD0024 CD0203 O-6-BTJP4-12-S99
#BD0024 CH0052 O-6-BTJP4-12-S99
#BD0024 CH0063 O-6-BTJP4-12-S99
```

*f3outgroup.par*
```bash
cat f3outgroup.par
# genotypename: $subset1.ped
# snpname: $subset1.map
# indivname: isolates.list
# popfilename: comparisons.txt
```

We ran *qp3Pop*
```bash
qp3Pop -p f3outgroup_$subset1 > f3outgroup_$subset1.out
```

The output file contains the calculates f3 values as well as it associates Z-scores.

## PCA-based method

We used *PLINK* to calculate Hamming distances from the .vcf file
```bash
plink --allow-extra-chr --geno 0 --distance square --geno 0 --vcf-filter --vcf $input.vcf --out $output_prefix
```

PCA was computed using the *R* function *prcomp*

```R
#R
distances <- read.table('$output_prefix.dist')
pca <- prcomp(distances, scale. = TRUE)
```
The robustness of clusters was ascertained using Silhouette Scores with the *R* package *cluster*
```R
library(cluster)

num_retained_pcs <- 10 # Number of retained PCs might change depending on the cummulative explained variance
coords_PCA <- pca$x[,x:num_retained_pcs]

for(number_of_clusters in 2:7){ # Iterates from the values 2 to 7
	clstr <- pam(coords_PCA, number_of_clusters, metric = "euclidean")
}
```
Finally, in order to generate confidence intervals, we sampled with replacement and recalculated the Silhouette scores
```R
iterations <- 1000 # Number of iterations
for(i in 1:iterations){
	# As the coordinates are kept in an object with dimensions i x j ;
	# where i is the number of individuals and j is the number of PCs, the sampling will be done row-ise
	sampled_coordinates <- coords_PCA[sample(nrow(coords_PCA), nrow(coords_PCA), replace = True), ]
	clstr_sampled <- pam(sampled_coordinates, number_of_clusters, metric = "euclidean")
}
```

