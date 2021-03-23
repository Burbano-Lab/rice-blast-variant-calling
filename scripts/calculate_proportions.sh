vcf=$1

nsnps=$(bcftools view -H $vcf | wc -l)

fixed=$(bcftools view $vcf -H | cut -f8 | awk -F "AC=" '{print $2}' | cut -f1 -d ";" | sort | uniq -c | sort -k2,2nr | head -n 1 | awk '{print $1}')
singletons=$(bcftools view $vcf -H | cut -f8 | awk -F "AC=" '{print $2}' | cut -f1 -d ";" | sort | uniq -c | sort -k2,2n | head -n 1 | awk '{print $1}')
interm=$(python -c "print($nsnps - ($fixed + $singletons))")

p_fixed=$(python -c "print( round($fixed / $nsnps, 2))")
p_interm=$(python -c "print(round($interm / $nsnps, 2))")
p_singletons=$(python -c "print(round($singletons / $nsnps, 2))")

echo -e "Total_Snps\t$nsnps\t1.0"
echo -e "Fixed\t$fixed\t$p_fixed"
echo -e "Segreg\t$interm\t$p_interm"
echo -e "Singletons\t$singletons\t$p_singletons"
