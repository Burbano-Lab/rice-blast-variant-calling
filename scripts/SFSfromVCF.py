# USAGE:
# Unfolded: python3 SFSfromVCF.py <file.vcf|.gz> <samples.txt> <outgroup1> ... <outgroup_n>
# Folded: python3 SFSfromVCF.py <file.vcf|.gz> -f <samples.txt>
# 2-Dimension SFS: python3 SFSfromVCF.py <file.vcf|.gz> -2d <samples1.txt> <samples2.txt> <outgroup1> ... <outgroup_n>

import vcfpytools
from numpy import arange
from sys import argv

vcf = argv[1]

if '-f' in argv:
    folded = True
    samples_infile = argv[-1]
    ancestral = []
else:
    folded = False
    samples_infile = argv[2]
    ancestral = argv[3:]

def SFS(vcf, samples, folded, ancestral):
    Nancestral = len(ancestral)
    freqs = []
    if folded == True:
        for position in vcfpytools.get_genotypes_hap(vcf, samples, binary=True):
            if '.' not in position:
                geno = ((position.count('0')), (position.count('1')))
                freqs.append(min(geno))
    else:
        for position in vcfpytools.get_genotypes_hap(vcf, samples + ancestral, binary = True):
            if '.' not in position: # Full info sites
                ancest_alleles = position[-Nancestral:]
                if ancest_alleles.count('0') == Nancestral:
                    freqs.append(position[:-Nancestral].count('1'))
                elif ancest_alleles.count('1') == Nancestral:
                    freqs.append(position[:-Nancestral].count('0'))
    if folded == True:
        out = {i:freqs.count(i) for i in range(1,int(len(samples)/2)+1)}
    else:
        out = {i:freqs.count(i) for i in range(1,len(samples))} #DAF 0 and 1 non-informative
    return(out)

def ESFS(segsites, numsamples, folded):
    '''This function estimates the Expected Site Frequency Spectrum (ESFS).
    ESFS(i) = theta / i ; where i is the given frequency (goes from 1 to n-1)
    theta can be replaced by its Watterson's estimator* (Watterson, G.A. 1975).
    * Assumes WF population at equilibrium and infinite sites mutation model)
    '''
    Eout = {i:0 for i in range(1,numsamples)}
    a = sum([1/i for i in range(1,numsamples)])
    for i in range(1,numsamples):
        Eout[i] = segsites / (i * a)
    if folded == True:
        Eout = {i:Eout[i] + Eout[len(samples)-i] for i in range(1, int(len(samples) * 0.5)+1)}
    return(Eout)

with open(samples_infile, 'r') as f:
    samples = [line.strip() for line in f]

out = SFS(vcf, samples, folded, ancestral)
total = sum(list(out.values()))
Eout = ESFS(total, len(samples), folded)

step = 0.1 # Change here the bin size
if folded == True:
    bins = arange(0,0.5,step)
else:
    bins = arange(0,1,step)
print('Bin_Start\tBin_End\tExp_Prop_Deriv_Sites\tObs_Prop_Deriv_Sites')
for bn in bins:
    binvals = []
    binEvals = []
    for count in out:
        if count / len(samples) > bn and count / len(samples) <= bn + step:
            binvals.append(out[count])
            binEvals.append(Eout[count])
    print(round(bn, 3), round(bn+step, 2), sum(binEvals) / sum(list(Eout.values())), sum(binvals) / sum(list(out.values())), sep = '\t')
