# USGAGE: python3 compute_all_invariant_sites.py <files_with_positions.list> <file.vcf> <reference.fasta> <N_threads>

import gzip
import vcfpytools
from Bio import SeqIO
from sys import argv
from itertools import islice
from multiprocessing import Pool

#### Parse argv
files = []
with open(argv[1], 'r') as f:
    for line in f:
        files.append(line.rstrip())

vcf_file = argv[2]
reference_file = argv[3]
threads = int(argv[4])

#### Functions
def compute(filelst):
    common = []
    with gzip.open(filelst[0], 'rt') as f:
            for line in f:
                common.append('_'.join(line.rsplit()))
    common = set(common)
    if len(filelst) > 1:
        rest = filelst[1:]
        for fl in rest:
            actual = []
            with gzip.open(fl, 'rt') as f:
                for line in f:
                    actual.append('_'.join(line.rsplit()))
            common = common.intersection(set(actual))
    return(common)

#### Define work load per thread

lst = iter(files)
filesall = list(iter(lambda: tuple(islice(lst, len(files) // threads)), ()))

#### All intersections (parallel)

p = Pool(threads)
out = p.map(compute, filesall)

#### Final intersections (single thread)

common = out[0]
for i in range(1,len(out)):
    common = common.intersection(out[i])

print('Common', len(common))
#### Substracting SNPs
snps = []
for record in vcfpytools.get_body(vcf_file):
    if record.rsplit()[6] == 'PASS':
        snps.append('_'.join(record.rsplit()[0:2]))

invar_sites = common.difference(snps)
print('Invar_Sites (Common - SNPs): ', len(invar_sites))

#### Assessing nucleotides from final list
reference = {}
for contig in SeqIO.parse(reference_file, 'fasta'):
    reference[contig.id] = contig.seq

output = {'A':0, 'C':0, 'G':0, 'T':0, 'N':0}
for record in invar_sites:
    contig = record.split('_')[0]
    pos = int(record.split('_')[1]) - 1 # Convert to python-based 0 values
    output[reference[contig][pos]] += 1

print('Base\tCount')
for i in output:
    print(i, output[i])
print('Ordered tag: ', ' '.join([str(output[i]) for i in ('A', 'C', 'G', 'T')]))
