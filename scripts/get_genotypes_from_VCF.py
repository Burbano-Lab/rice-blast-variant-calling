# 'USAGE: python3 get_genotypes_from_VCF.py <vcf_type> <input_file.vcf.gz>
#  vcf_type can be \'-s\' or \'--short\' for short read-based vcf files
#  or \'-l\'  or \'--long\' for long read-based vcf files
#
#  OUTPUT: REF    ALT    genotype_bin    genotype_base

from sys import argv, exit
import vcfpytools

vcf = argv[-1]
if '-s' in argv or '--short' in argv:
    vcftype = 'short'
elif '-l' in argv or '--long' in argv:
    vcftype = 'long'
else:
    exit('Specify \'-s\'hort reads or \'-l\'ong reads')

for record in vcfpytools.get_body(vcf):
    contig = record.split('\t')[0]
    position = record.split('\t')[1]
    REF = record.split('\t')[3][0]
    ALTER = record.split('\t')[4].split(',')
    if vcftype == 'short':
        DEPTH = int(record.split('\t')[9].split(':')[2])
        if DEPTH == 0:
            ALT = 'NA'
            genot = 'NA'
            print(contig, position, REF, ALT, genot, ALT, sep = '\t')
            continue
    else:
        DEPTH = int(record.split('\t')[7].split(';')[0].split('=')[1])
        if DEPTH != 1:
            ALT = 'NA'
            genot = 'NA'
            print(contig, position, REF, ALT, genot, ALT, sep = '\t')
            continue
    if len(ALTER) == 1:
        ALT = REF
        genot = 0
    else:
        if len(ALTER) == 2:
            if len(ALTER[0]) == 1:
                ALT = ALTER[0]
                genot = 1
            else:
                ALT = 'NA'
                genot = 'NA'
        else:
            ALT = 'NA'
            genot = 'NA'
    print(contig, position, REF, ALT, genot, ALT, sep = '\t')

