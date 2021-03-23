'''
USAGE: python3 compare_short_vs_long_genotypes.py <illumina_genotypes.gz> <nanopore_genotypes.gz>

OUTPUT
STDOUT > Table per position
STDERR > Summary

Codes:
A: Illumina 0 & Nanopore 0
B: Illumina 1 & Nanopore 1
C: Illumina 0 & Nanopore 1
D: Illumina 1 & Nanopore 0
E: Illumina X & Nanopore NA/abs
F: Illumina NA/abs & Nanopore X
N: Illumina NA & Nanopore NA
'''

import gzip
from sys import argv, stderr

f1 = gzip.open(argv[1], 'rt')
f2 = gzip.open(argv[2], 'rt')

f1_line = f1.readline()
f2_line = f2.readline()
out = {'A':0, 'B':0, 'C':0,'D':0,'E':0,'F':0,'N':0}

while f1_line != '' and f2_line != '':
    if f1_line.split('\t')[0:2] == f2_line.split('\t')[0:2]:
        if f1_line.split('\t')[4] == f2_line.split('\t')[4]:
            if f1_line.split('\t')[4] == '0':
                print('\t'.join(f1_line.split('\t')[0:2]), 'A', sep = '\t')
                out['A'] += 1
            elif f1_line.split('\t')[4] == '1':
                print('\t'.join(f1_line.split('\t')[0:2]), 'B', sep = '\t')
                out['B'] += 1
            else:
                print('\t'.join(f1_line.split('\t')[0:2]), 'N', sep = '\t')
                out['N'] += 1
        else:
            if f1_line.split('\t')[4] == '0' and f2_line.split('\t')[4] == '1':
                print('\t'.join(f1_line.split('\t')[0:2]), 'C', sep = '\t')
                out['C'] += 1
            elif f1_line.split('\t')[4] == '1' and f2_line.split('\t')[4] == '0':
                print('\t'.join(f1_line.split('\t')[0:2]), 'D', sep = '\t')
                out['D'] += 1
            elif f2_line.split('\t')[4] == 'NA':
                print('\t'.join(f1_line.split('\t')[0:2]), 'E', sep = '\t')
                out['E'] += 1
            elif f1_line.split('\t')[4] == 'NA':
                print('\t'.join(f1_line.split('\t')[0:2]), 'F', sep = '\t')
                out['F'] += 1
        f1_line = f1.readline()
        f2_line = f2.readline()
    else:
        c1 = int(f1_line.split('\t')[0].split('|')[1][-2:]) # Extract contig number from assembly notation
        p1 = int(f1_line.split('\t')[1])
        c2 = int(f2_line.split('\t')[0].split('|')[1][-2:])
        p2 = int(f2_line.split('\t')[1])
        if c1 == c2:
            if p1 < p2: # Absent in f2
                while p1 < p2:
                    print('\t'.join(f1_line.split('\t')[0:2]), 'E', sep = '\t')
                    out['E'] += 1
                    f1_line = f1.readline()
                    p1 = int(f1_line.split('\t')[1])
            elif p2 < p1: # Absent in f1
                while p2 < p1:
                    print('\t'.join(f2_line.split('\t')[0:2]), 'F', sep = '\t')
                    out['F'] += 1
                    f2_line = f2.readline()
                    p2 = int(f2_line.split('\t')[1])
        else:
            if c1 < c2:
                while not(c1 >= c2 and p1 >= p2):
                    print('\t'.join(f1_line.split('\t')[0:2]), 'E', sep = '\t')
                    out['E'] += 1
                    f1_line = f1.readline()
                    c1 = int(f1_line.split('\t')[0].split('|')[1][-2:])
                    p1 = int(f1_line.split('\t')[1])
            else:
                while not(c2 >= c1 and p2 >= p1):
                    print('\t'.join(f1_line.split('\t')[0:2]), 'F', sep = '\t')
                    out['F'] += 1
                    f2_line = f2.readline()
                    c2 = int(f2_line.split('\t')[0].split('|')[1][-2:])
                    p2 = int(f2_line.split('\t')[1])
        f1_line = f1.readline()
        f2_line = f2.readline()

while f1_line != '':
    print('\t'.join(f1_line.split('\t')[0:2]), 'E', sep = '\t')
    out['E'] += 1
    f1_line = f1.readline()
while f2_line != '':
    print('\t'.join(f2_line.split('\t')[0:2]), 'F', sep = '\t')
    out['F'] += 1
    f2_line = f2.readline()

f1.close()
f2.close()

total = sum(out.values())
nonmiss = out['A'] + out['B'] + out['C'] + out['D']
print('''A: Illumina 0 & Nanopore 0:\t{}\t{}%
B: Illumina 1 & Nanopore 1:\t{}\t{}%
C: Illumina 0 & Nanopore 1:\t{}\t{}%
D: Illumina 1 & Nanopore 0:\t{}\t{}%
E: Illumina ANY & Nanopore N/ABS:\t{}\t{}%
F: Illumina N/ABS & Nanopore ANY:\t{}\t{}%
N: Illumina N/ABS & Nanopore N/ABS:\t{}\t{}%
-----------------------------------
From non-missing data
Matching:\t{}\t{}%
Non-matching:\t{}\t{}%'''.format(
out['A'], round((out['A']/total)*100, 2),
out['B'], round((out['B']/total)*100, 2),
out['C'], round((out['C']/total)*100, 2),
out['D'], round((out['D']/total)*100, 2),
out['E'], round((out['E']/total)*100, 2),
out['F'], round((out['F']/total)*100, 2),
out['N'], round((out['N']/total)*100, 2),
out['A']+out['B'], round(((out['A']+out['B'])/nonmiss)*100, 2),
out['C']+out['D'], round(((out['C']+out['D'])/nonmiss)*100, 2)), file = stderr)

