# USAGE:
# python3 root_to_tip_distances.py <tree.newick>

from sys import argv
from Bio import Phylo

tree = Phylo.read(argv[1], format = 'newick')
tree.root_at_midpoint() # Root tree to Midpoint

RtoT = {}
for sample in tree.get_terminals():
    path = tree.get_path(sample)
    RT = sum([path[i].branch_length for i in range(len(path))])
    RtoT[sample.name] = RT

print('Sample\tRoot_to_Tip')
for i in RtoT:
    print(i, RtoT[i], sep = '\t')

