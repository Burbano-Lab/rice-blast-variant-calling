# USAGE: python3 tped2nexus.py <prefixfile.tped>

from sys import argv

tfam = argv[1] + '.tfam'
tped = argv[1] + '.tped'

ids = []
with open(tfam, 'r') as f:
    for line in f:
        ids.append(line.split(' ')[0]) 
sequences = {i:'' for i in range(len(ids))}
with open(tped, 'r') as f:
    for line in f: 
        for i in range(len(ids)): 
            column = (i * 2) + 4 
            sequences[i] += line.split(' ')[column] 

print('#NEXUS')
print('BEGIN DATA;')
print('DIMENSIONS NTAX={} NCHAR={};'.format(len(ids), len(sequences[0])))
print('FORMAT DATATYPE=DNA;')
print('MATRIX')
maxlenname = max([len(i) for i in ids])+1
for i in range(len(ids)):
    spaces = maxlenname - len(ids[i]) 
    print(ids[i], ' '*spaces, ''.join(sequences[i]), sep = '') 
print('\n;')
print('END;')
