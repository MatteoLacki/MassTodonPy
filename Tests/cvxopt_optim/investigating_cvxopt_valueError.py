import cPickle as pickle
from collections import Counter
import json


with open('/home/matteo/masstodon/repeatability_tests/MassTodonPy/errors_raw.matteo', 'r') as f:
    res = pickle.load(f)

X = res[(31, 150, 3000)]
X = [r for r in X if r['status'] != 'optimal'][0]

SG = X['SG']
SG.nodes(data=True)


Counter( SG.node[k]['type']  for k in SG )




mols = [ ([ SG.node[I]['mz'] for I in SG[M] ], [ SG.node[I]['intensity'] for I in SG[M] ], "_".join([ str(SG.node[M][x]) for x in ('molType', 'q', 'g', 'formula')]) ) for M in SG if SG.node[M]['type'] == 'M' ]

for N in SG:
    if SG.node[N]['type'] == 'M':
        mol['info'] = SG.node[N]
        mol['mz'] = [ I['mz'] for I in SG[N] ]
        mol['intensity'] = [ I['intensity'] for I in SG[N] ]
    

mz = [ (SG.node[N]['max_mz']+SG.node[N]['min_mz'])/2.0 for N in SG if SG.node[N]['type'] == 'G' ]
intensity = [ SG.node[N]['intensity'] for N in SG if SG.node[N]['type'] == 'G' ]


with open('error_structures', 'w') as f:
    json.dump([mols, mz, intensity], f)
