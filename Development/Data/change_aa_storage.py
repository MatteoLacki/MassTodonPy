%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_data import get_amino_acids

AAs = get_amino_acids()
AAs
AAs2 = {(aa, brick): count
        for aa, counts in AAs.items()
        for brick, count in counts.items()}
AAs
AAs23 = [(aa_name,
         [(brick, list(count.items())) for brick, count in counts.items()])
         for aa_name, counts in AAs.items()]


# tuple(AAs2[('A', 'C_alpha')].items())
# AAs3 = [ tuple(counts) for  ]

import json

path = "/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/"

with open(path+"amino_acids.json", "w") as f:
    json.dump(AAs23, f)

amino_acids = AAs23


lCnt({atom: count for atom, count in amino_acids[0][1][0][1]})

from linearCounter.linearCounter import linearCounter as lCnt
amino_acids


{(aa_name, brick):  lCnt({a: c for a, c in atom_cnt})
 for aa_name, bricks in amino_acids
 for brick, atom_cnt in bricks}
