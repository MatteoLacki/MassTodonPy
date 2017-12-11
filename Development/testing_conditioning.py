%load_ext autoreload
%autoreload 2

import numpy as np

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Formula.Formula import Formula
from MassTodonPy.Molecule.Molecule import Molecule

N = 5
base = {'C':378, 'H':628, 'N':105, 'O':118, 'S':1}
res = []
mz_precision = 3
shift = .5

for i in range(N):
    base['H'] += 1
    formula = Formula(base)
    mol = Molecule('ubiquitin', 'Michal', formula, 1)
    iso = mol.isotopologues(mz_precision=mz_precision+1)
    iso.mz += shift * 10**(-mz_precision)
    iso.round_mz(mz_precision)
    if mz_precision:
        iso.mz = np.round(iso.mz * 10**mz_precision).astype('int')
    else:
        iso.mz = iso.mz.astype('int')
    res.append(iso)

mzs = sorted(set(mz for r in res for mz, i in r))

M = []
for r in res:
    X = []
    for mz in mzs:
        L = list(r[mz-.1, mz+.1])
        if L:
            X.append(L[0][1])
        else:
            X.append(0.0)
    M.append(X)

X = np.matrix(M).T


print(np.linalg.cond(X.T * X))
