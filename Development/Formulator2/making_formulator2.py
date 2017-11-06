%load_ext autoreload
%autoreload 2

from collections import defaultdict
from linearCounter.linearCounter import linearCounter as lCnt
import re

from MassTodonPy.Formulator.formulator import get_formulas
from MassTodonPy.Data.get_data import get_dataset, get_amino_acids

mol = get_dataset("substanceP")
amino_acids = get_amino_acids()

modifications = mol['modifications']
# modifications = {}
modifications = defaultdict(lambda: defaultdict(lCnt),
                            {k: {name: lCnt(atomCnt)
                                 for name, atomCnt in v.items()}
                             for k, v in modifications.items()})

# if 11 in modifications:
#     print(modifications[F.11])

formulas = get_formulas(fasta=mol['fasta'],
                        Q=mol['Q'],
                        what_fragments="cz",
                        modifications=mol['modifications'])
