%load_ext autoreload
%autoreload 2

from collections import defaultdict
from linearCounter.linearCounter import linearCounter as lCnt
import re

from MassTodonPy.Formulator.formulator import get_formulas
from MassTodonPy.Data.get_data import get_dataset, get_amino_acids

mol = get_dataset("substanceP")
mol['fasta']
amino_acids = get_amino_acids()

modifications = mol['modifications']
# modifications = {}
modifications2 = defaultdict(lambda: defaultdict(lambda: lCnt()),
                             {k - 1: defaultdict(lambda: lCnt(),
                                                 {name: lCnt(atomCnt)
                                                  for name, atomCnt
                                                  in v.items()})
                             for k, v in modifications.items()})

modifications2
modifications2[10]['N']
modifications2[10]['C_carbo']
modifications2[10]['C_alpha']
modifications2[9]['N']
modifications2[9]['C_carbo']
modifications2[9]['C_alpha']

mol['fasta']
formulas = get_formulas(fasta=mol['fasta'],
                        Q=mol['Q'],
                        what_fragments="cz",
                        modifications=mol['modifications'])
