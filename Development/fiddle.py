%load_ext autoreload
%autoreload 2

import numpy as np
from collections import defaultdict
from MassTodonPy.Formulator.formulator import get_formulas
from MassTodonPy.Data.get_data import get_dataset, get_amino_acids
from MassTodonPy.IsotopeCalculator.isotopeCalculator import IsotopeCalculator
from MassTodonPy.Parsers.formula_parser import formulaParser
from MassTodonPy.Data.get_data import get_isotopic_masses_and_probabilities
from MassTodonPy.Spectra.operations import cdata2numpyarray,\
                                           aggregate_envelopes

from IsoSpecPy import IsoSpecPy

mol = get_dataset("substanceP")
mol['fasta']
amino_acids = get_amino_acids()
mol['modifications']

formulas = get_formulas(fasta=mol['fasta'],
                        Q=mol['Q'],
                        what_fragments="cz",
                        modifications=mol['modifications'])

atomCnt_str = "C100H200N10"
iso_masses, iso_probs = get_isotopic_masses_and_probabilities()


form_parser = formulaParser()

counts = []
isotope_masses = []
isotope_probs = []
atomCnt = form_parser.parse(atomCnt_str)

for el, cnt in atomCnt.items():
    counts.append(cnt)
    isotope_masses.append(iso_masses[el])
    isotope_probs.append(iso_probs[el])

joint_probability = .999
envelope = IsoSpecPy.IsoSpec(counts,
                             isotope_masses,
                             isotope_probs,
                             joint_probability)

masses, logprobs, _ = envelope.getConfsRaw()
masses = cdata2numpyarray(masses)
probs = np.exp(cdata2numpyarray(logprobs))

masses
probs

digits = 2
lists = defaultdict(list)
for mass, prob in zip(masses.round(digits), probs):
    lists[mass].append(prob)
help(lists)
from math import fsum
newMasses = np.array([k for k in lists])
newProbs = np.empty(len(newMasses))
for prob, mass in zip(np.nditer(newProbs, op_flags=['readwrite']),
                      newMasses):
    prob[...] = fsum(lists[mass])

iso_calc = IsotopeCalculator()
iso_calc.get_envelope(formulas[0], .999)
