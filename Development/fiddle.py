%load_ext autoreload
%autoreload 2

import numpy as np
from collections import defaultdict
from MassTodonPy.MoleculeMaker.MoleculeMaker import get_molecules
from MassTodonPy.Data.get_data import get_dataset, get_amino_acids
from MassTodonPy.IsotopeCalculator.isotopeCalculator import IsotopeCalculator
from MassTodonPy.Parsers.formula_parser import parse_formula
from MassTodonPy.Data.get_data import get_isotopic_masses_and_probabilities
from MassTodonPy.Spectra.operations import cdata2numpyarray,\
                                           aggregate_envelopes

from IsoSpecPy import IsoSpecPy


iso_masses, iso_probs = get_isotopic_masses_and_probabilities()
iso_probs['H']
iso_masses['H']

mol = get_dataset("substanceP")
amino_acids = get_amino_acids()
molecules = get_molecules(fasta=mol['fasta'],
                          Q=mol['Q'],
                          what_fragments="cz",
                          modifications=mol['modifications'])

# parse_formula(molecules[0].formula)

iso_calc = IsotopeCalculator()
iso_calc.get_envelope(molecules[0], .999)
iso_calc.get_monoisotopic_mz(molecules[0])
iso_calc.get_mean_mz(molecules[0])
iso_calc.get_mz_sd(molecules[0])

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

from math import fsum
newMasses = np.array([k for k in lists])
newProbs = np.empty(len(newMasses))
for prob, mass in zip(np.nditer(newProbs, op_flags=['readwrite']),
                      newMasses):
    prob[...] = fsum(lists[mass])
