import numpy as np
import json
from collections import namedtuple
try:
    import cPickle as pickle
except ImportError:
    import pickle

path = "/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/"


with open(path+"amino_acids.pickle", "rb") as f:
    bricks = pickle.load(f)


# class AminoAcid(object):
#     """Information on amino acids."""
#
#     __slots__ = ("N", "C_alpha", "C_carbo")
#
#     def __init__(self, N, C_alpha, C_carbo):
#         self.N = N
#         self.C_alpha = C_alpha
#         self.C_carbo = C_carbo
#
#     def __str__(self):
#         return "N\t= {}\nC_alpha\t= {}\nC_carbo\t= {}".format(self.N,
#                                                               self.C_alpha,
#                                                               self.C_carbo)
# new_bricks = {el: AminoAcid(N=v["L"], C_alpha=v["C"], C_carbo=v["R"])
#               for el, v in bricks.items()}

remapping_keys = {"L": "N", "C": "C_alpha", "R": "C_carbo"}
amino_acids = {el: {remapping_keys[aa_names]: atom_cnt
                    for aa_names, atom_cnt in v.items()}
               for el, v in bricks.items()}



# amino_acids = {el: N=v["L"], C_alpha=v["C"], C_carbo=v["R"]
#                for el, v in bricks.items()}
#
# AminoAcid = namedtuple("AminoAcid", ("N", "C_alpha", "C_carbo"))
# amino_acids = {el: AminoAcid(N=v["L"], C_alpha=v["C"], C_carbo=v["R"])
#                for el, v in bricks.items()}

with open(path + "amino_acids2.pickle", "wb") as f:
    pickle.dump(amino_acids, f, protocol=2)


import cPickle as pickle
from collections import namedtuple
AminoAcid = namedtuple("AminoAcid", ("N", "C_alpha", "C_carbo"))

path = "/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/"
with open(path + "amino_acids2.pickle", "rb") as f:
    amino_acids = pickle.load(f)

print amino_acids
print amino_acids["A"]



from MassTodonPy.AminoAcid.aminoAcid import AminoAcid

with open(path+"amino_acids2.pickle", "rb") as f:
    bricks = pickle.load(f)
