%load_ext autoreload
%autoreload
# from MassTodon import MassTodon
import numpy as np

fasta='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 9; modifications = {}; ionsNo  = 10000; P = .999


from Formulator import makeFragments

fragmentator = makeFragments(fasta, Q)
list(fragmentator.makeMolecules())
