from MassTodonPy.Precursor.Precursor import Precursor

substanceP_fasta = 'RPKPQQFFGLM'
ubiquitin_fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'


subPprec = Precursor(substanceP_fasta, 3)
subPmols = list(subPprec.molecules())
len(subPmols)

ubi_prec = Precursor(ubiquitin_fasta, 9)
ubi_mols = list(ubi_prec.molecules())
len(ubi_mols)
