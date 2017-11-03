from MassTodonPy.Formulator.formulator import make_formulas
from MassTodonPy.IsotopeCalculator.isotopeCalculator import IsotopeCalculator

fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
fasta = 'RPKPQQFFGLM'

F = make_formulas(fasta=fasta, Q=3)
precursor = list(F.makeMolecules())[0]

_, prec_str, K, q, g = precursor

IsoCalc = IsotopeCalculator(jP=.999, prec_digits=None)
spectrum = IsoCalc.isoEnvelope(prec_str)

mz, intensity = spectrum
max_intensity = max(intensity)
max_intensity_cnt = [i for i, j in enumerate(intensity) if j == max_intensity][0]
mz[max_intensity_cnt], max_intensity
