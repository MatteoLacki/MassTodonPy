from MassTodonPy import MassTodonize
from MassTodonPy.Formulator import make_formulas
from MassTodonPy.TestScripts.substancesP import substancesP

mol = substancesP[0] # WH = 0 WV = 300


formulas = make_formulas(
    fasta           = mol['fasta'],
    Q               = mol['precursorCharge'],
    modifications   = mol['modifications']    )
formulas = list(formulas.makeMolecules())


res  = MassTodonize(fasta           = mol['fasta'],
                    precursor_charge= mol['precursorCharge'],
                    mz_prec         = .1,
                    joint_probability_of_envelope = .999,
                    spectrum        = mol['spectrum'],
                    opt_P           = .99,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    max_times_solve = 10,
                    raw_data        = True,
                    for_plot        = True,
                    verbose         = False )

# [alpha for r in res['raw_estimates'] for alpha in r['alphas'] if alpha['molType'] == 'precursor']

from MassTodonPy.IsotopeCalculator import IsotopeCalculator

IC = IsotopeCalculator(verbose=True)
IC.isoEnvelope(
    atomCnt_str = 'C63H98N18O13S1',
    q = 1,
    g = 0
)

IC.isoEnvelope(
    atomCnt_str = 'C63H98N18O13S1',
    q = 1,
    g = 0
)


('precursor', 'C63H98N18O13S1', 11, 1, 0)
