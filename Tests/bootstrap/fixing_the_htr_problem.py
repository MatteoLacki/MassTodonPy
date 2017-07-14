import  os
os.environ['OMP_NUM_THREADS'] = "1"

from MassTodonPy.TestScripts.substancesP import substancesP
from MassTodonPy import MassTodonize
import cPickle as pickle

mol = None
for ID, mol in enumerate(substancesP):
    WH, WV = [mol['experimental_setting'][x] for x in ('WH','WV')]
    if WH == 150 and WV == 400:
        mol2 = mol

R = MassTodonize(   fasta       = mol2['fasta'],
                precursor_charge= mol2['precursorCharge'],
                mz_prec         = .065,                     # mz_prec = .05
                joint_probability_of_envelope= .999,
                modifications   = mol2['modifications'],
                spectrum        = mol2['spectrum'],
                opt_P           = .99, # cut_off = 500.0,
                solver          = 'multiprocessing',# solver  = 'sequential'
                multiprocesses_No = None,
                max_times_solve = 10,
                raw_data        = True,
                analyze_raw_data= False,
                verbose         = False )

from MassTodonPy.MatchMaker import match_cz_ions
R['raw_estimates'][0]['alphas']



match_cz_ions(  results_to_pair         = R['raw_estimates'],
                Q                       = mol['precursorCharge'],
                fasta                   = mol['fasta'],
                advanced_args           = {},
                min_acceptEstimIntensity= 0.0,
                analyzer                = 'intermediate',
                accept_nonOptimalDeconv = False,
                verbose                 = False  )

mol['fasta']
