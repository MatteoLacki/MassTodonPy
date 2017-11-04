from MassTodonPy import MassTodonize, get_data

substanceP = get_dataset('substanceP')
ubiquitin  = get_dataset('ubiquitin')
mol = substanceP.copy()

res = MassTodonize(fasta                           = mol['fasta'],
                   precursor_charge                = mol['Q'],
                   mz_prec                         = .05,
                   joint_probability_of_envelope   = .999,
                   modifications                   = mol['modifications'],
                   spectrum                        = mol['spectrum'],
                   opt_P                           = .99,
                   solver                          = 'multiprocessing',
                   multiprocesses_No               = None,
                   distance_charges                = 5,
                   max_times_solve                 = 10,
                   raw_data                        = True,
                   output_csv_path                 = '/Users/matteo/Documents/MassTodon/results/',
                   verbose                         = False )

res['summary']
