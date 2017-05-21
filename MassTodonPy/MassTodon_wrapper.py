from    MassTodonPy  import MassTodon
from    time import time

def run_masstodon(  fasta, Q, modifications, spectrum_path,
                    jP      =.999,
                    mzPrec  =.05,
                    precDigits  = 2,
                    M_minProb   = .7,
                    cutOff      = 100.,
                    cutOff2     = 0.0,
                    topPercent  = .999,
                    solver      = 'sequential',
                    solver_mode = 'MSE',
                    solver_max_T= 30,
                    L1_x=0.001,
                    L2_x=0.001,
                    L1_alpha=0.001,
                    L2_alpha=0.001,
                    verbose=False ):
    params = (fasta, Q, modifications, spectrum, jP, mzPrec, precDigits, M_minProb, cutOff, topPercent, max_times_solve, L1_x, L2_x, L1_alpha, L2_alpha)
    try:
        M = MassTodon(  fasta           = fasta,
                        precursorCharge = Q,
                        precDigits      = precDigits,
                        jointProbability= jP,
                        mzPrec          = mzPrec,
                        modifications   = modifications  )

        M.readSpectrum( path        = spectrum_path,
                        cutOff      = cutOff,
                        digits      = precDigits,
                        topPercent  = topPercent  )

        M.prepare_problems(M_minProb)
        T0_deconv = time()
        deconvolution_results = M.run(
                        solver          = solver,
                        method          = solver_mode,
                        max_times_solve = solver_max_T,
                        L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                        verbose = verbose )

        T1_deconv = time()
        T_deconv  = T1_deconv - T0_deconv
        results_analyzed = {}
        try:
            results_analyzed['base']    = M.analyze_reactions(analyzer='basic')
        try:
            results_analyzed['inter']   = M.analyze_reactions(analyzer='inter')
        try:
            results_analyzed['up_inter']= M.analyze_reactions(analyzer='up_inter')
        if verbose:
            res = deconvolution_results, results_analyzed, params
        else:
            res = deconvolution_results, results_analyzed
    except Exception as e:
        res = params, e
    return res

run_masstodon()
