_preprocessing_args =   spectrum_minimal_intensity=eps
                        spectrum_percent_top_peaks=1.0
                        
peak_picking_args = mz_tol
                    min_prob_per_molecule
                    m_over_z_precision=0.05,

   """Run a full session of MassTodon on your problem.

    Parameters
    ==========
    precursor_fasta : str
        The fasta of the input substance.
    precursor_charge : int
        The charge of the input substance.
    precursor_name : str
        The name of the input substance.
    precursor_modifications :
        The key in the modifications' dictionary should with of form
        **<N|Calpha|C><amino acid No>**, e.g. C10.
        The values contain dictionaries with diffs in elements,
        e.g. **{'H': 1, 'N': 1, 'O': -1}**. TODO: check it.
    m_over_z_precision : float
        The precision in the m/z domain: how far away [in Daltons] should
        one search for peaks from the theoretical mass.
    fragmentation_type : str
        Types of fragments you want to find.
        Warning
        =======
        Supports c and z fragmentations. So you can input "cz", "c", or "z".
    blockedFragments : set of strings
        Which fragments should not appear. Defaults to 'c0'.
        MoleculeMaker adds to this all proline-blocked fragments.
    minimal_distance_between_charges: int
        The minimal distance between charges on the protein.
        If set to 5, at least 4 amino acids must la
        between consecutive *charged* amino acids.
    isotopologues_joint_probability : float
        The joint probability threshold for generating
        theoretical isotopic distributions.
    min_prob_of_envelope_in_picking : float
        The minimal probability an envelope has to scoop
        to be included in the deconvolution graph.
    spectrum : path string or a tuple of numpy arrays (masses, intensities).
    spectrum_minimal_intensity : float
        Experimental peaks with lower height will be trimmed.
    spectrum_percent_top_peaks : float
        The percentage of the heighest experimental peaks used in the analysis.
    _max_times_solve : int
        How many times the CVXOPT solver should try to find the optimal
        deconvolution of the isotopic envelopes.
        Note
        ====
        CVXOPT sometimes does not meet its optimality criteria.
        However, it uses underneath the OPENBLAS library and this baby has
        its own dynamic scheduling that might sometimes influence the output
        of the calculations.
        Here, we try to use it to our advantage, in order to help CVXOPT
        to meet its optimisation criteria.
        This is clearly not the nice way to do it, but we will soon go bayesian
        anyway and use a different approach to perform the deconvolution.
    _isotope_masses : dict
        The isotopic masses, e.g. from IUPAC.
    _isotope_probabilities : dict
        The isotopic frequencies, e.g. from IUPAC.
    _fitting_criterion : str
        The minimization criterion used.
        Warning
        =======
            For now, only mean square error (MSE) considered.
    _L1_x : float
        The $L_1$ penalty for flows (intensities attributed to particular
        isotopologous and experimental groups) in the deconvolution problem.
    _L2_x : float
        The $L_2$ penalty for flows (intensities attributed to particular
        isotopologous and experimental groups) in the deconvolution problem.
    _L1_alpha : float
        The $L_1$ penalty for total intensities attributed to particular
        molecular species in the deconvolution problem.
    _L2_alpha : float
        The $L_2$ penalty for total intensities attributed to particular
        molecular species in the deconvolution problem.
    _verbose : boolean
        Do you want to check what MassTodon did undercover?
    Returns
    =======
    out : a results object
    """