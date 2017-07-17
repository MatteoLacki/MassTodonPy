import re

def parse_single_modification(mod):
    modification = {}
    for elem_no in mod.split(','):
        elem, no = elem_no.split(':')
        no = int(no)
        modification[elem] = no
    return modification

def parse_true_false(key, value):
    value = value.lower()
    if value in ('true', 'false'):
        return value == 'true'
    else:
        print('CONFIG FILE ERROR: value', value, 'not among possible logical value: True, False. This test is case-sensitive.')
        raise ValueError


def parse_plain_config_file(path):
    '''Parse a text file in a format similar to:

    fasta = RPKPQQFFGLM
    precursor_charge = 3
    modification C11 = H:1, O:-1, N:1
    cut_off = 100
    mz_prec = .05
    '''
    config_tmp = {}
    boolean_keys = ('forPlot', 'highcharts', 'raw_data', 'analyze_raw_data', 'verbose')
    float_keys = ('mz_prec', 'cut_off', 'opt_P', 'joint_probability_of_envelope',
        'min_prob_of_envelope_in_picking', 'L1_x', 'L2_x', 'L1_alpha', 'L2_alpha', 'output_deconvolution_threshold' )
    int_keys = ('precursor_charge', 'multiprocesses_No', 'max_times_solve')
    string_keys = ('fasta', 'frag_type', 'solver', 'method')


    with open(path, 'r') as f:
        for line in f:
            key, value = re.sub('[\s+]|[\n]','',line).split('=')
            config_tmp[key] = value
    config = {}
    modifs = {}
    for k in config_tmp:
        k_low = k.lower()
        if 'modification' == k_low[0:12]:
            aa_mod = k[12:]
            if not aa_mod in modifs:
                modifs[aa_mod] = parse_single_modification(config_tmp[k])
            else:
                print 'Watch out! There are at least two modifications on the same atom:', aa_mod
        elif k in boolean_keys:
            config[k] = parse_true_false(k, config_tmp[k])
        elif k in float_keys:
            config[k] = float(config_tmp[k])
        elif k in string_keys:
            config[k] = config_tmp[k]
        elif k in int_keys:
            config[k] = int(config_tmp[k])
        else:
            print('Key', k, 'not acceptable.')
            raise ValueError
    config['modifications'] = modifs
    return config


parse_plain_config_file('example_config_plain_text.txt')
