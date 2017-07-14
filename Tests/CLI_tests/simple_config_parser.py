def parse_single_modification(mod):
    modification = {}
    for elem_no in mod.split(','):
        elem, no = elem_no.split(':')
        no = int(no)
        modification[elem] = no
    return modification


def parse_plain_config_file(path):
    '''Parse a text file in a format similar to:

    fasta = RPKPQQFFGLM
    precursor_charge = 3
    modification C11 = H:1, O:-1, N:1
    cut_off = 100
    mz_prec = .05
    '''
    config_tmp = {}
    with open(path, 'r') as f:
        for line in f:
            key, value = re.sub('[\s+]|[\n]','',line).split('=')
            config_tmp[key] = value
    config = {}
    for k in config_tmp:
        if 'modification' == k[0:12]:
            value = parse_single_modification(config_tmp[k])
            k = k[12:]
        else:
            try:
                value = float(config_tmp[k])
            except ValueError:
                value = config_tmp[k]
        config[k] = value
    return config

parse_plain_config_file('example_config_plain_text.txt')
