import json

from MassTodonPy.MassTodon import MassTodon


def get_name(key):
    return "_".join(map(str,key))


def get_subsequence(fasta, name):
    if name[0]=='p':
        return fasta
    if name[0]=='z':
        return fasta[len(fasta) - int(name[1:]):]
    else:
        return fasta[0:int(name[1:])]


def gen_data(deconvolution_results, fasta, Q):
    for r in deconvolution_results:
        for x in r['alphas']:
            name = x['molType']
            if not (name == 'precursor' and x['q'] == Q):
                f = {'seq':get_subsequence(fasta, name), 'Q':x['q'],'G':x['g'],'fragName':name, 'intensity':x['estimate'] }
                yield f


def perform_calculations(spectrum_path, output_path, config):
    '''Run MassTodonPy on a given spectrum and save the results into a given output folder.

    Parameters
    ----------
    spectrum_path : str
        Path to the spectrum in a mzXml, mzml, or tsv or csv format.

    output_path : str
        Path to the output.

    config : dict
        A dictionary with the parsed configuration of MassTodon.
    '''
    print(config)
    # if 'csv' in config:
    #     del config['csv']
    #     config['output_csv_path'] = output_path
    #
    # results = MassTodon(spectrum=spectrum_path, **config)
    # #TODO how to output the Bokeh results?
