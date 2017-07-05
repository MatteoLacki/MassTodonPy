from    MassTodonPy  import MassTodonize
import  json

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


def perform_calculations(spectrum_path, output_path, file_name, config):
    '''Run MassTodonPy on a given spectrum into a given output folder.'''

    # todo: config should not contain an entry
    try:
        results = MassTodonize(spectrum_path= spectrum_path,
                               forPlot      = False,
                               highcharts   = True,
                               verbose      = False,
                               max_times_solve = 10,
                               **config )

        highcharts = results['highcharts']
        del results['highcharts']

        print output_path+'output.json'
        with open(output_path+'output.json', 'w') as f:
            json.dump(results,f)
    except Exception as e:
        print e
        print config
