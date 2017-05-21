from    MassTodonPy.CSI import run_masstodon
import  json
import  pandas as pd


spectrum_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/FRL-010513-SUBP-WH000-WV300.txt'

config_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/developing_csi/make_config_file.json'

with open(config_path, 'r') as f:
    config = json.load(f)


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

def save_results(spectrum_path, config):
    file_path, file_name, deconvolution_results, results_analyzed = run_masstodon(
        spectrum_path = spectrum_path, **config)

    with open(file_path+'/outputs/'+'probs.json', 'w') as f:
        json.dump(results_analyzed,f)

    ParsedOutput = pd.DataFrame(gen_data(deconvolution_results, config['fasta'], config['Q']))[[2,1,0,3,4]]

    ParsedOutput.to_csv(path_or_buf = file_path + "/outputs/"+ file_name + ".csv",
                        index   = False,
                        sep     = '\t' )
