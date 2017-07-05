import  json

path_to_data = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/developing_csi/make_config_file.json'

config = {  'fasta': 'RPKPQQFFGLM',
            'Q': 3,
            'modifications': { 'C11':{'H':1,'O':-1,'N':1}} }

with open(path_to_data, 'w') as fp:
    json.dump(config,fp)
