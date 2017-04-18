import  json

path_to_data = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/substanceP_results.json'

with open(path_to_data, 'r') as fp:
    data = json.load(fp)
