import json
from pprint import pprint

path='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/test_json_of_a_bitch.json'
with open(path) as data_file:
    data = json.load(data_file)

data['c'][0]
# Succcesss!!!
