import json
from pprint import pprint

with open('masstodon_input.json', 'r') as f:
    x = json.load(f)

pprint(x)
