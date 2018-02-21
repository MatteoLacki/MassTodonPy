import json

from MassTodonPy.Parsers.blocked_fragments import parse_blocked_fragments

def parse_json(path):
    """Parse a json settings file for MassTodon."""
    with open(path, 'r') as f:
        masstodon_args = json.load(f)

    # Renaming
    if u'deconvolution_method' in masstodon_args:
        method = masstodon_args.pop('deconvolution_method')
        masstodon_args['method'] = method

    # Parsing modifications
    modifications = masstodon_args.pop('modifications')
    for k in modifications.keys():
        v = modifications.pop(k)
        modifications[int(k)] = v
    masstodon_args['modifications'] = modifications

    # Parsing blocked_fragments:
    if 'blocked_fragments' in masstodon_args:
        masstodon_args['blocked_fragments'] = parse_blocked_fragments(masstodon_args['blocked_fragments'])

    return masstodon_args
