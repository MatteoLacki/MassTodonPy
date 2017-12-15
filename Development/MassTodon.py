%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MassTodon import MassTodon

mol = get_dataset('substanceP')
precursor = {'name': 'substanceP',
             'fasta': mol.precursor.fasta,
             'charge': 3}

%%time
masstodon = MassTodon(spectrum=mol.spectrum,
                      precursor=precursor,
                      mz_precision=.05,
                      _devel=True)

# masstodon.cz_match.intensities
# masstodon.cz_match.branching_ratio
# masstodon.cz_match.probabilities
