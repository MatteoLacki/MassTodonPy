import  argparse
import  json
from    MassTodonPy.CSI import run_masstodon, save_results



parser = argparse.ArgumentParser()
parser.add_argument( "spectrum_path",
    help="path to the spectrum file, with spectrum file extension either .txt or .mzxml, case insensitive." )
parser.add_argument( "config_path",
    help="path to the file with configuration" )
args = parser.parse_args()

with open(args.config_path, 'r') as f:
    config = json.load(f)

# results = run_masstodon(spectrum_path = args.spectrum_path, **config)
save_results(args.spectrum_path, config)
print 'Please, be kind to MassTodons.'
