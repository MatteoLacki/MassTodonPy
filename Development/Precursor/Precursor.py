from MassTodonPy.Data.get_dataset import get_dataset

subP = get_dataset('substanceP')
precursor = subP.precursor
precursor.get_AA(4, 'C_carbo')
precursor.get_AA(0, 'N')
precursor.get_AA(10, 'C_carbo')
precursor[10]
precursor.fasta[0]
precursor[0]
