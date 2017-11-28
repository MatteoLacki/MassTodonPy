%load_ext autoreload
%autoreload 2

from MassTodonPy.Formula.Formula import Formula

f = Formula('C100H202')
f.pattern
f.recompile_pattern('([A-Z][a-z]?)([4-9]*)')
f.pattern

g = Formula('C100H202')
g

list(Formula('C100H202').items())

repr(Formula('C100H202'))
str(Formula('C100H202'))
str(Formula('H202C100'))


f = Formula('C100H202').copy()
f['H'] = -10
f.check_positivity()
