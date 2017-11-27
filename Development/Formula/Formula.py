%load_ext autoreload
%autoreload 2

from MassTodonPy.Formula.Formula import Formula

list(Formula('C100H202').items())

repr(Formula('C100H202'))
str(Formula('C100H202'))
str(Formula('H202C100'))


f = Formula('C100H202').copy()
f
