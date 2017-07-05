import numpy as np
from operator import itemgetter
from itertools import izip as zip

if not None and not None:
    print 'lb'
else:
    print ';'

mz = np.linspace(1,10,100)
intensities = np.random.exponential(scale=1.0, size=mz.shape)



def aggregate( keys, values=None ):
    '''Aggregate values with the same keys.'''
    uniqueKeys, indices = np.unique( keys, return_inverse=True)
    return uniqueKeys, np.bincount( indices, weights=values )

def round_spec(mz, intensity, digits=2):
    '''Aggregate the spectrum so that intensities of masses with the same number of significant digits are summed.'''
    mz = np.round(mz, digits)
    mz, intensity = aggregate(mz, intensity)
    return mz, intensity

spectrum = round_spec(mz, intensities, digits=1)

a = zip(*spectrum)
list(a)

"test %s" % 10

var = '%(foo)s %(foo)s %(foo)s' % { 'foo': 'look_at_me_three_times' }
var = '%(bla)s %(get)s %(foo)s' % { 'foo': 'look_at_me_three_times', 'bla': 'haha', 'get': 'bu' }
var

def quantile_trim0(mz, intensities, perc = .95):
    mz_res, intensities_res = tuple( np.array(x) for x in zip(*sorted(zip(mz, intensities), key=itemgetter(1))) )
    i = 0
    S = 0.0
    total = intensities_res.sum()
    while True:
        S += intensities_res[i]/total
        if S < 1.0-perc:
            i += 1
        else:
            break

    mz_res = mz_res[ intensities_res >= intensities_res[i] ]
    intensities_res = intensities_res[ intensities_res >= intensities_res[i] ]
    return tuple( np.array(x) for x in zip(*sorted(zip(mz_res, intensities_res), key=itemgetter(0))) )

mz, intensities = quantile_trim0(mz, intensities, perc = .95)
