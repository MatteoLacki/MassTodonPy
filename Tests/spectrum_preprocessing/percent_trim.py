import numpy as np
from operator import itemgetter


if not None and not None:
    print 'lb'
else:
    print ';'

mz = np.linspace(1,10,10)
intensities = np.random.exponential(scale=1.0, size=mz.shape)


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
