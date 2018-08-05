from MassTodonPy.IsotopeCalculator.isospec_wrapper import isospec_numpy
from IsoSpecPy.IsoSpecPy import IsoSpec as isospec

import numpy as np

isospec_numpy((100,200),
                  ((1.1, 2.31), (12.4, 13.45)),
                  ((.9, .1), (.5,.5)),
                  .999)

%%timeit
isospec.getConfs()

IsoSpecify(str(mol.formula), .99)

tuple(i for i in range(10))
masses_npy = np.array(list(masses))

def getConfsNumpyBuffered(isospec):
    ffi = FFI()
    m, l, c = isospec.getConfsRaw()
    obj_cnt = len(m)
    obj_size = np.dtype('float64').itemsize
    m = np.frombuffer(ffi.buffer(m, obj_cnt * obj_size),
                      dtype = np.dtype('float64'))
    l = np.frombuffer(ffi.buffer(l, obj_cnt * obj_size),
                      dtype = np.dtype('float64'))
    return m, l

def getConfsNumpy(isospec):
    m, l, c = isospec.getConfsRaw()
    return np.array(list(m)), np.array(list(l))

%%timeit
m, l = getConfsNumpy(isospec)

%%timeit
m, l = getConfsNumpyBuffered(isospec)

%%timeit
m, l, c = isospec.getConfsRaw()

%%timeit
m, l, c = isospec.getConfs()

%%timeit
getConfsNumpy(isospec)
