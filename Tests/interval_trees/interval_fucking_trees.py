from intervaltree import Interval as II, IntervalTree as Itree


T = Itree()
T[1:2] = (1.5,10.)
T[1.5:2.5] = (2.0,12.)

T1 = [(1.5,10.), (2.0,12.), (2.1,13.)]

M = [ a for a,b in T1 ]
I = [ b for a,b in T1 ]

spec1 = (M[:-1],I[:-1])
spec2 = (M,I)

mzPrec = .4
T = Itree( II(mz-mzPrec, mz+mzPrec, (mz, intensity)) for mz, intensity in zip(*spec1) )

S = Itree( II(mz-mzPrec, mz+mzPrec, (mz, intensity)) for mz, intensity in zip(*spec2) )
S

S = Itree()
S[1:2] = (1.5,10.)
S[1.5:2.5] = (2.0,12.)
S[1.6:2.6] = (2.1,13.)


T
S
