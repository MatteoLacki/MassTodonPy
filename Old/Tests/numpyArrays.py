a = np.arange(6).reshape(2,3)
for x in np.nditer(a, order='F'):
    print x,

a

for x in np.nditer(a, order='C'):
    print x,

a

a = np.arange(6).reshape(2,3)
a
for x in np.nditer(a, flags=['external_loop']):
    print x,
