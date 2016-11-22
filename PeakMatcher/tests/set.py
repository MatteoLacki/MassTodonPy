A = set(['a','b',12])

print A
A.pop()
print A

D = {}

D['Matteo'] = ':)'
D['Yani'] = '8D'

print D

L = [ 'a' ]

print L
S = L
S.append('b')

print S
print L

CL = L[:]
print CL
CL.append('s')
print CL
print L

s = frozenset(['a','b','b'])
print s
s2=frozenset(['b','a'])
print s2
D[s] = 5