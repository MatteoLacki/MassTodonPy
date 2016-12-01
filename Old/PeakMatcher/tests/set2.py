from collections import defaultdict
A = frozenset(['a','b','c'])

print A

for a in A:
	print a

B = defaultdict( lambda: None )
B['a']= 10

if not B['c']:
	print B['c']
else:
	print B['c']

if not 0:
	print 'B[]'
else:
	print B['c']

print False == 0

C = set('a')
print C
print 'b' in C