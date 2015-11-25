A = {'a':2,'b':1,'c':13}

del A['b']
print A['a']

s = set(['a','b'])
print s
s.discard('a')
print s
# s.remove('c')

B = ['a','s','d']
B.pop()
print B

s.add('b')
print s

for node in set([]):
	print 'd'

if not set([]):
	print 'b'

print s
s |= set(['d'])
print s