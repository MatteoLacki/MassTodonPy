from multiprocessing import Pool

def f(x):
    return x+1

def g(x):
	return x[0]+x[1]

def h(x,y):
	return x+y

def iterator(k):
	for i in range(k):
		yield (i,i)



pool= Pool(processes=2)
res = pool.map(f, range(10))
res = pool.map(g, iterator(10))
print res


pool.close(); pool.join()



