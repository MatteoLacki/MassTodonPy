from six.moves import reduce

args = [(1,2), (3,4), (5,6)]

def foo(x, y):
    return x[0]+y[0], x[1]-y[1]

reduce(foo, args)

reduce(lambda x, y: (x[0]+y[0], x[1]-y[1]), args)
