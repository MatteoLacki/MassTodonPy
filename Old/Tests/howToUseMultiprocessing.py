from 	IsoSpecPy import IsoSpecPy
import 	cPickle as pickle
import 	time
from 	math import exp
import 	multiprocessing as mp
import  gc

# PATH = '/Volumes/doom/Users/matteo/Dropbox/Science/MassSpectrometry/IsoSpec/IsoSpecVSaltri/'
PATH = '/home/matteo/IsoSpecVSaltri/'

Data 			= pickle.load( open(PATH+'data/proteins.txt', 'rb') )
masses, probs 	= pickle.load( open(PATH+'data/isotopes.txt', 'rb') )
thresholds      = pickle.load( open( PATH+'data/proteins_thresh.txt', 'rb') )
FullData = {}
for tag, elements, counts in Data:
    FullData[tag] = [tag, elements, counts]

for tag, thr in thresholds:
    FullData[tag].append(thr)

FullData = [ FullData[tag] for tag in FullData ]

no_cores= 64
timeout = 20
tasks   = mp.Queue()
returns = mp.Queue()

for d in FullData[0:9]:
	tasks.put(d)

def subworker(q):
    ID, elements, counts, thr = q.get()
    M = [ masses[el] for el in elements]
    P = [ probs[el]  for el in elements]
    try:
        gc.collect()
        gc.disable()
        start   = time.time()
        res     = IsoSpecPy.IsoSpec( counts, M, P, thr, method = 'threshold_absolute' )
        end     = time.time()
        timeResult      = end-start
        confsNo         = len(res.getConfs())
        totalProbability= sum( exp(x[1]) for x in res.getConfs() )
        ret = (ID, timeResult, confsNo, totalProbability )
    except Exception as e:
        ret = (ID, e)
    q.put(ret)

def worker():
    while True:
        task= tasks.get()
        if task == "end":
            return
        else:
            ID  = task[0]
        q = mp.Queue()
        q.put(task)
        p = mp.Process( target = subworker, args = [q] )
        p.start()
        p.join( timeout )
        if p.is_alive():
            p.terminate()
            returns.put( (ID, 'timeout') )
        else:
            returns.put( q.get() )

procs = [ mp.Process( target=worker ) for _ in xrange(no_cores) ]
for proc in procs:
    proc.start()

for x in xrange(no_cores):
    tasks.put('end')

for proc in procs:
    proc.join()

results = []
while not returns.empty():
	results.append(returns.get())

with open(PATH + 'data/proteinsTrim.txt','w') as f:
    pickle.dump( results, f)
