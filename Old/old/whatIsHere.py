try:
   import cPickle as pickle
except:
   import pickle

AA = pickle.load(open('AA.txt', 'rb'))

print(AA)
print(AA['A'])



