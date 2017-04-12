
def get_number_of_explanations(P):
    return Counter( P.node[N]['q'] for N in P if P.node[N]['type']=='M' )

def get_average_charge(P):
    return Counter( P.node[N]['q'] for N in P if P.node[N]['type']=='M' )

from scipy import stats


def average_degree_of_Gs(P):
    degrees = [ len(P[G]) for G in P if P.node[G]['type']=='G' ]
    return np.array(degrees).mean()

def average_degree_of_Gs_weighted_with_intensity(P):
    T_intensity = 0.0
    degree = 0.0
    for G in P:
        if P.node[G]['type']=='G':
            intensity = P.node[G]['intensity']
            T_intensity += intensity
            degree += len(P[G])*intensity
    return degree/T_intensity

Y = np.array([ average_degree_of_Gs(P) for s,P in zip(sols,problems) if s['status']=='unknown' ])
X = np.array([ average_degree_of_Gs_weighted_with_intensity(P) for s,P in zip(sols,problems) if s['status']=='unknown' ])
XX= np.array([ average_degree_of_Gs_weighted_with_intensity(P) for s,P in zip(sols,problems) if s['status']=='optimal' ])

X.mean()
XX.mean()

X
XX
  
stats.mode(X)[0]
stats.mode(XX)[0]
X.max()
XX.max()
# ERGO: it might be that we are trying to explain the results with too many peaks.


def mode_charge(P):
    QQ = np.array([P.node[N]['q'] for N in P if P.node[N]['type']=='M'])
    return stats.mode(QQ)[0]

x = np.array([ mode_charge(P) for s,P in zip(sols,problems) if s['status']=='unknown' ]).mean()
xx= np.array([ mode_charge(P) for s,P in zip(sols,problems) if s['status']=='optimal' ]).mean()
x
xx

def average_inverse_charge(W):
    t_cnt = 0
    t_Q = 0
    for q,cnt in W.items():
        t_cnt += cnt
        t_Q += cnt/float(q)
    return t_Q/float(t_cnt)

np.array(map(average_inverse_charge,X)).mean()
np.array(map(average_inverse_charge,XX)).mean()

XX = [ get_number_of_explanations(P) for s,P in zip(sols,problems) if s['status']=='optimal' ]
