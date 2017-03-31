from    cvxopt          import matrix, spmatrix, sparse, spdiag, solvers

def get_incidence_matrix(frags):
    A = np.zeros((len(frags), len(frags.edges())))
    for i, N in enumerate(frags):
        for M in frags.edge[N]:
            j = frags.edge[N][M]['cnt']
            A[i,j] = 1.0
    return A


def diag(val, dim):
    return spdiag([spmatrix(val,[0],[0]) for i in xrange(dim)])


def get_G_h(var_No):
    '''Prepare for conditions Gx <= h'''
    G_mat = diag(-1.0, var_No)
    h_vec = matrix(0.0, size=(var_No, 1))
    return G_mat, h_vec


def collect_fragments(frags, Q, eps=0.0, c_penalty=0.0, solver='CVXOPT'):
    if len(frags)>1:
        for i, E in enumerate(frags.edges(data=True)):
            frags.edge[E[0]][E[1]]['cnt'] = i
        c = [] # costs
        for A, B in frags.edges():
            if A == B:
                frag_cost = Q-1-A[1]
                if A[0][0]=='c': # adds penalty for c frags
                    frag_cost += c_penalty
                c.append(frag_cost)
            else:
                c.append(Q-1-A[1]-B[1])
        c = np.array(c)
        b = np.array([ x[1]['intensity'] for x in frags.nodes(data=True) ])
        A = get_incidence_matrix(frags)
        I = len(frags)
        J = len(frags.edges())
        Amat, bmat = map(matrix,[A,b])
        cmat = matrix(c, tc='d')
        Gmat, hmat = get_G_h(J)
        hmat -= eps
        if solver=='CVXOPT':
            sol = solvers.conelp(c=cmat,G=Gmat,h=hmat,A=Amat, b=bmat)
        if solver=='scipy.linprog':
            sol = SimplexAlgorithm(c=c,A_eq=A,b_eq=b)
        return sol
    else:
        node = frags.nodes(data=True)
        return node[0][0], node[0][1]['intensity']


solvers.options['show_progress'] = False

ccs = [cc for cc in nx.connected_component_subgraphs(IDG)]
sols= [ (cc, collect_fragments(cc,Q)) for cc in ccs if len(cc)>1]
Counter([sol[1]['status'] for sol in sols])
wrong_sols = [ sol for sol in sols if sol[1]['status']!='optimal' ]

sols2 = [ (cc,collect_fragments(cc, Q, eps=100.0, c_penalty=0.0)) for cc,sol in wrong_sols]
print Counter([sol[1]['status'] for sol in sols2])
sols3 = [ (cc,collect_fragments(cc, Q, solver='scipy.linprog')) for cc,sol in wrong_sols]
[cc.nodes(data=True) for cc,sol in wrong_sols]

# sols2=[ collect_fragments(cc, Q, eps=0.0)
#     for cc,sol in zip(ccs,sols) if len(cc)>1 and sol['status'] != 'optimal' ]
# print Counter([sol['status'] for sol in sols2])
# strange = filter(lambda x: len(x[0])>1 and x[1]['status']!='optimal', zip(ccs,sols) )
# sG, sSol= strange[0]
# len(sG)
