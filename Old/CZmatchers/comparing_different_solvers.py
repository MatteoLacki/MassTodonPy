from    graph_compute   import max_cost_flaw

def etnod_ptr_on_missing_cofragment(nQ, nG, LogProbEtnod, LogProbPtr, Q):
    '''Get the number of ETnoD and PTR reactions on an edges with minimal cost.'''
    if LogProbEtnod > LogProbPtr:
        Netnod  = Q - 1 - nQ
        Nptr    = 0
    else:
        Netnod  = nG
        Nptr    = Q - 1 - nQ - nG
    return Netnod, Nptr


def etnod_ptr_on_c_z_pairing( q0, g0, q1, g1, Q ):
    '''Get the number of ETnoD and PTR reactions on a regular edge.'''
    Netnod  = g0 + g1
    Nptr    = Q - 1 - g0 - g1 - q0 - q1
    return Netnod, Nptr


def get_break_point( nType ):
    '''Get the amino acid number that was cleft.'''
    if nType[0] == 'c':
        bP = int(nType[1:])
    else:
        bP = L - int(nType[1:])
    return bP


def logBinomial(m,n):
    return lgamma(m+n+1.0)-lgamma(m+1.0)-lgamma(n+1.0)


def get_weight(C, Z, LogProb, Q):
    '''Weight for the weighted max flow optimization problem.'''
    (cT, cQ, cG), (zT, zQ, zG) = C, Z
    Netnod, Nptr= etnod_ptr_on_c_z_pairing( cQ, cG, zQ, zG, Q )
    w_e = logBinomial(Netnod, Nptr)

    logPptr     = LogProb['PTR']
    logPetnod   = LogProb['ETnoD']

    bP = get_break_point(cT)
    Cetnod, Cptr = etnod_ptr_on_missing_cofragment(cQ, cG, logPetnod, logPptr, Q)
    Zetnod, Zptr = etnod_ptr_on_missing_cofragment(zQ, zG, logPetnod, logPptr, Q)

    if logPetnod > logPptr:
        W_edge  = (logPptr-logPetnod) * Nptr - (Q-1)*logPetnod + w_e - LogProb[bP]
    else:
        w_cc    = logBinomial(Cetnod, Cptr)
        w_zz    = logBinomial(Zetnod, Zptr)
        W_edge  = -logPptr*(Q-1) + w_e - w_cc - w_zz - LogProb[bP]
    return W_edge


def initialize_flow_graph(G, Q, LogProb):
    '''Construct the flow graph corresponding to one pairing problem.'''

    totalIntensity = sum(G.node[N]['intensity'] for N in G )
    FG = nx.DiGraph()
    FG.add_node('S', demand= -totalIntensity) # start
    FG.add_node('T', demand=  totalIntensity) # terminus/sink
    for C in G:
        if C[0][0]=='c':
            FG.add_node(C)
            FG.add_edge( 'S', C, capacity=G.node[C]['intensity'] )
            for Z in G[C]:
                FG.add_node(Z)
                FG.add_edge( Z, 'T', capacity = G.node[Z]['intensity'] )
                FG.add_edge( C,  Z,   weight  = get_weight(C, Z, LogProb, Q) )
    max_cost_flaw(FG, 'S', 'T', cost="weight", capacity="capacity")
    return FG

# G.nodes(data=True)
# G.edges(data=True)

def initialize_flow_graph_nx(G, Q, LogProb):
    '''Construct the flow graph corresponding to one pairing problem.'''

    totalIntensity = sum(G.node[N]['intensity'] for N in G )
    FG = nx.DiGraph()
    FG.add_node('S', demand= -totalIntensity) # start
    FG.add_node('T', demand=  totalIntensity) # terminus/sink
    FG.add_edge('S','T')
    for C in G:
        if C[0][0]=='c':
            FG.add_node(C)
            FG.add_edge( 'S', C, capacity=G.node[C]['intensity'] )
            for Z in G[C]:
                FG.add_node(Z)
                FG.add_edge(Z, 'T', capacity = G.node[Z]['intensity'])
                weight = int(-get_weight(C,Z,LogProb, Q) * 10000)
                print weight
                FG.add_edge(C, Z, weight = weight)

    print FG.nodes(data=True)
    print FG.edges(data=True)

    minFlaw = nx.min_cost_flow(FG)
    return minFlaw
