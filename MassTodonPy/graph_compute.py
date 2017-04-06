import networkx as nx



def fill_in_attribs(BFG, attrname, defaultval):
    for n in BFG.node:
        for n2 in BFG[n].keys():
            if not attrname in BFG[n][n2]:
                BFG[n][n2][attrname] = defaultval


def get_empty_residual_graph(BFG, capacity="capacity", cost="cost"):
    BFRG = nx.DiGraph()
    BFRG.add_nodes_from(BFG.node)
    for v1, v2 in BFG.edges():
        if BFG.has_edge(v2, v1):
            raise NotImplemented("Case not implemented.")
        BFRG.add_edges_from([(v1, v2, {capacity:BFG[v1][v2][capacity], cost:BFG[v1][v2][cost]}),
                             (v2, v1, {capacity: 0, cost:-BFG[v1][v2][cost]})])
    return BFRG


def get_constrained_graph(BFRG, capacity="capacity"):
    BFCG = nx.DiGraph()
    BFCG.add_nodes_from(BFRG)
    for v1, v2 in BFRG.edges():
        if BFRG[v1][v2][capacity] > 0:
            BFCG.add_edge(v1, v2)
    return BFCG


def get_path_params(BFAG, path, parname):
    ret = []
    for i in xrange(len(path)-1):
        ret.append(BFAG[path[i]][path[i+1]][parname])
    return ret

    
def add_flaw_to_residual_graph(BFRG, path, flow, capacity="capacity"):
    for i in xrange(len(path)-1):
        BFRG[path[i]][path[i+1]][capacity] -= flow
        BFRG[path[i+1]][path[i]][capacity] += flow


def select_augmenting_path(BFRG, cost="cost", capacity="capacity"):
    cycles = list(nx.simple_cycles(get_constrained_graph(BFRG)))
    for cycle in cycles:
        cycle.append(cycle[0])
    cycles = cycles + [x[::-1] for x in cycles]
    augmenting_cycles = []
    for cycle in cycles:
        delta_cost = sum(get_path_params(BFRG, cycle, cost))
        if not delta_cost > 0:
            continue
        flow = min(get_path_params(BFRG, cycle, capacity))
        if not flow > 0:
            continue
        delta_cost *= flow
        augmenting_cycles.append((delta_cost, cycle))
    augmenting_cycles.sort(reverse=True)
    if len(augmenting_cycles) == 0:
        return None
    return augmenting_cycles[0]

    



def max_cost_flaw(BFG, src, sink, cost="cost", capacity="capacity"):
    assert(src in BFG)
    assert(sink in BFG)

    fill_in_attribs(BFG, cost, 0)
    fill_in_attribs(BFG, capacity, float("+inf"))

    BFG.add_edges_from([(sink, src, {cost:0, capacity:float("+inf")})])

    BFRG = get_empty_residual_graph(BFG, capacity=capacity, cost=cost)

    while True:
        ap = select_augmenting_path(BFRG, capacity=capacity, cost=cost)
        if ap is None:
            break
        flow = min(get_path_params(BFRG, ap[1], capacity))
        add_flaw_to_residual_graph(BFRG, ap[1], flow)

    BFG.remove_edge(sink, src)

    for v1, v2 in BFG.edges():
        BFG[v1][v2]['flaw'] = BFRG[v2][v1][capacity]







BFG = nx.DiGraph()
BFG.add_nodes_from(range(6))
BFG.add_edges_from([
                    (0, 2, {'capacity':4}), 
                    (0, 4, {'capacity':5}), 
                    (4, 3, {'cost':10}), 
                    (2, 3, {'cost':5}), 
                    (4, 5, {'cost':6}), 
                    (3, 1, {'capacity':5}), 
                    (5, 1, {'capacity':4})
                   ])




max_cost_flaw(BFG, 0, 1)

for v1, v2 in BFG.edges():
    print v1, v2, BFG[v1][v2]['flaw']
