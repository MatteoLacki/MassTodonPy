def _optimize(self, G):
    """Match the intensities in a cluster.

    Decorates the self.graph with flow values.
    """
    Q = self.precursor.q
    # lavish: all fragments lose cofragments
    lavish = sum((Q - 1 - N.q) * I for N, I in
                 G.nodes.data('intensity'))
    self.I_lavish += lavish
    if len(G) > 1:  # not trivial
        FG = nx.DiGraph()
        FG.add_node('S')  # start
        FG.add_node('T')  # terminus/sink
        for C, C_intensity in G.nodes.data("intensity"):
            if C.type is 'c':
                FG.add_node(C)
                FG.add_edge('S', C, capacity=C_intensity)
                for Z in G[C]:
                    Z_intensity = G.node[Z]['intensity']
                    FG.add_node(Z)
                    FG.add_edge(C, Z)
                    FG.add_edge(Z, 'T', capacity=Z_intensity)
        max_flow_min_cost = nx.max_flow_min_cost(FG, 'S', 'T')
        total_flow, flows = nx.maximum_flow(FG, 'S', 'T')
        self.I_ETnoD_PTR_fragments += lavish - (Q - 1) * total_flow
        for N in G:
            self.graph.add_edge(N, N)
        # no double count: flows = { start: {end: {value}, .. }, .. }
        for N in flows:
            for M in flows[N]:
                if N is 'S':  # M is a C fragment
                    self.graph[M][M]['flow'] = G.node[M]['intensity'] - flows[N][M]
                elif M is 'T':  # N is a Z fragment
                    self.graph[N][N]['flow'] = G.node[N]['intensity'] - flows[N][M]
                else:  # N is a C and M a Z fragment
                    self.graph[N][M]['flow'] = flows[N][M]
    else:  # trivial
        self.I_ETnoD_PTR_fragments += lavish
        N, N_intensity = list(G.nodes.data('intensity'))[0]
        self.graph.add_edge(N, N, flow=N_intensity)




#            # Old QP code
        #     TotalFrags = sum(I)
        #     TotalPTR   = 0.0
        #     TotalETnoD = 0.0
        #     for i, (N0, N1) in enumerate(G.edges_iter()):
        #         if N0 != N1:
        #             TotalETnoD  += I[i]*G.edge[N0][N1]['ETnoD']
        #             TotalPTR    += I[i]*G.edge[N0][N1]['PTR']
        # else:
        #     (nType, nQ, nG), Data = G.nodes(data=True)[0]
        #     I   = Data['intensity']
        #     TotalPTR   = 0
        #     TotalETnoD = 0
        #     TotalFrags = I
        #     S = {'status':'trivial', 'x':I}
