import igraph as ig
from collections import Counter

class MissingAminoAcid(Exception):
    pass

class WrongArgument(Exception):
    pass

class OH_already_deleted(Exception):
    pass

def getAttr(seq, attr_name):
    try:
        return seq[attr_name]
    except KeyError:
        return [None] * len(seq)

def join(g1, g2, vertex_attributes, edge_attributes):
    '''Joins two graphs so that their attributes are not cancelled out.'''
    g = g1+g2
    for attr in vertex_attributes:
        g.vs[attr] = getAttr(g1.vs, attr) + getAttr(g2.vs, attr)
    for attr in edge_attributes:
        g.es[attr] = getAttr(g1.es, attr) + getAttr(g2.es, attr)
    return g

class AminoAcids(object):
    '''Makes a graph of an amino acid.'''
    def __init__(self ):
        self.aas = set(('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'))

    def addBackbone(self, A):
        B = self.Backbone()
        C = join( A, B, vertex_attributes=['elem','name'], edge_attributes=['Roep'] )
        C += ('Calpha', 'Cbeta')
        return C

    def makeAtoms(self,atomCnt):
        vAtt = {'name': [], 'elem': [] }
        for atom in atomCnt:
            vAtt['elem'].extend( [ atom ]*atomCnt[atom] )
            for aC in range(atomCnt[atom]):
                vAtt['name'].append( atom+str(aC) )
            G = ig.Graph( sum( atomCnt[a] for a in atomCnt ), vertex_attrs=vAtt)
        return G

    def Backbone(self):
        B = self.makeAtoms({ 'C':2, 'H':2, 'N':1, 'O':1 })
        B.add_edge( 'N0', 'C0', Roep='cz' )
        B.add_edge( 'C0', 'C1', Roep='ax' )
        B = B + ('H0','N0')+('C0','H1')+('C1','O0')+('C1','O0')
        B.vs['name'] = ['HNalpha', 'Halpha', 'Calpha', 'Ccarbo', 'Ocarbo', 'Nalpha']
        return B

    def A(self):
        A = self.makeAtoms({'C':1, 'H':3})+('C0','H0')+('C0','H1')+('C0','H2')
        A.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hbeta3', 'Cbeta']
        A = self.addBackbone(A)
        A['name'] = 'Alanine'
        return A

    def R(self):
        R = self.makeAtoms({ 'C':4, 'H':10, 'N':3 })+\
        ('C0','H0')+('C0','H1')+('C0','C1')+\
        ('C1','H2')+('C1','H3')+('C1','C2')+\
        ('C2','H4')+('C2','H5')+('C2','N0')+\
        ('N0','H6')+('N0','C3')+\
        ('C3','N1')+('C3','N1')+('C3','N2')+\
        ('N1','H7')+\
        ('N2','H8')+('N2','H9')
        R.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'Hdelta1', 'Hdelta2', 'HNdelta', 'HNomega1', 'HNomega2_1', 'HNomega2_2', 'Cbeta', 'Cgamma', 'Cdelta', 'guadinino-C', 'Ndelta', 'Nomega1', 'Nomega2']
        R = self.addBackbone(R)
        R['name'] ='Arginine'
        return R

    def N(self):
        N = self.makeAtoms({ 'C':2, 'H':4, 'N':1, 'O':1 })+\
        ('C0','H0')+('C0','H1')+('C0','C1')+\
        ('C1','O0')+('C1','O0')+('C1','N0')+\
        ('N0','H2')+('N0','H3')
        N.vs['name'] = ['Hbeta1', 'Hbeta2', 'HNgamma1', 'HNgamma2', 'Cbeta', 'Cgamma', 'Ogamma', 'Ngamma']
        N = self.addBackbone(N)
        N['name'] ='Asparagine'
        return N

    def D(self):
        D = self.makeAtoms({ 'C':2, 'H':3, 'O':2 })+\
        ('C0','H0')+('C0','H1')+('C0','C1')+\
        ('C1','O0')+('C1','O0')+('C1','O1')+\
        ('O1','H2')
        D.vs['name'] = ['Hbeta1', 'Hbeta2', 'HOgamma2', 'Cbeta', 'Cgamma', 'Ogamma1', 'Ogamma2']
        D = self.addBackbone(D)
        D['name'] ='Aspartic acid'
        return D

    def C(self):
        C = self.makeAtoms({ 'C':1, 'H':3, 'S':1 })+\
        ('C0','H0')+('C0','H1')+('C0','S0')+\
        ('S0','H2')
        C.vs['name'] = ['Hbeta1', 'Hbeta2', 'HSbeta', 'Cbeta', 'Sbeta']
        C = self.addBackbone(C)
        C['name'] = 'Cysteine'
        return C

    def Q(self):
        Q = self.makeAtoms({ 'C':3, 'H':6, 'N':1, 'O':1 })+\
        ('C0','H0')+('C0','H1')+('C0','C1')+\
        ('C1','H2')+('C1','H3')+('C1','C2')+\
        ('C2','O0')+('C2','O0')+('C2','N0')+\
        ('N0','H4')+('N0','H5')
        Q.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'HNdelta1', 'HNdelta2', 'Cbeta', 'Cgamma', 'Cdelta', 'Odelta', 'Ndelta']
        Q = self.addBackbone(Q)
        Q['name'] ='Glutamine'
        return Q

    def E(self):
        E = self.makeAtoms({ 'C':3, 'H':5, 'O':2 })+\
        ('C0','H0')+('C0','H1')+('C0','C1')+\
        ('C1','H2')+('C1','H3')+('C1','C2')+\
        ('C2','O0')+('C2','O0')+('C2','O1')+\
        ('O1','H4')
        E.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'HOdelta', 'Cbeta', 'Cgamma', 'Cdelta', 'Odelta1', 'Odelta2']
        E = self.addBackbone(E)
        E['name'] ='Glutamic acid'
        return E

    def G(self):
        G = self.Backbone()
        G.add_vertex('Halpha2',elem='H')
        G.add_edge('Halpha2','Calpha')
        G.vs.find('Halpha')['name'] = 'Halpha1'
        G['name'] ='Glycine'
        return G

    def H(self):
        H = self.makeAtoms({ 'C':4, 'H':5, 'N':2 })+\
        ('C0','H0')+('C0','H1')+('C0','C1')+\
        ('C1','C2')+('C1','C2')+('C1','N0')+\
        ('N0','C3')+('N0','C3')+\
        ('C3','H3')+('C3','N1')+('N1','H2')+\
        ('N1','C2')+\
        ('C2','H4')
        H.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hpi', 'Htau', 'H4', 'Cbeta', 'C4', 'C5', 'C2', 'Npi', 'Ntau']
        ['Hbeta1', 'Hbeta2', 'Htau', 'H3', 'H2', 'Cbeta', 'C4', 'C5', 'C2', 'Npi', 'Ntau']
        H = self.addBackbone(H)
        H['name'] ='Histidine'
        return H

    def I(self):
        I = self.makeAtoms({ 'C':4, 'H':9 })+\
        ('C0','H0')+('C0','C1')+('C0','C2')+\
        ('C1','H1')+('C1','H2')+('C1','H3')+\
        ('C2','H4')+('C2','H5')+('C2','C3')+\
        ('C3','H6')+('C3','H7')+('C3','H8')
        I.vs['name'] = ['Hbeta', 'Hbeta1_1', 'Hbeta1_2', 'Hbeta1_3', 'Hgamma1', 'Hgamma2', 'Hdelta1', 'Hdelta2', 'Hdelta3', 'Cbeta', 'Cbeta1', 'Cgamma', 'Cdelta']
        I = self.addBackbone(I)
        I['name'] ='Isoleucine'
        return I

    def L(self):
        L = self.makeAtoms({ 'C':4, 'H':9 })+\
        ('C0','H0')+('C0','H1')+('C0','C1')+\
        ('C1','H2')+('C1','C2')+('C1','C3')+\
        ('C2','H3')+('C2','H4')+('C2','H5')+\
        ('C3','H6')+('C3','H7')+('C3','H8')
        L.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hgamma', 'Hgamma1_1', 'Hgamma1_2', 'Hgamma1_3', 'Hgamma2_1', 'Hgamma2_2', 'Hgamma2_3', 'Cbeta', 'Cgamma', 'Cgamma1', 'Cgamma2']
        L = self.addBackbone(L)
        L['name'] ='Leucine'
        return L

    def K(self):
        K = self.makeAtoms({ 'C':4, 'H':10, 'N':1 })+\
        ('C0','H0')+('C0','H1')+('C0','C1')+\
        ('C1','H2')+('C1','H3')+('C1','C2')+\
        ('C2','H4')+('C2','H5')+('C2','C3')+\
        ('C3','H6')+('C3','H7')+('C3','N0')+\
        ('N0','H8')+('N0','H9')
        K.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'Hdelta1', 'Hdelta2', 'Hepsilon1', 'Hepsilon2', 'HNepsilon1', 'HNepsilon2', 'Cbeta', 'Cgamma', 'Cdelta', 'Cepsilon', 'Nepsilon']
        K = self.addBackbone(K)
        K['name'] ='Lysine'
        return K

    def M(self):
        M = self.makeAtoms({ 'C':3, 'H':7, 'S':1 })+\
        ('C0','H0')+('C0','H1')+('C0','C1')+\
        ('C1','H2')+('C1','H3')+('C1','S0')+\
        ('S0','C2')+\
        ('C2','H4')+('C2','H5')+('C2','H6')
        M.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'Hepsilon1', 'Hepsilon2', 'Hepsilon3', 'Cbeta', 'Cgamma', 'Cepsilon', 'Sdelta']
        M = self.addBackbone(M)
        M['name'] = 'Methionine'
        return M

    def F(self):
        F = self.makeAtoms({ 'C':7, 'H':7 })+\
        ('C0','H0')+('C0','H1')+('C0','C1')+\
        ('C1','C2')+('C1','C2')+('C1','C6')+\
        ('C2','H2')+('C2','C3')+\
        ('C3','H3')+('C3','C4')+('C3','C4')+\
        ('C4','H4')+('C4','C5')+\
        ('C5','H5')+('C5','C6')+('C5','C6')+\
        ('C6','H6')
        F.vs['name']=['Hbeta1', 'Hbeta2', 'H6', 'H5', 'H4', 'H3', 'H2', 'Cbeta', 'C1', 'C6', 'C5', 'C4', 'C3', 'C2']
        F = self.addBackbone(F)
        F['name'] = 'Phenylalanine'
        return F

    def P(self): # This is special: links to backbone in two places.
        P = self.makeAtoms({ 'C':5, 'H':7, 'N':1, 'O':1 })
        P.add_edge( 'C0', 'C1', Roep='ax' )
        P.add_edge( 'N0', 'C1', Roep='cz' )
        P = P + ('C1','H0')+\
        ('C0','O0')+('C0','O0')+\
        ('C1','C2')+\
        ('C2','H1')+('C2','H2')+('C2','C3')+\
        ('C3','H3')+('C3','H4')+('C3','C4')+\
        ('C4','N0')+('C4','H5')+('C4','H6')
        P.vs['name'] = ['Halpha', 'H3_1', 'H3_2', 'H4_1', 'H4_2', 'H5_1', 'H5_2', 'Ccarbo', 'C2', 'C3', 'C4', 'C5', 'Ocarbo', 'N1']
        P['name'] = 'Proline'
        return P

    def S(self):
        S = self.makeAtoms({ 'C':1, 'H':3, 'O':1 })+\
        ('C0','H0')+('C0','H1')+('C0','O0')+\
        ('O0','H2')
        S.vs['name']=['Hbeta1', 'Hbeta2', 'HObeta', 'Cbeta', 'Obeta']
        S = self.addBackbone(S)
        S['name'] = 'Serine'
        return S

    def T(self):
        T = self.makeAtoms({ 'C':2, 'H':5, 'O':1 })+\
        ('C0','H0')+('C0','O0')+('C0','C1')+\
        ('O0','H1')+\
        ('C1','H2')+('C1','H3')+('C1','H4')
        T.vs['name']=['Hbeta', 'HObeta', 'Hgamma1', 'Hgamma2', 'Hgamma3', 'Cbeta', 'Cgamma', 'Obeta']
        T = self.addBackbone(T)
        T['name'] = 'Threonine'
        return T

    def W(self):
        W = self.makeAtoms({ 'C':9, 'H':8, 'N':1 })+\
        ('C0','H0')+('C0','H1')+('C0','C1')+\
        ('C1','C2')+('C1','C3')+('C1','C3')+\
        ('C2','C4')+('C2','C4')+('C2','C5')+\
        ('C3','N0')+('C3','H2')+\
        ('N0','H3')+('N0','C5')+\
        ('C5','C6')+('C5','C6')+\
        ('C6','H4')+('C6','C7')+\
        ('C4','C8')+('C4','H5')+\
        ('C7','H6')+('C7','C8')+('C7','C8')+\
        ('C8','H7')
        W.vs['name']= ['Hbeta1', 'Hbeta2', 'H2', 'HN1', 'H7', 'H4', 'H6', 'H5', 'Cbeta', 'C3', 'C3a', 'C2', 'C4', 'C7a', 'C7', 'C6', 'C5', 'N1']
        W = self.addBackbone(W)
        W['name'] = 'Tryptophan'
        return W

    def Y(self):
        Y = self.makeAtoms({ 'C':7, 'H':7, 'O':1 })+\
        ('C0','H0')+('C0','H1')+('C0','C1')+\
        ('C1','C2')+('C1','C2')+('C1','C3')+\
        ('C2','H2')+('C2','C4')+\
        ('C3','H3')+('C3','C5')+('C3','C5')+\
        ('C4','H4')+('C4','C6')+('C4','C6')+\
        ('C5','H5')+('C5','C6')+\
        ('C6','O0')+\
        ('O0','H6')
        Y.vs['name']=['Hbeta1', 'Hbeta2', 'H6', 'H2', 'H5', 'H3', 'HO4', 'Cbeta', 'C1', 'C6', 'C2', 'C5', 'C3', 'C4', 'O4']
        Y = self.addBackbone(Y)
        Y['name'] = 'Tyrosine'
        return Y

    def V(self):
        V = self.makeAtoms({ 'C':3, 'H':7 })+\
        ('C0','C1')+('C0','C2')+('C0','H0')+\
        ('C1','H1')+('C1','H2')+('C1','H3')+\
        ('C2','H4')+('C2','H5')+('C2','H6')
        V.vs['name']=['Hbeta', 'Hbeta1_1', 'Hbeta1_2', 'Hbeta1_3', 'Hgamma1', 'Hgamma2', 'Hgamma3', 'Cbeta', 'Cbeta1', 'Cgamma']
        V = self.addBackbone(V)
        V['name'] = 'Valine'
        return V

    def plot(self,target=''):
        elem_color = {'C':'grey', 'H':'red', 'N':'blue', 'O':'white', 'S':'pink' }
        visual_style = {}
        visual_style['vertex_color']= [ elem_color[gender] for gender in self.G.vs['elem']]
        visual_style['vertex_label']= self.G.vs['name']
        visual_style['edge_label']  = self.G.es['Roep']
        if len(target) == 0:
            ig.plot( self.G, **visual_style )
        else:
            ig.plot( self.G, target=target, **visual_style )

    def add_OH(self):
        self.G.vs.find('Ocarbo')['name'] = 'Ocarbo1'
        self.G.add_vertex(name='Ocarbo2', elem='O')
        self.G.add_vertex(name='HOcarbo2', elem='H')
        self.G.add_edge('Ocarbo2','HOcarbo2')
        self.G.add_edge('Ccarbo','Ocarbo2',Roep='by')

    def add_H(self):
        if self.G['name'] == 'Proline':
            self.G.add_vertex(name='HN1', elem='H')
            self.G.add_edge('N1','HN1')
        else:
            self.G.vs.find('HNalpha')['name'] = 'HNalpha1'
            self.G.add_vertex(name='HNalpha2', elem='H')
            self.G.add_edge('Nalpha','HNalpha2') # Might this be another by bond? Assume not

    def get(self):
        aminoAcids = {}
        for aa in self.aas:
            G = getattr(self, aa)()
            C = G.vs.find('Ccarbo').index
            if aa == 'P':
                N = G.vs.find('N1').index
            else:
                N = G.vs.find('Nalpha').index

            aminoAcids[aa] = { 'graph' : G, 'NalphaIDX' : N, 'CcarboIDX' : C }
        return aminoAcids

    def are_encoded(self, fastas):
        return all( set(f).issubset(self.aas) for f in fastas )
