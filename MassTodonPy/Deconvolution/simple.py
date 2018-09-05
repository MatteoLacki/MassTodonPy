from    networkx.linalg.attrmatrix  import attr_matrix
import  numpy                       as     np
import  matplotlib.pyplot           as     plt

from    MassTodonPy.models.nnls     import nnls


class DeconvolutionProblem(object):
    def fit(self, 
            connected_component_of_G,
            total_intensities,
            min_mz,
            max_mz,
            mean_mz,
            include_zero_intensities = False):
        self.cc     = connected_component_of_G
        mol_columns = np.array([N < 0  for N in self.cc])
        peak_rows   = np.array([N >= 0 for N in self.cc])
        X, ordering = attr_matrix(self.cc, edge_attr='prob')
        X           = X[:,mol_columns][peak_rows,:]
        ordering    = np.array(ordering)
        self.idx    = ordering[ordering >= 0]
        self.total_intensities = total_intensities[self.idx]
        self.mz_s   = min_mz[self.idx]
        self.mz_e   = max_mz[self.idx]
        self.mean_mz= mean_mz[self.idx]
        Y           = self.total_intensities
        if include_zero_intensities:
            Y = np.concatenate((Y, np.zeros(X.shape[1])))
            x = 1.0 - np.array(X.sum(axis=0)).flatten()
            X = np.concatenate((X, np.diag(x)))
        self.model  = nnls(X, Y)

    def l1(self):
        return sum(np.abs(self.model.res()))

    def l2(self):
        """This is copying what I once complained to Ludi about :)"""
        return self.model.l2()

    def spectrum(self):
        pred = self.model.fitted()
        Y    = self.model.Y
        return self.mz_s, self.mz_e, self.mean_mz, self.pred, Y[Y > 0]

    def plot(self, plt_style = 'fast',
                   bar_color = 'grey',
                   bar_alpha = 0.5,
                   show      = True):
        plt.style.use(plt_style)
        plt.bar(x       = self.mz_s,
                height  = self.model.fitted(),
                width   = self.mz_e - self.mz_s,
                align   = 'edge',
                alpha   = bar_alpha,
                color   = bar_color)
        Y = self.model.Y
        plt.scatter(self.mean_mz, Y[Y > 0], c = 'red', s = 8)
        if show:
            plt.show()


# def deconvolution_problem(connected_component_of_G,
#                           total_intensities,
#                           min_mz,
#                             max_mz,
#                             mean_mz,
#                             include_zero_intensities = False):
#     dp = DeconvolutionProblem()
#     dp.fit()

