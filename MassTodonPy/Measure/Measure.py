from __future__ import absolute_import, division, print_function
from bisect import bisect_left
from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool
import csv
import numpy as np
from operator import itemgetter
from six.moves import range, zip

from MassTodonPy.Data.Constants import infinity
from MassTodonPy.Misc.strings import repr_long_list
# from MassTodonPy.Misc.sorting import sort_by_first
from MassTodonPy.Parsers.Paths import parse_path
from MassTodonPy.Reporter.buffers import buffers


class Measure(object):
    """Store a discrete finite measure with atoms on the real line."""

    def __init__(self, atoms=np.array([]),
                       masses=np.array([]),
                       sort=True):
        """Initialize a measure.

        Parameters
        ==========
        atoms : numpy array
            The atoms upon which the measure holds the mass.
        masses : numpy array
            The masses on atoms.

        Remarks
        =======
        If you have sorted the measure beforehand, do set 'is_sorted' to True.

        """
        self.atoms = atoms
        self.masses = masses
        if sort:
            self.sort()
        self._store_names = ('atom', 'mass')

    def sort(self):
        """Sort measure by atomic values."""       
        self.atoms, self.masses = (np.array(x) for x in \
                                   sort_by_first(self.atoms, self.masses))


    def __has_type_of(self, other):
        """Assert that 'self' and 'other' have the same type."""
        assert self.__class__.__name__ == other.__class__.__name__,\
             "\tIllegal to add class {0} to class {1}.\n".format(
                other.__class__.__name__,\
                self.__class__.__name__)

    def __add__(self, other):
        """Add two measures.

        Parameters
        ----------
        other : Measure
            A measure we want to stack on top of this one.

        """
        atoms = np.concatenate((self.atoms, other.atoms))
        masses = np.concatenate((self.masses, other.masses))
        self.__has_type_of(other)
        new_measure = self.__class__(atoms, masses)
        new_measure.__aggregate()
        return new_measure

    def __radd__(self, other):
        """Add two measures.

        Parameters
        ----------
        other : Measure
            A measure we want to stack on top of this one.

        """
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __iadd__(self, other):
        """Add two measures.

        Parameters
        ----------
        other : Measure
            A measure we want to stack on top of this one.

        """
        self.__has_type_of(other)
        self.atoms = np.concatenate((self.atoms, other.atoms))
        self.masses = np.concatenate((self.masses, other.masses))
        self.__aggregate()
        return self

    def copy(self):
        """Make a deep copy of me."""
        out = self.__class__(self.atoms, self.masses)
        return out

    def __mul__(self, scalar):
        """Multiply by a scalar."""
        if scalar == 0:
            return self.__class__()
        elif scalar == 1:
            return self.copy()
        else:
            return self.__class__(self.atoms, scalar * self.masses)

    def __rmul__(self, scalar):
        """Multiply by a scalar."""
        if scalar == 1:
            return self
        else:
            return self.__mul__(scalar)

    def __imul__(self, scalar):
        """Multiply by a scalar."""
        if scalar != 1:
            self.masses = self.masses * scalar
        return self

    def __aggregate(self):
        """Aggregate masses with the same atoms."""
        self.atoms, indices = np.unique(self.atoms, return_inverse=True)
        self.masses = np.bincount(indices, weights=self.masses)

    def round_atoms(self, precision=infinity):
        """Round the atoms of the measure to a given precision.

        Parameters
        ----------
        precision : integer
            The number of digits after which the atoms' masses get rounded.
            E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
            Defaults to 'inf', which prevents any rounding.

        """
        if precision != infinity:
            self.atoms = np.around(self.atoms, precision)
            self.__aggregate()

    def trim(self, cut_off):
        """Trim masses below the provided cut off.

        Parameters
        ----------
        cut_off : float

        """
        if cut_off > 0:
            self.atoms = self.atoms[self.masses >= cut_off]
            self.masses = self.masses[self.masses >= cut_off]

    def split_measure(self, cut_off):
        """Split measure into two according to the cut off on masses.

        Retain the measure with masses greater or equal to the cut off.
        Parameters
        ----------
        cut_off : float
        Returns
        ----------
        other : Measure
            A measure with masses strictly below the cut off.
        """
        other = self.__class__(self.atoms[self.masses < cut_off],
                               self.masses[self.masses < cut_off])
        self.trim(cut_off)
        return other

    def get_P_set_cut_off(self, P=.99):
        """Get the cut off resulting in optimal P-set.

        Parameters
        ----------
        P : float
            The percentage of the initial total value that the new measure
            will contain. The new measure contains only atoms with
            heighest masses.

        """
        # TODO replace this with something brighter, as linearly operating heap.
        assert 0.0 <= P <= 1.0, "Wrong P for P-optimal set."
        total_value = self.masses.sum()
        i = 0
        S = 0.0
        masses = np.sort(self.masses)[::-1]
        for intensity in masses:
            S += intensity
            if S < P:
                break
        return intensity

    def __repr__(self):
        """Represent the measure."""
        if len(self.atoms) > 0:
            out = "{0}:\n\t{1} = {2}\n\t{3} = {4}\n".format(self.__class__.__name__,
                                                            self._store_names[0],
                                                            repr_long_list(self.atoms)[1:-1],
                                                            self._store_names[1],
                                                            repr_long_list(self.masses)[1:-1])
        else:
            out = self.__class__.__name__ + " (empty)"
        return out

    def __len__(self):
        """Get size of the measure: the number of atoms."""
        return len(self.atoms)

    def __iter__(self):
        """Iterate over pairs (atom, mass)."""
        return zip(self.atoms, self.masses)

    def __getitem__(self, key):
        """Filter atoms between 'L' and 'R'.

        Parameters
        ==========
        key : tuple
            Either (left_mz, right_mz) or
            (left_mz, right_mz, provide_idx),
            where 'provide_idx' is boolean.
        Returns
        =======
        out : generator
            Generate tuples '(mz, intensity)'
            or '(idx, mz, intensity)',
            where 'idx' is the unique ID of the atom.

        """
        try:
            if len(key) == 2:
                L, R = key
                provide_idx = False
            elif len(key) == 3:
                L, R, provide_idx = key
            idx = bisect_left(self.atoms, L)
            while self.atoms[idx] <= R:
                if not provide_idx:
                    yield self.atoms[idx], self.masses[idx]
                else:
                    yield idx, self.atoms[idx], self.masses[idx]
                idx += 1
        except IndexError:
            return

    def total_mass(self):
        return self.masses.sum()

    def write(self, path):
        """Write the spectrum to a csv or tsv file.

        Parameters
        ==========
        path : str
            A path to the file to write to.
        """
        file_path, file_name, file_ext = parse_path(path)
        delimiter = ',' if file_ext == '.csv' else '\t'
        with open(path, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=delimiter)
            writer.writerow(self._store_names)
            for atom, mass in self:
                writer.writerow([atom, mass])

    def plot(self,
             path="",
             mode="inline",
             bar_width=.01,
             width=800,
             height=600,
             show_plot=True,
             _simple=False,
             _max_buffer_len=2,
             **kwds):
        """Make an interactive Bokeh barplot.

        Parameters
        ==========
        path : string
            Path to where to save the output html file.
            If not provided, MassTodon will use the folder it is called from.
        mode : string
            The mode of plotting a bokeh plot, e.g. 'inline'.
        bar_width : float
            The width of one peak.
        width : int
            The width of the plot.
        height : int
            The height of the plot.
        show_plot : bool
            Show the plot in the default web browser.
        _simple : bool
            Show simpler version of the plot.
        _max_buffer_len : float
            The maximal length of invisible buffer space between peaks.
            These extend the area that triggers the hover tool on.
        """
        output_file(path, mode=mode)

        # this could be better handled in Python
        bar_width = float(bar_width)
        width = int(width)
        height = int(height)
        show_plot = bool(show_plot)
        _simple = bool(_simple)
        _max_buffer_len = int(_max_buffer_len)

        if len(self.atoms) > 0:
            if _simple:  # a simple plot
                if bar_width == 1:  # get minimal width - the minimal space between atoms
                    prev_atom = self.atoms[0]
                    bw = infinity
                    for i in range(1, len(self.atoms)):
                        atom = self.atoms[i]
                        bw = min(atom - prev_atom, bw)
                        prev_atom = atom
                    bar_width = min(bw, bar_width)  # prevent infinite width
                hover = HoverTool(tooltips=[(self._store_names[0], "@x{0,0.000}"),
                                            (self._store_names[1], "@top{0,0}")],
                                  mode='vline')

            TOOLS = "crosshair pan wheel_zoom box_zoom undo redo reset box_select lasso_select save".split(' ')

            max_mass = self.masses.max() * 1.05
            plot = figure(tools=TOOLS,
                          width=width,
                          height=height,
                          y_range=(0, max_mass))

            raw_measure = plot.vbar(x=self.atoms,
                                    top=self.masses,
                                    width=bar_width,
                                    color='black')
            plot.sizing_mode = 'scale_both'

            if not _simple:  # nice, complex plot
                buffers_L, buffers_R = buffers(self.atoms,
                                               self.atoms,
                                               _max_buffer_len)
                source = ColumnDataSource(data={'atoms': self.atoms,
                                                'top': self.masses,
                                                'left': buffers_L,
                                                'right': buffers_R})
                invisible_buffers = plot.quad(left='left',
                                              right='right',
                                              top='top',
                                              bottom=0.0,
                                              alpha=0.0,
                                              color='black',
                                              source=source)
                hover = HoverTool(renderers=[invisible_buffers],
                                  tooltips=[(self._store_names[1], "@top{0,0}"),
                                            (self._store_names[0], "@atoms{0,0.000}")],
                                  mode='vline')
            plot.add_tools(hover)
            if show_plot:
                show(plot)
            return plot
        else:
            print('You try to plot emptiness: look deeper into your heart.')
