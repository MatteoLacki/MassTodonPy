# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
#   MassTodon is free software: you can redistribute it and/or modify
#   it under the terms of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3.
#
#   MassTodon is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3 along with MassTodon.  If not, see
#   <https://www.gnu.org/licenses/agpl-3.0.en.html>.

from IsoSpecPy      import IsoSpecPy
from formulaParser  import formulaParser
from math           import exp, floor, fsum
from collections    import Counter, defaultdict
try:
  import cPickle as pickle
except:
  import pickle
from    numpy.random   import multinomial, normal
import  scipy.stats  as ss
import  numpy        as np
import  pkg_resources
from    time import time

def cdata2numpyarray(x):
    """Turn c-data into a numpy array.
    Parameters
    ----------
    x : cdata table
        A table of cdata from cffi.

    Returns
    -------
    res : array
        A numpy array of numbers.
    """
    res = np.empty(len(x))
    for i in xrange(len(x)):
        res[i] = x[i]
    return res


def agg_spec_proper(masses, probs, digits=2):
    """Aggregate values with the same keys."""
    lists = defaultdict(list)
    for mass, prob in zip(masses.round(digits), probs):
        lists[mass].append(prob)
    newMasses = np.array(lists.keys())
    newProbs  = np.empty(len(newMasses))
    for prob, mass in zip(np.nditer(newProbs,op_flags=['readwrite']), newMasses):
        prob[...] = fsum(lists[mass])
    return newMasses, newProbs


def aggregate( keys, values=None ):
    """Aggregate values with the same keys.

    Parameters
    ----------
    keys : array
        Keys, usually m/z values.
    values : array
        Values to aggregate, usually intensities.

    Returns
    -------
    out : tuple
        A tuple containing unique keys and aggregated values.
    """
    uniqueKeys, indices = np.unique( keys, return_inverse=True)
    return uniqueKeys, np.bincount( indices, weights=values )


def merge_runs(spec1, spec2):
    mz = np.concatenate((spec1[0], spec2[0]))
    I  = np.concatenate((spec1[1], spec2[1]))
    return aggregate(mz, I)


# To do: change how the IsoSpec is called to make it faster.
class IsotopeCalculator:
    """A class for isotope calculations."""

    def __init__(   self,
                    jP = .999,
                    prec_digits = 2,
                    iso_masses = None,
                    iso_probs  = None,
                    verbose    = False     ):
        """Initiate class with information on isotopes. Calculates basic statistics of isotope frequencies: mean masses and their standard deviations."""

        if iso_masses==None or iso_probs==None:
            path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
            iso_masses, iso_probs = pickle.load(open(path+'isotopes.txt', 'rb'))
        self.iso_masses = iso_masses
        self.iso_probs  = iso_probs
        self.jP = jP
        self.elementsMassMean = dict(
            (el, sum( pr*m for pr, m in zip(self.iso_probs[el], self.iso_masses[el]) ) )
            for el in self.iso_masses.keys() )
        self.elementsMassVar  = dict(
            (el, sum( pr*m**2 for pr, m in zip(self.iso_probs[el], self.iso_masses[el])) - self.elementsMassMean[el]**2 )
            for el in self.iso_masses.keys() )
        self.isotopicEnvelopes = {}
        self.prec_digits = prec_digits
        self.formParser = formulaParser()
        self.verbose = verbose
        self.stats = Counter()

    def __getMonoisotopicMass(self, atomCnt):
        """Calculate monoisotopic mass of an atom count."""
        return sum( self.iso_masses[el][0]*elCnt for el, elCnt in atomCnt.items() )


    def __getMassMean(self, atomCnt):
        """Calculate average mass of an atom count."""
        return sum( self.elementsMassMean[el]*elCnt for el, elCnt in atomCnt.items() )


    def __getMassVar(self, atomCnt):
        """Calculate mass variance of an atom count."""
        return sum( self.elementsMassVar[el]*elCnt for el, elCnt in atomCnt.items() )


    def __getSummary(self, atomCnt_str):
        atomCnt = self.formParser.parse(atomCnt_str)
        return (    self._getMonoisotopicMass(atomCnt),
                    self._getMassMean(atomCnt),
                    self._getMassVar(atomCnt)    )


    def __getOldEnvelope(self, atomCnt_str, jP, prec_digits):
        masses, probs = self.isotopicEnvelopes[(atomCnt_str, jP, prec_digits)]
        return masses.copy(), probs.copy()


    def __getNewEnvelope(self, atomCnt_str, jP, prec_digits):
        T0 = time()
        counts          = []
        isotope_masses  = []
        isotope_probs   = []
        atomCnt = self.formParser.parse(atomCnt_str)
        for el, cnt in atomCnt.items():
            counts.append(cnt)
            isotope_masses.append(self.iso_masses[el])
            isotope_probs.append(self.iso_probs[el])
        envelope = IsoSpecPy.IsoSpec( counts, isotope_masses, isotope_probs, jP )
        masses, logprobs, _ = envelope.getConfsRaw()
        masses  = cdata2numpyarray(masses)
        probs   = np.exp(cdata2numpyarray(logprobs))
        masses, probs = agg_spec_proper(masses, probs, prec_digits)
        # memoization#TODO get rid of it when IsoSpec using stats convolution
        self.isotopicEnvelopes[ (atomCnt_str, jP, prec_digits) ] = ( masses, probs )
        T1 = time()
        self.stats['Envelopes Generation Total T'] += T1-T0

        return masses.copy(), probs.copy()


    def __getEnvelope(self, atomCnt_str, jP, prec_digits):
        if (atomCnt_str, jP, prec_digits) in self.isotopicEnvelopes:
            masses, probs = self.__getOldEnvelope(atomCnt_str,jP,prec_digits)
        else:
            masses, probs = self.__getNewEnvelope(atomCnt_str,jP,prec_digits)
        return masses, probs


    def isoEnvelope(self, atomCnt_str, jP=None, q=0, g=0, prec_digits=None):
        """Get an isotopic envelope consisting of a numpy array of masses and numpy array of probabilities.

        Parameters
        ----------
        atomCnt_str : str
            The chemical formula of a molecular species.
        jP : float
            The joint probability of the theoretical isotopic envelope.
        q : int
            The charge-state of the molecular species.
        g : int
            The quenched charge of teh molecular species.
        prec_digits : float
            The number of digits after which the floats get rounded.

        Returns
        -------
        out : tuple
            A tuple containing the theoretical spectrum: mass over charge values and intensities, both numpy arrays.
        """
        if jP is None:
            jP = self.jP
        if prec_digits is None:
            prec_digits = self.prec_digits
        masses, probs  = self.__getEnvelope(atomCnt_str, jP, prec_digits)
        if q is not 0:
            masses = np.around( (masses + g + q)/q, decimals=prec_digits )
        masses, probs = aggregate(masses, probs)
        return masses, probs


    def makeRandomSpectrum(self, mols, quants, sigma, jP=None, prec_digits=None):
        """Simulate a mixture of isotopic envelopes.

        Parameters
        ----------
        mols : list
            A list of molecular species: tuples containing (id, chemical formula string, something, charge, quenched charge)
        quants : list
            A list of total intensities of each molecular species.
        sigma : float
            The standard deviation of the masses of isotopologues - theoretical equivalent of the mass resolution.
        jP : float
            The joint probability of the theoretical isotopic envelope.
        prec_digits : float
            The number of digits after which the floats get rounded.

        Returns
        -------
        spectrum : tuple
            A tuple containing the theoretical spectrum: mass over charge values and intensities.
        """
        x0 = sum(quants)
        if not prec_digits:
            prec_digits = self.prec_digits
        if not jP:
            jP = self.jP

        def get_intensity_measure(mols, quants):
            for mol, quant in zip(mols, quants):
                _, atomCnt_str, _, q, g = mol
                ave_mz, ave_intensity = self.isoEnvelope(atomCnt_str=atomCnt_str, jP=jP, q=q, g=g, prec_digits=2)
                ave_intensity = quant * ave_intensity
                yield ave_mz, ave_intensity

        mz_average, intensity = reduce(merge_runs, get_intensity_measure(mols, quants))
        probs   = intensity/sum(intensity)
        counts  = np.array( multinomial(x0, probs), dtype='int')

        if sigma>0.0:
            spectrum = Counter()
            for m_average,cnt in zip(mz_average, counts):
                if cnt > 0:
                    m_over_z = np.round(normal(loc=m_average, scale=sigma, size=cnt), prec_digits)
                    spectrum.update(m_over_z)

            spectrum = np.array(spectrum.keys()), np.array([ float(spectrum[k]) for k in spectrum ])
        else:
            spectrum = (mz_average, counts)
        return spectrum
