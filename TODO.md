### MassTodon todo list

# important
* Update to IsoSpec 2.0
* Get rid of the Pandas dependency for csv saving (sic!)
* Get rid of highcharts -> turn it into Bokeh
* check different modes of isotopic calculations under IsoSpec2.0
* add csv input spectrum
* add assertions to read_mzxml_spectrum
* Replace the some arguments for MassTodonize with dictionaries.
    * e.g. solver_args, deconvolutor_args
* Replace ubiquitin dataset
    * now, the intensity if monotonically increasing
* RE must parse H-200 -> {'H': -200}
* Drawing of
    * spectra
    * deconvolion
* Cannot pass different non-default args to mol.isotopologues()
    * e.g. mol.isotopologues(5, .99) seems to loop like hell.
    * might be because of some issues with aggregation.

# elegant
* Automate devel set up
* Get lavish estimates

# Don't forget
* The bloody Proline has precursor.get_AA(4,'C_alpha') == lCnt()
* The +1H makes part of the c-fragment definition
    * think about the possible products.
* Get rid of spurious dependencies in setup.py

# IsoSpec:
* Get a version that does not need string parsing.
* no need to be wrapped in 'cdata2numpyarray'.
* rounded m/z to some precision.
* with probabilities not logprobabilities.



# Some long term stuff:
* when more than one precursor
    * exception for the same precursor tags in make_molecules
        * Rationale: otherwise we will not be able to trace the origin of a fragment
        * Problem: what if two substances point to the same formula?
            * check hashes



============================================================================
# Done:
* Both:
    * D = Deconvolutor(molecules, spectrum, L1_x=0.002, L2_x=0.002, L1_alpha=0.002, L2_alpha=0.002)
    * D = Deconvolutor(molecules, spectrum)
    * work
* subclass linearCounter into the project and rename it as atom_cnt.
    * add monoisotopic mass
    * add IsotopeCalculator as a subroutine here: not bad! This can be a class field.
        * The Formula class should initiate a field called isotopic generator.