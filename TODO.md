### MassTodon todo list

# important
* Update the networkx module usage while retaining backward compatibility
* Update the IsoSpecPy module usage
* Get rid of the Pandas dependency for csv saving (sic!)
* Get rid of highcharts -> turn it into Bokeh
* exception for the same precursor tags in make_molecules
    * Rationale: otherwise we will not be able to trace the origin of a fragment
    * Problem: what if two substances point to the same formula?
* check different modes of isotopic calculations under IsoSpec2.0
* apply the joint probability thresholding to real spectrum.
* add csv input spectrum
* add assertions to read_mzxml_spectrum
* compare the new spectra with the older ones.
* Replace the some arguments for MassTodonize with dictionaries.
    * e.g. solver_args, deconvolutor_args
* Replace the current ubiquitin dataset with one that is not corrupted
    * this one has an ever growing intensity, for some reason.
* Check if adding works well in the LinearCounter and iclude it as a module.

# elegant
* Import linearCounter into the project and rename it as atom_cnt.
    * add __str__ method (replacement for atom_cnt_2_string)
    * add monoisotopic mass
    * add IsotopeCalculator as a subroutine here.
    * add possibility to by-pass memoization in IsotopeCalculator

# Less important


# Don't forget
* The bloody Proline has precursor.get_AA(4,'C_alpha') == lCnt()
* The +1H makes part of the c-fragment definition
    * think about the possible products.
