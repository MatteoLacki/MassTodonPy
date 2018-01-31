### MassTodon todo list

# important
+ add reading from '.json'
+ update docs and readmes.
+ implement back the multiprocessing version of the software.
    + it is an argument for the graph representation
+ eliminate the alphas and raw estimates
    + we are saving them with the reporter class
+ csv:
    + input spectrum
    * save results
    * no pandas
+ Update to IsoSpec 2.0
+ check different modes of isotopic calculations under IsoSpec2.0
+ add assertions to read_mzxml_spectrum
+ Replace ubiquitin dataset
    * now, the intensity if monotonically increasing
* RE must parse H-200 -> {'H': -200}
+ Drawing of
    * spectra
    + deconvolution
        + Add the simple plot (observed-predicted) to the module
        + Add the complex plot (observed-predicted-per-molecule)
+ Have a look at IsotopeCalculator/simulator.py
    + DO WE USE SPECTRUM CLASS HERE? We certainly could.
+ Spectrum reader must be run without sorting the spectra all the time.
+ Bash scripts:
    + plot_spectrum to open any mzXml, mzMl, txt spectra file
    + IsoSpecPy plot tools
+ Export data to ETDetective
+ Support multiple input precursors!!!

# elegant
+ implement the additional 1D-regression-test for fragment inclusion.
+ Memoization of isotopic distributions should be an option, not a must.
+ Automate devel set up
+ Add the invisible buffers to Measure.plot()
+ add another intensity-based criterion here.
    + basically, check how much of a substance could there be, if the was the only possible source of these ions. and accept if it is more than some number. This way silly solutions should be eliminated.

# Don't forget
+ The bloody Proline has precursor.get_AA(4,'C_alpha') == lCnt()
+ The +1H makes part of the c-fragment definition
    + think about the possible products.
+ Get rid of spurious dependencies in setup.py

# IsoSpec:
+ Get a version that does not need string parsing.
+ no need to be wrapped in 'cdata2numpyarray'.
+ rounded m/z to some precision.
+ with probabilities not logprobabilities.

# Some long term stuff:
+ when more than one precursor
    + exception for the same precursor tags in make_molecules
        + Rationale: otherwise we will not be able to trace the origin of a fragment
        + Problem: what if two substances point to the same formula?
            + check hashes



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
* Cannot pass different non-default args to mol.isotopologues()
    * e.g. mol.isotopologues(5, .99) seems to loop like hell.
    * might be because of some issues with aggregation.
* Get rid of highcharts -> turn it into Bokeh
