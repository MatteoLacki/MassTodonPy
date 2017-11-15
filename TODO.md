### MassTodon todo list

# important
* Update the networkx module usage while retaining backward compatibility
* Update the IsoSpecPy module usage
* Get rid of the Pandas dependency for csv saving (sic!)
* Get rid of highcharts -> turn it into Bokeh
* exception for the same precursor tags in make_molecules
    * Rationale: otherwise we will not be able to trace the origin of a fragment
    * Problem: what if two substances point to the same formula?

# elegant
* Get rid of the superatoms list. Should be a generator.
* replace make_superatoms with amino_acids(left or right)
* Precursor should supply the get_super_atoms method.
    * Two iterators actually: left to right and right to left
* Change the storage of amino acid bricks.

# Less important
