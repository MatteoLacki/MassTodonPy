# MassTodonPy

![Alt](/Webpage/masstodon2.jpg)

A Python module investigating the Electron Transfer Dissociation in Mass Spectrometry

# Prerequisites

In order to play with MassTodonPy on your computer shall need:
1. UNIX based operating system, such as Linux (eg. Ubuntu, Fedora, Gentoo, or similar) or macOS (eg. macOS Sierra)
2. python2.7 interpreter

To check if you have the correct interpreter, open the terminal and simply type

```{bash}
python2.7
```
and you should get:

![Alt](/Webpage/python_terminal.png)

To leave the terminal press cntr+d

To install globally MassTodonPy, run in terminal
```{bash}
pip2 install MassTodonPy
```

This will install MassTodonPy and the required dependencies.
To check the installation, simply write now:

```{bash}
masstodon_example_call
```

This will run an example session of the program to check that you are able to get any output. By no means is it necessary to run it more than once, just after installation.

Great, you seem to have all the tools in place for some massive calculations!

There are now two ways to run the program:

1. In terminal.
2. As part of another Python script.

# Terminal Call

To run MassTodonPy in terminal, simply type

```{bash}
masstodon <spectrum_path> <config_path> -o <results_path>
```
where <spectrum_path> and config file is in <config_path> and output needs to be written to <results_path>.

The spectrum path can either lead to an mzXml spectrum, or to a simple tab separated text file containing the spectrum that in the first column contains m/z values and in second columns -- intensities. If you follow the latter option, the file could look somewhat like this:

```{bash}
191.932 17.36
271.183 98.33
415.8 17.23
425.948 15.21
444.232 208.4
444.359 6.41
444.568 117.6
445.236 44.26
449.284 19.72
... ...
```

The config file should look like this:

```{bash}
fasta = RPKPQQFFGLM
precursor_charge = 3
modification C11 = H:1, O:-1, N:1
cut_off = 100
mz_prec = .05
verbose = True
```

If you want to save a csv file with results for later inspection, add the following line to the config:

```{bash}
csv = True
```


# Python Scripting

```{python}
from MassTodonPy import MassTodonize
from MassTodonPy.TestScripts.substanceP import substanceP

mol = substanceP.copy() # some data to compare

res = MassTodonize( fasta           = mol['fasta'],
                    precursor_charge= mol['Q'],
                    mz_prec         = .05,
                    joint_probability_of_envelope = .999,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = .99,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    max_times_solve = 10,
                    raw_data        = True,
                    output_csv_path = '/Users/matteo/Documents/MassTodon/results/',
                    highcharts      = False,
                    verbose         = False )

```

# Web Service

This is coming up soon!

Remember: be nice to MassTodons.
