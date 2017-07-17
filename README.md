# MassTodonPy

![Alt](/Webpage/masstodon2.jpg)

A Python module investigating the Electron Transfer Dissociation in Mass Spectrometry

# Prerequisites

In order to play with MassTodonPy on your computer shall need:
1. UNIX based operating system, such as Linux (Ubuntu, Fedora, Gentoo, ...) or macOS (Sierra)

# assume you work on a unix with python2.7

pip2 install MassTodonPy    # install MassTodonPy in the virtual environment
masstodon_example_call   	# run an example call of the program to check that you are able to get any output

# assume your spectrum is in <spectrum_path> and config file is in <config_path> and output needs to be written to <results_path>, call
masstodon <spectrum_path> <config_path> -o <results_path>

# for csv files to appear, the config needs to contain additionally line
csv = True




### VIRTUAL ENVIRONMENT

# assume you work on a unix with python2.7 and virtualenv installed.
# assume you want to put MassTodon in folder <path>

virtualenv -p $(which python2.7) masstodonve    # prepare a python virtual environment with python2.7 interpreter
masstodonve/bin/pip install MassTodonPy			# install MassTodonPy in the virtual environment
<path>/masstodonve/bin/masstodon_example_call   # run an example call of the program

# assume your spectrum is in <spectrum_path> and config file is in <config_path> and output needs to be written to <results_path>, call
<path>/masstodonve/bin/masstodon <spectrum_path> <config_path> -o <results_path>

# for csv files to appear, the config needs to contain


./masstodon /Users/matteo/Documents/MassTodon/transport/subP_spectrum.txt /Users/matteo/Documents/MassTodon/transport/example_config_plain_text.txt
