PATH = /Users/matteo/Documents/MassTodon/MassTodonPy/
PATH_VISUAL = /Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/

install:
	pip install -e $(PATH)

reinstall:
	pip uninstall -y MassTodonPy
	pip install -e $(PATH)

example_call:
	python2 $(PATH)/Tests/calls/example_call.py

compare_spectra_plots:
	python2	$(PATH_VISUAL)sub_P_plot_data.py
	Rscript $(PATH_VISUAL)spectrum_fitting.R