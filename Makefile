PATH_TO_PROJECT = /Users/matteo/Documents/MassTodon/MassTodonPy/
PATH_VISUAL 	= /Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/
PATH_BOOTSTRAP 	= /Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/


install:
	pip install -e $(PATH_TO_PROJECT)

reinstall:
	pip uninstall -y MassTodonPy
	pip install -e $(PATH_TO_PROJECT)

example_call:
	python2 $(PATH_TO_PROJECT)Tests/calls/example_call.py

compare_spectra_plots:
	python2	$(PATH_VISUAL)sub_P_plot_data.py
	Rscript $(PATH_VISUAL)spectrum_fitting.R

run_bootstrap_substance_P:
	python2 $(PATH_BOOTSTRAP)bootstrap_subP.py

analyze_bootstrap_substance_P:
	Rscript $(PATH_BOOTSTRAP)merge_ETnoD_PTR_plots.R

run_real_substance_P:
	python2 $(PATH_BOOTSTRAP)real_subP.py
