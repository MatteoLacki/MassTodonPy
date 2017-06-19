PROJ_PATH = ./
VE_PATH   = ../MassTodonVE/bin/
PATH_VISUAL 	= $(PROJ_PATH)/Tests/visual/
PATH_BOOTSTRAP 	= $(PROJ_PATH)/Tests/bootstrap/

### Installing
install:
	virtualenv ../MassTodonVE
	$(VE_PATH)pip install -e $(PROJ_PATH)

reinstall:
	$(VE_PATH)pip uninstall -y MassTodonPy
	$(VE_PATH)pip install -e $(PROJ_PATH)


### Running
example_call:
	$(VE_PATH)python2 $(PROJ_PATH)Tests/calls/example_call.py

compare_spectra_plots: # Relies on Rscript
	$(VE_PATH)python2 $(PATH_VISUAL)sub_P_plot_data.py
	Rscript $(PATH_VISUAL)spectrum_fitting.R

run_bootstrap_substance_P:
	$(VE_PATH)python2 $(PATH_BOOTSTRAP)bootstrap_subP.py

analyze_bootstrap_substance_P:
	$(VE_PATH)Rscript $(PATH_BOOTSTRAP)merge_ETnoD_PTR_plots.R

run_real_substance_P:
	$(VE_PATH)python2 $(PATH_BOOTSTRAP)real_subP.py


### Cleaning
clean_ve:
	rm -rf ../MassTodonVE
