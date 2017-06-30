PATH_VISUAL 	= ./Tests/visual
PATH_BOOTSTRAP 	= ./Tests/bootstrap
PATH_INSILICO   = ./Tests/in_silico
PYTHON   	= ../MassTodonVE/bin/python2
PIP 		= ../MassTodonVE/bin/pip2

### Installing
install_linux:
	virtualenv -p /usr/bin/python2.7 ../MassTodonVE
	$(PIP) install -e .

install_mac:
	virtualenv ../MassTodonVE
	$(PIP) install -e .

reinstall:
	$(PIP) uninstall -y MassTodonPy
	$(PIP) install -e .

check_installation:
	$(PYTHON) ./Tests/calls/check_installation.py

### Running
example_call:
	$(PYTHON) ./Tests/calls/example_call.py

compare_spectra_plots: # Relies on Rscript
	$(PYTHON) $(PATH_VISUAL)/sub_P_plot_data.py
	Rscript $(PATH_VISUAL)/spectrum_fitting.R

run_bootstrap_substance_P:
	$(PYTHON) $(PATH_BOOTSTRAP)/bootstrap.py ./Tests/bootstrap/RESULTS_CSV_27_06_2017/

analyze_bootstrap_substance_P: # Relies on Rscript
	Rscript $(PATH_BOOTSTRAP)/merge_ETnoD_PTR_plots.R

run_real_substance_P:
	$(PYTHON) $(PATH_INSILICO)/real_subP.py

# prepare_ciach_spectra:
# 	mkdir $(PATH_INSILICO)/ready_spectra/
# 	$(PYTHON) $(PATH_INSILICO)/prep_Ciach_data_4_sims.py
# 	rm -rf $(PATH_INSILICO)/ready_spectra

run_in_silico_analysis_mac:
	nice -n 10 $(PYTHON) $(PATH_INSILICO)/deconvolution.py /Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico 5

run_in_silico_analysis_wloczykij:
	nice -n 10 $(PYTHON) $(PATH_INSILICO)/deconvolution.py /home/matteo/masstodon/deconvolution/MassTodonPy/Tests/in_silico 65

run_in_silico_analysis_czczmiel:
	nice -n 10 $(PYTHON) $(PATH_INSILICO)/deconvolution.py  /home/matteo/masstodon/MassTodonPy/Tests/in_silico 25
### Cleaning
clean_ve:
	rm -rf ../MassTodonVE
