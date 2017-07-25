PATH_VISUAL 	= ./Tests/visual
PATH_INSILICO   = ./Tests/in_silico
PYTHON   	= ../MassTodonVE/bin/python2
PIP 		= ../MassTodonVE/bin/pip2

### Installing
install_linux:
	virtualenv -p /usr/bin/python2.7 ../MassTodonVE
	$(PIP) install -e .

install_quaqua:
	virtualenv -p /usr/bin/python2.7 ../MassTodonVE
	../MassTodonVE/bin/pip install -e .

install_mac:
	virtualenv ../MassTodonVE
	$(PIP) install -e .

reinstall:
	$(PIP) uninstall -y MassTodonPy
	$(PIP) install -e .

check_installation:
	$(PYTHON) Tests/calls/check_installation.py

### Running
example_call: 			## run an example session of the algorithm
	$(PYTHON) ./bin/masstodon_example_call

run_ubiquitins:			## run on all available ubiquitins
	$(PYTHON) Tests/calls/run_ubiquitins.py

run_ubiquitins_plots:			## run on all available ubiquitins for plots
	$(PYTHON) Tests/calls/run_ubiquitins_plots.py

run_substancesP_plots:			## run on all available substances P for plots
	$(PYTHON) Tests/calls/run_substancesP_plots.py


compare_spectra_plots:
	$(PYTHON) $(PATH_VISUAL)/sub_P_plot_data.py
	Rscript $(PATH_VISUAL)/spectrum_fitting.R

run_bootstrap_substance_P: 	## run statistical bootstrap analysis on all substance P spectra
	$(PYTHON) Tests/bootstrap/bootstrap_subP.py Tests/bootstrap/boot_subP_24_07_2017/
	# $(PYTHON) Tests/bootstrap/bootstrap_analysis_subP.py
	# Rscript Tests/bootstrap/analyze_bootstrap_subP.R
	# Rscript Tests/bootstrap/analyze_frag_probs_subP.R
	# Rscript Tests/bootstrap/analyze_fit_error_subP.R


run_bootstrap_ubiquitin: 	## run statistical bootstrap analysis on all ubiquitin spectra
	$(PYTHON) Tests/bootstrap/bootstrap_ubi.py Tests/data/ubiquitins.example  Tests/bootstrap/ubi_14_07_2017_mzPrec-065/
	# $(PYTHON) Tests/bootstrap/bootstrap_analysis_ubi.py
	# Rscript Tests/bootstrap/analyze_bootstrap_ubi.R
	# Rscript Tests/bootstrap/analyze_frag_probs_ubi.R
	# Rscript Tests/bootstrap/analyze_fit_error_ubi.R

SPECTRUM_PATH =
test_CLI_mac:
	node ./Tests/CLI_tests/test_input.js
	masstodon ./Tests/data/FRL-010513-SUBP-WH000-WV300.txt  ./Tests/CLI_tests/config.json -o ./Tests/CLI_tests/

# run_real_substance_P:
# 	$(PYTHON) $(PATH_INSILICO)/real_subP.py

# prepare_ciach_spectra:
# 	mkdir $(PATH_INSILICO)/ready_spectra/
# 	$(PYTHON) $(PATH_INSILICO)/prep_Ciach_data_4_sims.py
# 	rm -rf $(PATH_INSILICO)/ready_spectra

run_in_silico_analysis_mac:
	nice -n -10 $(PYTHON) $(PATH_INSILICO)/deconvolution.py /Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico 5

run_in_silico_analysis_wloczykij:
	nice -n 10 $(PYTHON) $(PATH_INSILICO)/deconvolution.py /home/matteo/masstodon/deconvolution/MassTodonPy/Tests/in_silico 65

run_in_silico_analysis_czczmiel:
	nice -n 10 $(PYTHON) $(PATH_INSILICO)/deconvolution.py  /home/matteo/masstodon/MassTodonPy/Tests/in_silico 25

### Tests
test_terminal_masstodon:  	## Test MassTodon's CLI.
	../MassTodonVE/bin/masstodon /Users/matteo/Documents/MassTodon/transport/subP_spectrum.txt /Users/matteo/Documents/MassTodon/transport/example_config_plain_text.txt -o /Users/matteo/Documents/MassTodon/transport/output


### Cleaning
clean_ve: 			## remove the python virtual environment
	rm -rf ../MassTodonVE



### Mikolaj Stuff :D
### Trolololo!!!!

projectName = masstodon

version = `cat server/package.json \
  | grep version \
  | head -1 \
  | awk -F: '{ print $2 }' \
  | sed 's/[",]//g' \
  | sed 's/[version:]//g' \
  | tr -d '[[:space:]]'`

appName = masstodon-server
appBaseName = masstodon-server-base

# docker cloud stuff
dockerAccount = masstodon

# PrebuildNames
localLatestName = $(appName):latest

remoteLatestTagName = $(dockerAccount)/$(appName):latest

remoteTestTagName     = $(dockerAccount)/$(appName):test
remoteDateTestTagName = $(dockerAccount)/$(appName):$(version)-$(shell date +"%b-%d")-test

remoteStagingTagName     = $(dockerAccount)/$(appName):staging
remoteDateStagingTagName = $(dockerAccount)/$(appName):$(version)-$(shell date +"%b-%d")-staging

remoteProdTagName     = $(dockerAccount)/$(appName):prod
remoteDateProdTagName = $(dockerAccount)/$(appName):$(version)-$(shell date +"%b-%d")-prod

run: ## run the app
	@echo 'Don't know how to run it

# --no-cache  Do not use cache when building the image.
buildCompose: ## Build / download docker images for compose DEV
	- docker-compose -f dockerfiles/dev-compose.yml build --pull --force-rm

compose: ## Start app DEV containers
	- docker-compose -f dockerfiles/dev-compose.yml -p $(projectName) up

composeProd: ## Start app PROD containers
	- docker-compose -f docker-compose-prod-test.yml -p $(projectName) up

buildBase: ## Build base app base dev image
	- docker build -t $(appBaseName) --no-cache -f dockerfiles/app-prod-base .

buildFinal: ## Build release version image
	- docker build -t $(appName) -f dockerfiles/app-prod .


pushTest: ## Push the latest test release of a prod Image
	- @echo 'Pushing as TEST. Are you sure? 3 seconds to abort. WILL AUTOREDEPLOY'
	- sleep 3
	- docker tag $(localLatestName) $(remoteLatestTagName)
	- docker tag $(localLatestName) $(remoteTestTagName)
	- docker tag $(localLatestName) $(remoteDateTestTagName)

	- sleep 5
	- docker push $(remoteLatestTagName)
	- docker push $(remoteTestTagName)
	- docker push $(remoteDateTestTagName)

	- @echo 'Done. There is no going back now'



pushStaging: ## Push the latest staging release of a prod Image
	- @echo 'Pushing as STAGING. Are you sure? 5 seconds to abort'
	- sleep 5

	- docker tag $(localLatestName) $(remoteLatestTagName)
	- docker tag $(localLatestName) $(remoteStagingTagName)
	- docker tag $(localLatestName) $(remoteDateStagingTagName)

	- docker push $(remoteLatestTagName)
	- docker push $(remoteStagingTagName)
	- docker push $(remoteDateStagingTagName)

	- @echo 'Done. There is no going back now'


pushProduction: ## Push the latest production release of a prod Image
	- @echo 'Pushing as PRODUCTION. Are you sure? 8 seconds to abort'
	- sleep 8

	- docker tag $(localLatestName) $(remoteLatestTagName)
	- docker tag $(localLatestName) $(remoteProdTagName)
	- docker tag $(localLatestName) $(remoteDateProdTagName)

	- docker push $(remoteLatestTagName)
	- docker push $(remoteProdTagName)
	- docker push $(remoteDateProdTagName)

	- @echo 'Done. There is no going back now'

attachApp: ## Attach to getty CMS backend
	- docker exec --interactive --tty $(projectName)_backend_1 /bin/bash


# -----------------------------------------------------------
# -----  EVERYTHING BELOW THIS LINE IS NOT IMPORTANT --------
# -----       (Makefile helpers and decoration)      --------
# -----------------------------------------------------------
#
# Decorative Scripts - Do not edit below this line unless you know what it is

.PHONY: help
.DEFAULT_GOAL := help

NO_COLOR    = \033[0m
INDENT      = -30s
BIGINDENT   = -50s
GREEN       = \033[36m
BLUE        = \033[34m
DIM         = \033[2m
help:
	@printf '\n\n$(DIM)Commands:$(NO_COLOR)\n'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "$(GREEN) % $(BIGINDENT)$(NO_COLOR) - %s\n", $$1, $$2}'
