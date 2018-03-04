PYTHON_2   	= ../MassTodonVEpy2/bin/python2
PYTHON_3   	= ../MassTodonVEpy3/bin/python3
NOSE_2 		= ../MassTodonVEpy2/bin/nosetests
NOSE_3 		= ../MassTodonVEpy3/bin/nosetests

masstodon_2: ## Install MassTodon virtual environment with python2
	virtualenv -p python2 ../MassTodonVEpy2
	../MassTodonVEpy2/bin/pip2 install -e .
	../MassTodonVEpy2/bin/python2 install ipykernel

masstodon_3:
	virtualenv -p python3 ../MassTodonVEpy3
	../MassTodonVEpy3/bin/pip3 install -e .
	../MassTodonVEpy3/bin/python3 install ipykernel

reinstall_2:
	../MassTodonVEpy2/bin/pip2 uninstall -y MassTodonPy
	pip2 install -e .

reinstall_3:
	../MassTodonVEpy3/bin/pip3 uninstall -y MassTodonPy
	pip3 install -e .

call_2: ## run an example session of the algorithm
	$(PYTHON_2) ./bin/masstodon_example_call

call_3: ## run an example session of the algorithm
	$(PYTHON_3) ./bin/masstodon_example_call

tests:  ## run all possible tests.
	$(NOSE_2) MassTodonPy
	$(NOSE_3) MassTodonPy

tests2:
	$(NOSE_2)

tests3:
	$(NOSE_3)

clean_masstodon_2: ## remove virtual environment with python2
	rm -rf ../MassTodonVEpy2

clean_masstodon_3: ## remove virtual environment with python2
	rm -rf ../MassTodonVEpy3

pypi: ## send to pypi
	 git tag -d GutenTag
	 git push gh :refs/tags/GutenTag

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
