bin = ./node_modules/.bin

# Test stuff

NODETESTS = `find src -name 'spec-node.coffee'`
# spec, dot, progress, tap, landing, nyan, min (mocha --reporters)
REPORTER = spec

COFFEE =  coffee:coffee-script/register
# COFFEE =  coffee:coffeescript/register
COFFEE2 =  coffee:coffeescript/register


test: ## Run funny tests
	@ NODE_ENV=test \
	  $(bin)/mocha-webpack --require ./test-setup.js --ui exports --webpack-config webpack.config.coffee --glob \"*spec.coffee\" --recursive src --interactive false

run: ## Just run a server
	@ $(bin)/coffee \
		src/index.coffee

# easyNewTest: ## Run basic tests in an easy node setup
# 	@ echo $(NODETESTS)
# 	@ $(bin)/mocha \
# 		--require should \
# 		--require ignore-styles \
# 		--compilers $(COFFEE2) \
# 		--reporter $(REPORTER) \
# 		$(NODETESTS)

# oldTest: ## Run basic tests in an easy node setup
# 	@ echo $(NODETESTS)
# 	@ $(bin)/mocha \
# 		--require should \
# 		--require ignore-styles \
# 		--compilers $(COFFEE) \
# 		--reporter $(REPORTER) \
# 		$(NODETESTS)


# watchEasy: ## Run basic tests in an easy node setup
# 	@ echo $(NODETESTS)
# 	@ $(bin)/mocha \
# 		--require should \
# 		--require ignore-styles \
# 		--compilers $(COFFEE) \
# 		--watch \
# 		--reporter spec \
# 		$(NODETESTS)



# start: ## Start a normal developer environment
# 	@ NODE_ENV=develop \
# 	  $(bin)/coffee server.coffee

# develop: ## Start a normal developer environment (like start)
# 	@ NODE_ENV=develop \
# 	  $(bin)/coffee server.coffee

# dashboard: ## Develop with a dashboard
# 	@ NODE_ENV=develop \
# 	  $(bin)/webpack-dashboard -- coffee server.coffee --env=develop


# build: ## production build
# 	- rimraf dist/*
# 	@ NODE_ENV=production \
# 	  $(bin)/webpack

# storybook: ## Run storybook on a 9001 port
# 	@ $(bin)/start-storybook -p 9001 -c .storybook

# lint: ## Run a coffeelinter
# 	@ $(bin)/coffeelint -f coffeelint.json src

# vendorStuff: ## Prepare vendor stuff as a DLL :D
# 	@ $(bin)/webpack --config webpack.vendor.config.js












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
YELLOW      = \033[33m

# @printf '$(DIM)Variables$(NO_COLOR):\n'
# @printf "$(DIM)$(BLUE) % $(BIGINDENT)$(NO_COLOR) %s \n" "Secrets folder name:" ./$(secretsFolderName)
help:
	@printf '$(DIM)Commands:$(NO_COLOR)\n'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "$(GREEN) % $(BIGINDENT)$(NO_COLOR) - %s\n", $$1, $$2}'
