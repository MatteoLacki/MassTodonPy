_ = require 'lodash'
needle = require 'needle'


projectName = 'masstodon-test'
promisify = require 'es6-promisify'

myGet = promisify needle.get



# https://masstodon-test.firebaseio.com/work/work12.json

# /work/test


dataUrl = (path) ->
	"https://#{projectName}.firebaseio.com#{path}.json"

fetchData = (path) ->
	myGet dataUrl path


testPath = '/work/work12'

try
	await res = fetchData testPath
	console.log res
	console.log "Downloading data from: #{dataUrl testPath}"
	console.log 'Loaded config: ', res.body
catch e
	console.log 'some err', e


# a.catch (e) ->

# 	console.log 'EERRRR', e
