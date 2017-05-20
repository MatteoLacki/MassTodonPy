_ = require 'lodash'
axios = require 'axios'


projectName = 'masstodon-test'
dataUrl = (path) ->
	"https://#{projectName}.firebaseio.com#{path}.json"

fetchData = (path) ->
	await axios.get dataUrl path


testPath = '/work/work12'
try
	fetchData(testPath).then (c) ->
		console.log 'Got c: ', c.data
catch e
	console.log 'Failed to load data'
