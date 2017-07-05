_ = require 'lodash'
axios = require 'axios'

# --- test code
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

# --- Worker
path = require 'path'
mkdirp = require 'mkdirp'

workDir = path.resolve process.cwd(), '.tmp'
console.log workDir
mkdirp workDir, (e) ->
	console.log 'Ensured', workDir

gcs = require './gcs'

myBucket = gcs.defaultBucket

# --- Runner

firebase = require './firebase'
jobs = firebase.database().ref 'jobs'


jobs.on 'value', (snap) ->
	allJobs = snap.val()
	console.log allJobs

	_.forEach allJobs, (job, k) ->
		if !job.processed
			jobDir = path.resolve workDir, job.owner, 'input'

			# Ensure folder and download input file
			mkdirp jobDir, (e) ->
				console.log 'job dir ensured: ', jobDir

				myBucket


				jobs.child(k).child('processed').set(true).then (e) ->
					console.log 'UPDATED JOB'
		else
			console.log 'job is done', k
