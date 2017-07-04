storage = require '@google-cloud/storage'

currentEnv = 'test' || 'test'

conf = require "./envs/#{currentEnv}"




gcs = storage
	projectId: conf.gCloud.projectId
	credentials: conf.gCloud.credentials


# module.exports = gcs

module.exports =
	defaultBucket: gcs.bucket conf.gCloud.defaultBucket
	gcs: gcs
