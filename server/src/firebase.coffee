
# module.exports =
# 	fireService: {//firebase service json key}

currentEnv = 'test' || 'test'

conf = require "./envs/#{currentEnv}"


admin = require 'firebase-admin'
admin.initializeApp
	credential: admin.credential.cert conf.firebaseService
	databaseURL: conf.fireUri


module.exports = admin
