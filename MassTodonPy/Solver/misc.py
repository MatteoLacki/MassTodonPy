class Indexer(dict):
	def get(self, v1, v2):
		args = frozenset( (v1,v2) )
		return super(Indexer, self).__getitem__(args)
	def add(self, v1, v2, val):
		args = frozenset( (v1,v2) )
		super(Indexer, self).__setitem__(args, val)
