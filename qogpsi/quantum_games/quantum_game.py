class Quantum_Game(object):
	def __init__(self, **kwargs):
		self.initial_state = kwargs.get('state',None)
		self.state = self.initial_state