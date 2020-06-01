'''

Intelligent agent

Link to derivations
https://en.wikipedia.org/wiki/Intelligent_agent

'''

# Some dummy agents to use in experimentation

class Agent(object):
	GLOBAL_DUMMY_AGENT_NAME_LIST = [
		["Alice", "The first agent"],
		["Bob", "The second agent"],
		["Charlie", "The third agent"],
		["Derrick", "The fourth agent"],
	]
	GLOBAL_DUMMY_AGENT_NAME_LIST_COUNTER = 0
	def __init__(self, **kwargs):
		self.identity = Agent.GLOBAL_DUMMY_AGENT_NAME_LIST[Agent.GLOBAL_DUMMY_AGENT_NAME_LIST_COUNTER%len(Agent.GLOBAL_DUMMY_AGENT_NAME_LIST)]
		self.name = kwargs.get('name',self.identity[0])
		self.description = kwargs.get('description',self.identity[1])
		self.pot = kwargs.get('pot', 0)
		self.trace_pot = [0]
		self.expected_utility = {}
		self.last_move = None
		self.response_function = kwargs.get('response_function',None)
		if (kwargs.get('name',None) == None):
			Agent.GLOBAL_DUMMY_AGENT_NAME_LIST_COUNTER += 1