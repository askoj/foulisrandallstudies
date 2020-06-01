
import math as mt

'''

	Returns the ID of the option that has the highest payoff

'''
class highPay(object):
	def __init__(self, **kwargs):
		self.name = "Maximisation Of Utility"

	def run(self,**kwargs):
		game = kwargs.get('game', None)
		agent = kwargs.get('agent', None)
		option_payoffs = {}
		for option in game.options:
			option_payoffs.update({ option.id : 0 })

		for k,v in game.payoffs.items():
			option_payoffs[v[game.agent_index(agent)].id] += agent.expected_utility[k]
		highest_payoff = -mt.inf
		nominated_option_id = None
		for k, v in option_payoffs.items():
			if v > highest_payoff:
				highest_payoff = v
				nominated_option_id = k
		nominated_option = None
		for option in game.options:
			if (nominated_option_id == option.id):
				nominated_option = option
		return nominated_option

maximiseUtility = highPay()

'''

	Returns the ID of the option that has the highest payoff with respect to K

'''
class highPayByK(object):
	def __init__(self, **kwargs):
		self.name = "Maximisation Of Utility (K Sensitive)"

	def run(self,**kwargs):
		game = kwargs.get('game', None)
		agent = kwargs.get('agent', None)
		agent_k = agent.k
		if kwargs.get('fitted', False):
			agent_k = agent.fitted_k
		option_payoffs = {}
		for option in game.options:
			option_payoffs.update({ option.id : 0 })
		for k,v in game.payoffs.items():
			option_payoffs[v[game.agent_index(agent)].id] += agent.expected_utility[k]
		nominated_option_id = None
		sorted_option_payoffs = sorted(option_payoffs.items(), key=lambda kv: kv[1])
		chosen_index = round((len(sorted_option_payoffs)-1)*agent_k)
		nominated_option_id = sorted_option_payoffs[chosen_index][0]
		nominated_option = None
		for option in game.options:
			if (nominated_option_id == option.id):
				nominated_option = option
		return nominated_option

maximiseUtilityByK = highPayByK()