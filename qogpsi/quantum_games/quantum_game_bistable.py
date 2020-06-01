from qutip import *
import numpy as np
import math as mt
import cmath as cm
import sympy as sm
import matplotlib
import re
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt

from qogpsi.stylization.themes import display_theme
from qogpsi.reporting.logs import log_title, log_subtitle, log_dictionary, log_text
from qogpsi.quantum_games.quantum_game import Quantum_Game
from qogpsi.multi_agent_systems.agent import Agent
from qogpsi.multi_agent_systems.option import Option

from random import choice, randrange
from string import ascii_lowercase

ii = cm.sqrt(-1)
CONST_K_UPPER_BOUND = 1
CONST_K_LOWER_BOUND = 0
CONST_THETA_UPPER_BOUND = mt.pi
CONST_THETA_LOWER_BOUND = 0
CONST_PHI_UPPER_BOUND = mt.pi
CONST_PHI_LOWER_BOUND = 0
CONST_SIGMA_UPPER_BOUND = mt.pi
CONST_SIGMA_LOWER_BOUND = 0
CONST_RHO_UPPER_BOUND = mt.pi
CONST_RHO_LOWER_BOUND = 0

CONST_FAIR_GAME = (ket("00")+(ii*ket("11")))+(ket("01")+(ii*ket("10")))
CONST_CORRELATED_GAME = (ket("00")+(ii*ket("11")))

def randomize_parameter(val):
	if val == "k":
		return ((randrange(0,101,1)/100)*(CONST_K_UPPER_BOUND-CONST_K_LOWER_BOUND))+CONST_K_LOWER_BOUND
	if (val == "theta"):
		return ((randrange(0,101,1)/100)*(CONST_THETA_UPPER_BOUND-CONST_THETA_LOWER_BOUND))+CONST_THETA_LOWER_BOUND
	if (val == "phi"):
		return ((randrange(0,101,1)/100)*(CONST_PHI_UPPER_BOUND-CONST_PHI_LOWER_BOUND))+CONST_PHI_LOWER_BOUND
	if (val == "sigma"):
		return ((randrange(0,101,1)/100)*(CONST_SIGMA_UPPER_BOUND-CONST_SIGMA_LOWER_BOUND))+CONST_SIGMA_LOWER_BOUND
	if (val == "rho"):
		return ((randrange(0,101,1)/100)*(CONST_RHO_UPPER_BOUND-CONST_RHO_LOWER_BOUND))+CONST_RHO_LOWER_BOUND

def determine_entanglement(ketstate):
	rho = ket2dm(ketstate)
	sysy = tensor(sigmay(), sigmay())
	rho_tilde = Qobj(np.array(rho))*Qobj(np.array(sysy)) * Qobj(np.array(rho.conj()))*Qobj(np.array(sysy))
	evals = abs(np.sort(np.real(rho_tilde.eigenenergies())))
	lsum = np.sqrt(evals[3]) - np.sqrt(evals[2]) - np.sqrt(evals[1]) - np.sqrt(evals[0])
	return max(0, lsum)

def correlation_to_distribution(corr):
    diag_a = 0.0
    diag_b = 0.0
    if corr > 0:
        diag_a = abs(corr)*0.5
        diag_b = abs(1-abs(corr))*0.5
    else:
        diag_b = abs(corr)*0.5
        diag_a = abs(1-abs(corr))*0.5
    distribution = np.array([	[float('%.3f' % diag_a), float('%.3f' % diag_b)],
                                [float('%.3f' % diag_b), float('%.3f' % diag_a)] ])
    distribution /= (sum(sum(distribution)))
    return distribution

def extract_correlation(q_obj):
    return np.sign(q_obj[0][0][0])*((q_obj[0][0][0].real**2) + (q_obj[0][0][0].imag**2))

def contextuality_dz(Aa_Bb, Aa_Bb_pr, Aa_pr_Bb, Aa_pr_Bb_pr):
    p1 = Aa_Bb[0][0]
    p2 = Aa_Bb[0][1]
    p3 = Aa_Bb[1][0]
    p4 = Aa_Bb[1][1]
    p5 = Aa_Bb_pr[0][0]
    p6 = Aa_Bb_pr[0][1]
    p7 = Aa_Bb_pr[1][0]
    p8 = Aa_Bb_pr[1][1]
    p9  = Aa_pr_Bb[0][0]
    p10 = Aa_pr_Bb[0][1]
    p11 = Aa_pr_Bb[1][0]
    p12 = Aa_pr_Bb[1][1]
    p13 = Aa_pr_Bb_pr[0][0]
    p14 = Aa_pr_Bb_pr[0][1]
    p15 = Aa_pr_Bb_pr[1][0]
    p16 = Aa_pr_Bb_pr[1][1]

    A_11 = (2*(p1+p2)) - 1
    A_12 = (2*(p5+p6)) - 1
    A_21 = (2*(p9+p10)) - 1
    A_22 = (2*(p13+p14)) - 1

    B_11 = (2*(p1+p3)) - 1
    B_12 = (2*(p9+p11)) - 1
    B_21 = (2*(p5+p7)) - 1
    B_22 = (2*(p13+p15)) - 1

    delta_0 = (abs(A_11 - A_12) + abs(A_21 - A_22) + abs(B_11 - B_21) + abs(B_12 - B_22))/2

    A_11_B_11 = (p1 + p4) - (p2 + p3)
    A_12_B_12 = (p5 + p8) - (p6 + p7)
    A_21_B_21 = (p9 + p12) - (p10 + p11)
    A_22_B_22 = (p13 + p16) - (p14 + p15)

    l1 = abs(+ A_11_B_11 + A_12_B_12 + A_21_B_21 - A_22_B_22)
    l2 = abs(+ A_11_B_11 + A_12_B_12 - A_21_B_21 + A_22_B_22)
    l3 = abs(+ A_11_B_11 - A_12_B_12 + A_21_B_21 + A_22_B_22)
    l4 = abs(- A_11_B_11 + A_12_B_12 + A_21_B_21 + A_22_B_22)

    sd = 2*(1+delta_0)

    va = (l1 - sd) if (l1 >= sd) else 0
    vb = (l2 - sd) if (l2 >= sd) else 0
    vc = (l3 - sd) if (l3 >= sd) else 0
    vd = (l4 - sd) if (l4 >= sd) else 0
    return (va+vb+vc+vd)

def correlative_measure(Φ_state):
	Sx = sigmax()
	Sz = sigmaz()
	A_a1 = tensor( Sz, identity(1))
	A_a2 = tensor( Sx, identity(1))
	B_b1 = (-1/mt.sqrt(2))*tensor(identity(1),(Sz+Sx))
	B_b2 = (1/mt.sqrt(2))*tensor(identity(1),(Sz-Sx))
	A_a1_B_b1 = Qobj(Φ_state.data).dag()*Qobj(tensor(A_a1,B_b1).data)*Qobj(Φ_state.data)
	A_a2_B_b1 = Qobj(Φ_state.data).dag()*Qobj(tensor(A_a2,B_b1).data)*Qobj(Φ_state.data)
	A_a1_B_b2 = Qobj(Φ_state.data).dag()*Qobj(tensor(A_a1,B_b2).data)*Qobj(Φ_state.data)
	A_a2_B_b2 = Qobj(Φ_state.data).dag()*Qobj(tensor(A_a2,B_b2).data)*Qobj(Φ_state.data)
	Aa_Bb = correlation_to_distribution(extract_correlation(A_a1_B_b1))
	Aa_Bb_pr = correlation_to_distribution(extract_correlation(A_a1_B_b2))
	Aa_pr_Bb = correlation_to_distribution(extract_correlation(A_a2_B_b1))
	Aa_pr_Bb_pr = correlation_to_distribution(extract_correlation(A_a2_B_b2))
	return ((A_a1_B_b1+A_a2_B_b1+A_a1_B_b2+A_a2_B_b2)[0][0][0]).real

def contextuality_measure(Φ_state):
	Sx = sigmax()
	Sz = sigmaz()
	A_a1 = tensor( Sz, identity(1))
	A_a2 = tensor( Sx, identity(1))
	B_b1 = (-1/mt.sqrt(2))*tensor(identity(1),(Sz+Sx))
	B_b2 = (1/mt.sqrt(2))*tensor(identity(1),(Sz-Sx))
	A_a1_B_b1 = Qobj(Φ_state.data).dag()*Qobj(tensor(A_a1,B_b1).data)*Qobj(Φ_state.data)
	A_a2_B_b1 = Qobj(Φ_state.data).dag()*Qobj(tensor(A_a2,B_b1).data)*Qobj(Φ_state.data)
	A_a1_B_b2 = Qobj(Φ_state.data).dag()*Qobj(tensor(A_a1,B_b2).data)*Qobj(Φ_state.data)
	A_a2_B_b2 = Qobj(Φ_state.data).dag()*Qobj(tensor(A_a2,B_b2).data)*Qobj(Φ_state.data)
	Aa_Bb = correlation_to_distribution(extract_correlation(A_a1_B_b1))
	Aa_Bb_pr = correlation_to_distribution(extract_correlation(A_a1_B_b2))
	Aa_pr_Bb = correlation_to_distribution(extract_correlation(A_a2_B_b1))
	Aa_pr_Bb_pr = correlation_to_distribution(extract_correlation(A_a2_B_b2))
	return contextuality_dz(Aa_Bb, Aa_Bb_pr, Aa_pr_Bb, Aa_pr_Bb_pr)

def generate_string(n=10):
	return "".join(choice(ascii_lowercase) for i in range(n))

def scalar(v):
	# omitted this
    #return mt.sqrt(v.real**2 + v.imag**2)
    return mt.sqrt(v.real**2 + v.imag**2)*np.sign(v.real)

class Feasibility():
	def __init__(self, **kwargs):
		self.id = generate_string()
		self.agent_names = []
		self.agent_response_functions = []
		self.agent_fitted_option_names = []
		self.snapshot_options = kwargs.get('snapshot_options', [])
		self.snapshot_state_of_game = kwargs.get('snapshot_state_of_game', (((ket("00")+(ii*ket("11")))+(ket("01")+(ii*ket("10"))))/1))
		self.agents_list = kwargs.get('agents', [])
		for agent in self.agents_list:
			self.agent_names.append(agent.name)
			if agent.fitted_option is not None:
				self.agent_fitted_option_names.append(agent.fitted_option.name)
			self.agent_response_functions.append(agent.response_function.name)

		self.snapshot_payoffs = {}
		for k, v in kwargs.get('snapshot_payoffs', {}).items():
			self.snapshot_payoffs.update({ k : "%s chooses %s & %s chooses %s \n\t\t\t %s is paid %s and %s is paid %s" % (
				self.agent_names[0],
				v[0].name,
				self.agent_names[1],
				v[1].name,
				self.agent_names[0],
				v[2],
				self.agent_names[1],
				v[3]) })

		self.fitted_agent_ks = kwargs.get('fitted_agent_ks', [])
		self.fitted_agent_thetas = kwargs.get('fitted_agent_thetas', [])
		self.fitted_agent_phis = kwargs.get('fitted_agent_phis', [])
		self.fitted_agent_rhos = kwargs.get('fitted_agent_rhos', [])
		self.fitted_agent_sigmas = kwargs.get('fitted_agent_sigmas', [])
		self.agent_parameter_settings = kwargs.get('agent_parameter_settings', {})
		self.feasibilities = kwargs.get('feasibilities', {})
	def get(self, argument):
		if (argument == "number_of_feasibilities"):
			return self.feasibilities.count(1)
		if (argument == "game_state_entanglement_a"):
			aaa = (Qobj(np.array(tensor(ket("00"),bra("00"))))+Qobj(np.array(tensor(ket("11"),bra("11")))))*Qobj(np.array(self.snapshot_state_of_game))
			return determine_entanglement(aaa)
		if (argument == "game_state_entanglement_b"):
			bbb = (Qobj(np.array(tensor(ket("01"),bra("01"))))+Qobj(np.array(tensor(ket("10"),bra("10")))))*Qobj(np.array(self.snapshot_state_of_game))
			return determine_entanglement(bbb)
		if (argument == "game_state_entanglement_c"):
			return determine_entanglement(self.snapshot_state_of_game)

class Bistable_Quantum_Game(Quantum_Game):
	def __init__(self, **kwargs):
		self.id = generate_string()
		self.initial_state = kwargs.get('state', (ket("00")+(ii*ket("11"))))#(((ket("00")+(ii*ket("11")))+(ket("01")+(ii*ket("10"))))/1))
		self.trace_entanglement = [determine_entanglement(self.initial_state)]
		self.trace_correlation = [correlative_measure(self.initial_state)]
		self.trace_contextuality = [contextuality_measure(self.initial_state)]
		
		self.state = self.initial_state
		self.agents = kwargs.get('agents', [])
		self.options = kwargs.get('options', [])
		self.payoffs = kwargs.get('payoffs', { "α" : 3, "β" : 0, "γ" : 5, "δ" : 1 })
		self.DEFAULT_PRECISION_IN_GRID_SEARCH = 100
		self.break_fitting = False
		self.agent_parameter_settings = []
		self.fitted_conditions = []
		self.agent_parameter_settings_traversal = {}
		self.fitted_agent_thetas = []
		self.fitted_agent_phis = []
		self.fitted_agent_rhos = []
		self.fitted_agent_sigmas = []
		self.fitted_agent_ks = []
		self.feasibilities = []
		Agent.GLOBAL_DUMMY_AGENT_NAME_LIST_COUNTER = 0

	'''

		Set the game's state entanglement (beta function)

		Theta needs to be in radians

		Presently, only entangles kets of 00 and 11

		Concurrence is used as the measure of entanglement

	

	def entanglement(self, theta=None):
		if theta is not None:
			self.state = ((((mt.cos(theta/2)*ket("00"))+(mt.sin(theta/2)*ket("11")))+((mt.cos(theta/2)*ket("10"))-(mt.sin(theta/2)*ket("01"))))/1)
		else:
			return determine_entanglement(self.state)
	'''

	def agent_index(self, agent):
		this_agent_index = -1
		for i in range(0, len(self.agents)):
			if (agent == self.agents[i]):
				this_agent_index = i
		return this_agent_index

	def option_from_id(self, argument):
		for option in self.options:
			if (option.id == argument):
				return option

	def kraus_projections(self, kA, kB, options):
		# The current Kraus projection structure restricts a game to two agents (and two options)
		# (we will need to build this theory further)
		a1 = (kA*kB)+((1-kA)*(1-kB))
		a2 = (kB*(1-kA))+(kA*(1-kB))
		a3 = ((2*kA)-1)*((2*kB)-1)
		kPcc = Qobj(np.array([
			[     a1,      0,      0, -ii*a3 ],
			[      0,     a2,      0,      0 ],
			[      0,      0,     a2,      0 ],
			[  ii*a3,      0,      0,     a1 ]]))/2
		kPcd = Qobj(np.array([
			[     a2,      0,      0,      0 ],
			[      0,     a1,  ii*a3,      0 ],
			[      0, -ii*a3,     a1,      0 ],
			[      0,      0,      0,     a2 ]]))/2
		kPdc = Qobj(np.array([
			[     a2,      0,      0,      0 ],
			[      0,     a1, -ii*a3,      0 ],
			[      0,  ii*a3,     a1,      0 ],
			[      0,      0,      0,     a2 ]]))/2
		kPdd = Qobj(np.array([
			[     a1,      0,      0,  ii*a3 ],
			[      0,     a2,      0,      0 ],
			[      0,      0,     a2,      0 ],
			[ -ii*a3,      0,      0,     a1 ]]))/2
		projection_map = {
			options[0].id+options[0].id : kPcc,
			options[0].id+options[1].id : kPcd,
			options[1].id+options[0].id : kPdc,
			options[1].id+options[1].id : kPdd
		}
		return projection_map

	def compile_payoffs_map(self, arg_payoffs):
		payoffs_map = {}
		payoff_identities = []
		for k,v in arg_payoffs.items():
			payoffs_map.update({ v[0].id+v[1].id : [k,v[2],v[3]] })
			payoff_identities.append(v[0].id+v[1].id)
		return [payoffs_map, payoff_identities]

	def determine_outcome(self, arg_agents, payoffs_map, payoff_identities, projection_map, **kwargs):
		game_state = kwargs.get('game_state', self.state)
		outcome = []
		for agent in arg_agents:
			agent.expected_utility = {}
			for payoff_id in payoff_identities:
				# elements are multiplied by their complex conjugates, hence why we can simply reduce the complex numbers to their real parts
				agent.expected_utility.update({ payoffs_map[payoff_id][0] : 
					scalar(np.array(payoffs_map[payoff_id][self.agent_index(agent)+1]*Qobj(game_state.dag().data)*projection_map[payoff_id]*Qobj(game_state.data))[0][0]) })
			agent_decision = agent.response_function.run(game=self, agent=agent, fitted=kwargs.get('fitted', False))
			outcome.append(agent_decision.id)
			# Only apply the last move if this is not fitted
			if (not kwargs.get('fitted', False)):
				agent.last_move = agent_decision
			else:
				agent.fitted_last_move = agent_decision
		return outcome

	def relaxed_parameter_traversal(self, arg_starting_index=0):

		# If there is an upper and lower bound (there should always be one - because we fill it in if it isn't there)
		if len(self.agent_parameter_settings) > 0:
			this_agent_parameter_settings = self.agent_parameter_settings[arg_starting_index]
			lower_bound = int(this_agent_parameter_settings[1]*self.DEFAULT_PRECISION_IN_GRID_SEARCH)
			upper_bound = int(this_agent_parameter_settings[2]*self.DEFAULT_PRECISION_IN_GRID_SEARCH)
			gs_interval = int(this_agent_parameter_settings[3]*self.DEFAULT_PRECISION_IN_GRID_SEARCH)
			adjusted_interval = int((upper_bound-lower_bound)/float(gs_interval)*self.DEFAULT_PRECISION_IN_GRID_SEARCH)
			grid_search_range = [x for x in range(lower_bound, upper_bound, adjusted_interval)]
			# Add the upper bound if it is not in the grid search range
			if grid_search_range[-1] != upper_bound:
				if (abs(grid_search_range[-1]-upper_bound) < adjusted_interval):
					grid_search_range[-1] = upper_bound 
				else:
					grid_search_range.append(upper_bound)
			'''
			if (grid_search_range[-1] != upper_bound):
				print(grid_search_range[-1])
				print(upper_bound)
				grid_search_range.append(upper_bound)
				#grid_search_range[-1] = upper_bound
			'''
			for value in grid_search_range:
				acted_value = (value/self.DEFAULT_PRECISION_IN_GRID_SEARCH)
				self.agent_parameter_settings_traversal.update({ this_agent_parameter_settings[0] : acted_value })
				if (arg_starting_index != (len(self.agent_parameter_settings)-1)):
					if (not self.break_fitting):
						self.relaxed_parameter_traversal(arg_starting_index+1)
				else:

					# Preset the agents fitted values if possible
					for i in range(0, len(self.agents)):
						# Wipe the fitted values - as a precaution for mixing values
						self.agents[i].fitted_k = None
						self.agents[i].fitted_strategy_theta = None
						self.agents[i].fitted_strategy_phi = None
						self.agents[i].fitted_strategy_rho = None
						self.agents[i].fitted_strategy_sigma = None
						if (("k_%s" % (i)) in self.agent_parameter_settings_traversal):
							self.agents[i].fitted_k = self.agent_parameter_settings_traversal[("k_%s" % (i))]
						#
						if (("theta_%s" % (i)) in self.agent_parameter_settings_traversal):
							self.agents[i].fitted_strategy_theta = self.agent_parameter_settings_traversal[("theta_%s" % (i))]
						#
						if (("phi_%s" % (i)) in self.agent_parameter_settings_traversal):
							self.agents[i].fitted_strategy_phi = self.agent_parameter_settings_traversal[("phi_%s" % (i))]
						#
						if (("rho_%s" % (i)) in self.agent_parameter_settings_traversal):
							self.agents[i].fitted_strategy_rho = self.agent_parameter_settings_traversal[("rho_%s" % (i))]
						#
						if (("sigma_%s" % (i)) in self.agent_parameter_settings_traversal):
							self.agents[i].fitted_strategy_sigma = self.agent_parameter_settings_traversal[("sigma_%s" % (i))]

					# If there are no fitted values, then set them
					for agent in self.agents:
						if (agent.fitted_k == None):
							agent.fitted_k = agent.k
						if (agent.fitted_strategy_phi == None):
							agent.fitted_strategy_phi = agent.strategy_phi
						if (agent.fitted_strategy_theta == None):
							agent.fitted_strategy_theta = agent.strategy_theta
						if (agent.fitted_strategy_rho == None):
							agent.fitted_strategy_rho = agent.strategy_rho
						if (agent.fitted_strategy_sigma == None):
							agent.fitted_strategy_sigma = agent.strategy_sigma

					# For a fitted round - this does not edit game content
					game_state_temporal = self.quantum_state_update(
						self.state, 
						self.agents[0].fitted_strategy_theta,
						self.agents[0].fitted_strategy_phi,
						self.agents[0].fitted_strategy_rho,
						self.agents[0].fitted_strategy_sigma,
						self.agents[1].fitted_strategy_theta,
						self.agents[1].fitted_strategy_phi,
						self.agents[1].fitted_strategy_rho,
						self.agents[1].fitted_strategy_sigma)
					#
					projection_map = self.kraus_projections(self.agents[0].fitted_k, self.agents[1].fitted_k, self.options)
					#
					payoffs_map, payoff_identities = self.compile_payoffs_map(self.payoffs)
					#
					outcome = self.determine_outcome(self.agents, payoffs_map, payoff_identities, projection_map, fitted=True, game_state=game_state_temporal)
					#
					passed_fitted_conditions = True
					for fitted_condition in self.fitted_conditions:
						if not eval(fitted_condition):
							passed_fitted_conditions = False
					feasibility = 0
					if passed_fitted_conditions:
						#self.break_fitting = True
						#print("PASSED fitting conditions")
						#for i in range(0, len(self.agents)):
						#	print("   name: %s k: %s theta: %s phi: %s" % (self.agents[i].name,self.agents[i].k,self.agents[i].strategy_theta,self.agents[i].strategy_phi))
						#print("   OUTCOME: %s:%s %s:%s" % (self.agents[0].name,self.option_from_id(outcome[0]).name, self.agents[1].name,self.option_from_id(outcome[1]).name))
						'''
						print("T | A: k=%s theta=%s phi=%s | B: k=%s theta=%s phi=%s" % 
							(self.agents[0].k, self.agents[0].strategy_theta, self.agents[0].strategy_phi,
								self.agents[1].k, self.agents[1].strategy_theta, self.agents[1].strategy_phi))
						'''
						feasibility = 1
					else:
						#print("DID NOT PASS fitting conditions")
						#for i in range(0, len(self.agents)):
						#	print("   name: %s k: %s theta: %s phi: %s" % (self.agents[i].name,self.agents[i].k,self.agents[i].strategy_theta,self.agents[i].strategy_phi))
						#print("   OUTCOME: %s:%s %s:%s" % (self.agents[0].name,self.option_from_id(outcome[0]).name, self.agents[1].name,self.option_from_id(outcome[1]).name))
						'''
						print("F | A: k=%s theta=%s phi=%s | B: k=%s theta=%s phi=%s" % 
							(self.agents[0].k, self.agents[0].strategy_theta, self.agents[0].strategy_phi,
								self.agents[1].k, self.agents[1].strategy_theta, self.agents[1].strategy_phi))
						'''
					for i in range(0, len(self.agents)):
						self.fitted_agent_thetas[i].append(self.agents[i].fitted_strategy_theta)
						self.fitted_agent_phis[i].append(self.agents[i].fitted_strategy_phi)
						self.fitted_agent_rhos[i].append(self.agents[i].fitted_strategy_rho)
						self.fitted_agent_sigmas[i].append(self.agents[i].fitted_strategy_sigma)
						self.fitted_agent_ks[i].append(self.agents[i].fitted_k)
					self.feasibilities.append(feasibility)
		else:
			for agent in self.agents:
				agent.fitted_k = agent.k
				agent.fitted_strategy_phi = agent.strategy_phi
				agent.fitted_strategy_theta = agent.strategy_theta
				agent.fitted_strategy_rho = agent.strategy_rho
				agent.fitted_strategy_sigma = agent.strategy_sigma

			# For a fitted round - this does not edit game content
			game_state_temporal = self.quantum_state_update(
				self.state, 
				self.agents[0].fitted_strategy_theta,
				self.agents[0].fitted_strategy_phi,
				self.agents[0].fitted_strategy_rho,
				self.agents[0].fitted_strategy_sigma,
				self.agents[1].fitted_strategy_theta,
				self.agents[1].fitted_strategy_phi,
				self.agents[1].fitted_strategy_rho,
				self.agents[1].fitted_strategy_sigma)
			#
			projection_map = self.kraus_projections(self.agents[0].fitted_k, self.agents[1].fitted_k, self.options)
			#
			payoffs_map, payoff_identities = self.compile_payoffs_map(self.payoffs)
			#
			outcome = self.determine_outcome(self.agents, payoffs_map, payoff_identities, projection_map, fitted=True, game_state=game_state_temporal)
			#
			passed_fitted_conditions = True
			for fitted_condition in self.fitted_conditions:
				if not eval(fitted_condition):
					passed_fitted_conditions = False
			feasibility = 0
			if passed_fitted_conditions:
				feasibility = 1
			else:
				pass
			for i in range(0, len(self.agents)):
				self.fitted_agent_thetas[i].append(self.agents[i].fitted_strategy_theta)
				self.fitted_agent_phis[i].append(self.agents[i].fitted_strategy_phi)
				self.fitted_agent_rhos[i].append(self.agents[i].fitted_strategy_rho)
				self.fitted_agent_sigmas[i].append(self.agents[i].fitted_strategy_sigma)
				self.fitted_agent_ks[i].append(self.agents[i].fitted_k)
			self.feasibilities.append(feasibility)


	def feasibility(self, agent_a, agent_b):

		# Validate that both of the provided agents are within the game to run the feasibility check
		agent_a_cleared = False
		agent_b_cleared = False
		for agent in self.agents:
			if (agent.id == agent_a.id):
				agent_a_cleared = True
			if (agent.id == agent_b.id):
				agent_b_cleared = True

		if (agent_a_cleared and agent_b_cleared):

			# Determine if there are fitting conditions
			self.fitted_conditions = []
			for i in range(0, len(self.agents)):
				if (self.agents[i].fitted_option != None):
					self.fitted_conditions.append((('self.agents[%s].fitted_last_move.id == self.agents[%s].fitted_option.id') % (i,i)))

			# If there are fitting conditions, this is a fitted round
			if (len(self.fitted_conditions) > 0):

				# Flush the fitting conditions, and declare the temporal conditions
				self.fitted_agent_thetas = []
				self.fitted_agent_phis = []
				self.fitted_agent_rhos = []
				self.fitted_agent_sigmas = []
				self.fitted_agent_ks = []
				self.feasibilities = []
				upper_bound_k = 1.0
				lower_bound_k = 0.0
				upper_bound_theta = mt.pi
				lower_bound_theta = 0.0
				upper_bound_phi = mt.pi
				lower_bound_phi = 0.0
				upper_bound_rho = mt.pi
				lower_bound_rho = 0.0
				upper_bound_sigma = mt.pi
				lower_bound_sigma = 0.0
				self.break_fitting = False

				# Set each of the parameter settings by analysing the accompanying settings as regex
				self.agent_parameter_settings = []
				for i in range(0, len(self.agents)):
					for relaxed_parameter_string in self.agents[i].relaxed_parameters:
						regex_string_parameter_name = r"(?=^).*?(?=\[|$)"
						regex_string_parameter_lower_bound = r"(?<=\[).*(?=->)"
						regex_string_parameter_upper_bound = r"(?<=->).*?(?=\||\])"
						regex_string_parameter_gs_interval = r"(?<=\|).*?(?=\])"
						parameter_name = None
						parameter_lower_bound = 0.0
						parameter_upper_bound = 0.0
						parameter_gs_interval = 1
						for regex_string_parameter in [regex_string_parameter_name, regex_string_parameter_lower_bound, regex_string_parameter_upper_bound, regex_string_parameter_gs_interval]:
							
							# Retrieve the parameter settings from the 'relax' string
							try:
								relaxed_parameter_string_regex_results = [m for m in enumerate(re.finditer(regex_string_parameter, 
									relaxed_parameter_string, re.MULTILINE), start=1)]
								regex_returned_string = relaxed_parameter_string_regex_results[0][1][0]
							except:
								relaxed_parameter_string_regex_results = []
								regex_returned_string = None

							# If the 'relax' string is well formed
							if ((len(relaxed_parameter_string_regex_results) > 0) and (len(regex_returned_string) > 0)):
								if (regex_string_parameter == regex_string_parameter_name):

									parameter_name = regex_returned_string
									parameter_name.strip()
									# Set the default upper and lower bounds, depending on the parameter
									if (parameter_name == 'k'):
										parameter_lower_bound = 0.0
										parameter_upper_bound = 1.0
									if (parameter_name == 'theta'):
										parameter_lower_bound = 0.0
										parameter_upper_bound = mt.pi
									if (parameter_name == 'phi'):
										parameter_lower_bound = 0.0
										parameter_upper_bound = mt.pi
									if (parameter_name == 'rho'):
										parameter_lower_bound = 0.0
										parameter_upper_bound = mt.pi
									if (parameter_name == 'sigma'):
										parameter_lower_bound = 0.0
										parameter_upper_bound = mt.pi

								try:

									# Set the lower bound
									if (regex_string_parameter == regex_string_parameter_lower_bound):
										parameter_lower_bound = float(regex_returned_string)
										parameter_gs_interval = 10 # We set this for cases where there is no step procedure

									# Set the upper bound
									if (regex_string_parameter == regex_string_parameter_upper_bound):
										parameter_upper_bound = float(regex_returned_string)
										parameter_gs_interval = 10 # We set this for cases where there is no step procedure

									# Set the interval between the bounds
									if (regex_string_parameter == regex_string_parameter_gs_interval):
										parameter_gs_interval = int(regex_returned_string)
								except:
									pass

						# Sanity test of parameters
						if (parameter_name != None and (not (parameter_lower_bound == 0.0 and parameter_upper_bound == 0.0)) and parameter_gs_interval > 0):
							
							# Drop the parameter setting into the dictionary of the corresponding agent
							self.agent_parameter_settings.append([ 
								"%s_%s" % (parameter_name,i),  # index of agent included in parameter name
								parameter_lower_bound, 
								parameter_upper_bound, 
								parameter_gs_interval ])

				self.agent_parameter_settings_traversal = {}

				# Sort the parameter settings (This is very important for traversal)
				self.agent_parameter_settings.sort()

				# Set the graphing variables
				for i in range(0, len(self.agents)):
					self.fitted_agent_thetas.append([])
					self.fitted_agent_phis.append([])
					self.fitted_agent_rhos.append([])
					self.fitted_agent_sigmas.append([])
					self.fitted_agent_ks.append([])
				self.feasibilities = []

				self.relaxed_parameter_traversal()
				
				return Feasibility( 
					snapshot_options=self.options, 
					snapshot_payoffs=self.payoffs, 
					snapshot_state_of_game=self.state, 
					agents=self.agents,
					fitted_agent_ks=self.fitted_agent_ks, 
					fitted_agent_thetas=self.fitted_agent_thetas,
					fitted_agent_phis=self.fitted_agent_phis, 
					fitted_agent_rhos=self.fitted_agent_rhos,
					fitted_agent_sigmas=self.fitted_agent_sigmas, 
					agent_parameter_settings=self.agent_parameter_settings, 
					feasibilities=self.feasibilities )
						
				'''
				for i in range(int(lower_bound_k*10), int(upper_bound_k*10)+1):
					for l in range(int(lower_bound_k*10), int(upper_bound_k*10)+1):
						if (not break_fitting):
							self.agents[0].k = (i/10.0)
							self.agents[0].strategise_fitted()#theta=(j/10.0), phi=(k/10.0))
							self.agents[1].k = (l/10.0)
							self.agents[1].strategise_fitted()#theta=(m/10.0), phi=(n/10.0))
							# For a fitted round (this is also editting consistent game content)
							self.combined_strategies_unitary_operator = self.combine_strategies(self.agents)
							self.state = self.combined_strategies_unitary_operator*self.state
							#
							projection_map = self.kraus_projections(self.agents[0].k, self.agents[1].k, self.options)
							#
							payoffs_map, payoff_identities = self.compile_payoffs_map(self.payoffs)
							#
							outcome = self.determine_outcome(self.agents, payoffs_map, payoff_identities, projection_map)
							#
							passed_fitted_conditions = True
							for fitted_condition in fitted_conditions:
								if not eval(fitted_condition):
									passed_fitted_conditions = False
							if passed_fitted_conditions:
								break_fitting = True
				if (break_fitting):
					print('successfully found configuration that fits')
					print('Alice k : %s' % (self.agents[0].k))
					print('Bob k : %s' % (self.agents[1].k))
				else:
					print('did not find configuration that fits')

				# NEED TO ADD PAYOFFS soon
				'''
			else:
				pass
				# TODO Error - no fitted conditions
		else:
			pass
			# TODO Error - agents are not linked to game

	def graph_feasibility(self, feasibility):

		'''
			Feasibility Graph
		'''

		# Graph settings
		feasibility_graph = plt.figure()
		feasibility_graph.subplots_adjust(hspace=1)
		feasibility_graph.set_dpi(80)
		gs = GridSpec((7*len(feasibility.agent_parameter_settings))+5,1)
		green_palette = feasibility_graph.add_subplot(gs[2:4,:])
		plt.bar([x for x in range(0,len(feasibility.feasibilities))],feasibility.feasibilities, 
			label=("Feasibilities"), 
			color=display_theme()["dark"]["color_bank"][4])
		plt.margins(0)
		title = plt.title('Feasibility Of %s Choosing %s And %s Choosing %s' % 
			(feasibility.agent_names[0], feasibility.agent_fitted_option_names[0], feasibility.agent_names[1], feasibility.agent_fitted_option_names[1]), 
			fontsize=display_theme()["dark"]["font_size_title"], 
			family=display_theme()["dark"]["font_family"], 
			color=display_theme()["dark"]["font_color"])
		title.set_position([.5, 1.7])
		feasibility_graph.set_size_inches(11, 5*(((2+((3*len(feasibility.agent_parameter_settings))))/8)))
		
		# Green palette
		green_palette.set_facecolor(display_theme()["dark"]["face_color_backdrop"])
		green_palette.set_yticklabels(['0.00','0.00','0.00'])
		for tick in green_palette.yaxis.get_major_ticks():
			tick.label.set_fontsize(display_theme()["dark"]["font_size_ticks"])
		green_palette.tick_params(axis='x', colors=(0,0,0,0))
		green_palette.tick_params(axis='y', colors=(0,0,0,0))
		green_palette.spines['top'].set_visible(False)
		green_palette.spines['right'].set_visible(False)
		green_palette.spines['bottom'].set_visible(False)
		green_palette.spines['left'].set_visible(False)
		colors = [display_theme()["dark"]["color_bank"][4]]
		texts = ["Feasibility"]
		patches = [ plt.plot([], marker="o", ms=12, ls="", mec=None, color=colors[i], label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
		l = green_palette.legend(loc='center left', bbox_to_anchor=(1, 0.5), 
			facecolor=None, edgecolor=None, framealpha=0, prop={
				"size":display_theme()["dark"]["font_size_legend"], 
				"family":display_theme()["dark"]["font_family"]}, handles=patches, handletextpad=0)
		for text in l.get_texts():
		    text.set_color(display_theme()["dark"]["font_color"])

		# For each tested parameters, draw the accompanying feasibility graph
		for h in range(0, len(feasibility.agent_parameter_settings)):
			plot_data = None
			plot_data_agent_name = None
			plot_data_variable = None
			for i in range(0, len(feasibility.agent_names)):
				term = feasibility.agent_parameter_settings[h][0]
				if (("k_%s" % (i)) in term):
					plot_data = feasibility.fitted_agent_ks[i]
					plot_data_variable = "K"
					plot_data_agent_name = feasibility.agent_names[i]
				#
				if (("theta_%s" % (i)) in term):
					plot_data = feasibility.fitted_agent_thetas[i]
					plot_data_variable = "Theta"
					plot_data_agent_name = feasibility.agent_names[i]
				#
				if (("phi_%s" % (i)) in term):
					plot_data = feasibility.fitted_agent_phis[i]
					plot_data_variable = "Phi"
					plot_data_agent_name = feasibility.agent_names[i]
				#
				if (("rho_%s" % (i)) in term):
					plot_data = feasibility.fitted_agent_rhos[i]
					plot_data_variable = "Rho"
					plot_data_agent_name = feasibility.agent_names[i]
				#
				if (("sigma_%s" % (i)) in term):
					plot_data = feasibility.fitted_agent_sigmas[i]
					plot_data_variable = "Sigma"
					plot_data_agent_name = feasibility.agent_names[i]
			parameter_lower_bound = feasibility.agent_parameter_settings[h][1]
			parameter_upper_bound = feasibility.agent_parameter_settings[h][2]
			parameter_middle_bound = ((parameter_upper_bound+parameter_lower_bound)/2.0)
			display_color = display_theme()["dark"]["color_bank"][(h+5)%len(display_theme()["dark"]["color_bank"])]
			ax1 = feasibility_graph.add_subplot(gs[5+((h)*5):5+((h+1)*5)-1,:])
			ax1.set_ylim(parameter_lower_bound, parameter_upper_bound)
			ax1.set_yticks([parameter_lower_bound, parameter_middle_bound, parameter_upper_bound])
			ax1.set_yticklabels(["%0.2f" % parameter_lower_bound, "%0.2f" % parameter_middle_bound, "%0.2f" % parameter_upper_bound])
			ax1.set_facecolor(display_theme()["dark"]["face_color_backdrop"])
			plt.plot(plot_data, linewidth=1.5, label=("%s's %s" % (plot_data_agent_name, plot_data_variable)), color=display_color)
			plt.bar([x for x in range(0,len(feasibility.feasibilities))],np.array(feasibility.feasibilities)*parameter_upper_bound, 
				label=("Feasibilities"), color=display_theme()["dark"]["color_bank_alpha_5"][4])
			plt.margins(0)
			for tick in ax1.yaxis.get_major_ticks():
				tick.label.set_fontsize(display_theme()["dark"]["font_size_ticks"])
			ax1.tick_params(axis='y', colors=display_theme()["dark"]["font_color_alpha"])
			ax1.spines['top'].set_visible(False)
			ax1.spines['right'].set_visible(False)
			ax1.spines['bottom'].set_visible(False)
			ax1.spines['left'].set_visible(True)
			ax1.spines['left'].set_color(display_theme()["dark"]["font_color_alpha"])
			ax1.tick_params(axis='x', colors=(0,0,0,0))
			ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
			texts = ["%s's %s" % (plot_data_agent_name, plot_data_variable)]
			colors = [ display_color ]
			patches = [ plt.plot([], marker="o", ms=12, ls="", mec=None, color=colors[i], label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
			l = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), facecolor=None, edgecolor=None, 
				framealpha=0, prop={
				"size":display_theme()["dark"]["font_size_legend"], 
				"family":display_theme()["dark"]["font_family"]}, handles=patches, handletextpad=0)
			for text in l.get_texts():
				text.set_color(display_theme()["dark"]["font_color"])
		
		plt.show()

	def summarise_feasibility(self, feasibility):

		is_this_instance_feasible = (1 in feasibility.feasibilities)
		number_of_feasible_occurrences = feasibility.feasibilities.count(1)
		# Summarise the feasibility

		log_title('Feasibility Summary Report')
		log_subtitle("Parametrics")
		log_dictionary({
				"Agents" : ""
			})

		for i in range(0, len(feasibility.agent_names)):
			log_dictionary({
					"%s" % (feasibility.agent_names[i]) : ""
				},indent="\t", new_line=False)
			log_dictionary({
					"K" : feasibility.agents_list[i].k,
					"Strategy - Theta" : feasibility.agents_list[i].strategy_theta,
					"Strategy - Phi" : feasibility.agents_list[i].strategy_phi,
					"Strategy - Rho" : feasibility.agents_list[i].strategy_rho,
					"Strategy - Sigma" : feasibility.agents_list[i].strategy_sigma,
					"Decision Function" : feasibility.agent_response_functions[i]
				},indent="\t\t", new_line=False)

		game_state_ket = np.array(feasibility.snapshot_state_of_game)
		ket_val_1 = game_state_ket[0][0]
		ket_val_2 = game_state_ket[1][0]
		ket_val_3 = game_state_ket[2][0]
		ket_val_4 = game_state_ket[3][0]

		log_dictionary({
				"Quantum State Of Game" : ("Ket ( %0.2f + %0.2fj, %0.2f + %0.2fj, %0.2f + %0.2fj, %0.2f + %0.2fj )" % (ket_val_1.real, ket_val_1.imag, ket_val_2.real, ket_val_2.imag, ket_val_3.real, ket_val_3.imag, ket_val_4.real, ket_val_4.imag)),
				"Entanglement Of State" : determine_entanglement(feasibility.snapshot_state_of_game),
				"Payoffs" : "" # Added for further dictionary
			})
		log_dictionary(feasibility.snapshot_payoffs,indent="\t", new_line=False)

		log_subtitle("Relaxed Parametrics")

		for parameter in feasibility.agent_parameter_settings:
			parameter_code = parameter[0].split("_")[0]
			if (parameter_code == "k"):
				parameter_name = "K"
			if (parameter_code == "theta"):
				parameter_name = "Theta"
			if (parameter_code == "phi"):
				parameter_name = "Phi"
			if (parameter_code == "rho"):
				parameter_name = "Rho"
			if (parameter_code == "sigma"):
				parameter_name = "Sigma"
			parameter_agent_name = feasibility.agent_names[int(parameter[0].split("_")[1])]
			parameter_lower_bound = parameter[1]
			parameter_upper_bound = parameter[2]
			parameter_intervals = parameter[3]
			log_text("%s of %s :" % (parameter_name, parameter_agent_name), indent="\t")
			log_text("Tested range of [%s - %s] at interval of %s" % (parameter_lower_bound, parameter_upper_bound, parameter_intervals), indent="\t\t")

		#TODO

		log_subtitle("Outcome")
		log_dictionary({
				"Requirements" : ""
			})
		log_dictionary({
				"%s's Outcome" % (feasibility.agent_names[0]) : "%s" % (feasibility.agent_fitted_option_names[0]),
				"%s's Outcome" % (feasibility.agent_names[1]) : "%s" % (feasibility.agent_fitted_option_names[1]),
			},indent="\t", new_line=False)
		log_dictionary({
				"Feasible" : str(is_this_instance_feasible),
				"No. Of Tested Configurations" : len(feasibility.feasibilities),
				"No. Of Feasible Configurations" : number_of_feasible_occurrences
			})
		

	def play_round(self):
		# Update the game state
		self.state = self.quantum_state_update(
			self.state, 
			self.agents[0].strategy_theta,
			self.agents[0].strategy_phi,
			self.agents[0].strategy_rho,
			self.agents[0].strategy_sigma,
			self.agents[1].strategy_theta,
			self.agents[1].strategy_phi,
			self.agents[1].strategy_rho,
			self.agents[1].strategy_sigma)

		# Calculate the Kraus projections, and compile the projections dictionary
		projection_map = self.kraus_projections(self.agents[0].k, self.agents[1].k, self.options)

		# Compile the payoffs dictionary
		payoffs_map, payoff_identities = self.compile_payoffs_map(self.payoffs)

		# Use the agent decision functions to calculate the outcome of the game
		# updates:
		# 	agent expected utility
		# 	agent last move
		outcome = self.determine_outcome(self.agents, payoffs_map, payoff_identities, projection_map, fitted=False)
		
		# Attribute payoffs to players (currently, this is only formalised for two players)
		outcome_objects = []
		for outcome_i in outcome:
			for option in self.options:
				if (option.id == outcome_i):
					outcome_objects.append(option)
		for agent in self.agents:
			this_agent_index = self.agent_index(agent)
			for k,v in self.payoffs.items():
				if (v[:len(self.options)] == outcome_objects):
					agent.pot += v[len(self.options)+this_agent_index]
					agent.trace_pot.append(agent.pot)

		# Trace the entanglement
		self.trace_entanglement.append(determine_entanglement(self.state))
		self.trace_correlation.append(correlative_measure(self.state))
		self.trace_contextuality.append(contextuality_measure(self.state))

	# While this function appears to edit the state of the game, it actually only edits the state of the SUPPLIED game
	def quantum_state_update(self, initial_state, theta_a, phi_a, rho_a, sigma_a, theta_b, phi_b, rho_b, sigma_b):
		ket00_t = ((cm.e**(ii*(phi_a+phi_b)))*mt.cos(theta_a/2)*mt.cos(theta_b/2))+(ii*mt.sin(theta_a/2)*mt.sin(theta_b/2))
		ket11_t = (mt.sin(theta_a/2)*mt.sin(theta_b/2))+(ii*cm.e**(-ii*(phi_a+phi_b))*mt.cos(theta_a/2)*mt.cos(theta_b/2))
		#ket01_t = ((-cm.e**(ii*rho_a))*mt.cos(sigma_a/2)*mt.sin(sigma_b/2)) + ((ii*cm.e**(-ii*rho_b))*mt.sin(sigma_a/2)*mt.cos(sigma_b/2))
		#ket10_t = ((-cm.e**(ii*rho_b))*mt.sin(sigma_a/2)*mt.cos(sigma_b/2)) + ((ii*cm.e**(-ii*rho_b))*mt.cos(sigma_a/2)*mt.sin(sigma_b/2))
		ket10_t = ((cm.e**(ii*(rho_a+rho_b)))*mt.cos(sigma_a/2)*mt.cos(sigma_b/2))+(ii*mt.sin(sigma_a/2)*mt.sin(sigma_b/2))
		ket01_t = (mt.sin(sigma_a/2)*mt.sin(sigma_b/2))+(ii*cm.e**(-ii*(rho_a+rho_b))*mt.cos(sigma_a/2)*mt.cos(sigma_b/2))
		return 	((( (ket00_t*Qobj(np.array(tensor(ket("00"),bra("00"))))*Qobj(np.array(initial_state)))+
					(ket11_t*Qobj(np.array(tensor(ket("11"),bra("11"))))*Qobj(np.array(initial_state))) )+
				( (ket10_t*Qobj(np.array(tensor(ket("01"),bra("01"))))*Qobj(np.array(initial_state)))+
					(ket01_t*Qobj(np.array(tensor(ket("10"),bra("10"))))*Qobj(np.array(initial_state))) ))/1)

	def plot_payoffs(self):
		fig, ax = plt.subplots()
		min_v = mt.inf
		max_v = -mt.inf
		for agent in self.agents:
			plt.plot(agent.trace_pot, linewidth=1,label=("%s's Payoffs" % (agent.name)))
			this_max = np.array(agent.trace_pot).max()
			this_min = np.array(agent.trace_pot).min()
			if (this_max > max_v):
				max_v = this_max
			if (this_min > min_v):
				min_v = this_min

		#trace_entanglement_plot = np.array(self.trace_entanglement)/np.array(self.trace_entanglement).max()
		#plt.plot((trace_entanglement_plot*this_max)+this_min, linewidth=1,label="Entanglement")

		#trace_correlation_plot = np.array(self.trace_correlation)/np.array(self.trace_correlation).max()
		#plt.plot((trace_correlation_plot*this_max)+this_min, linewidth=1,label="Correlation")

		
		plt.title('Cumulative Payoffs')
		plt.ylabel("Payoff")
		plt.xlabel("Round")
		plt.margins(.1)
		fig.set_size_inches(10, 4)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		return plt.show()





class Bistable_Quantum_Agent(Agent):
	def __init__(self, **kwargs):
		super().__init__(**kwargs)
		self.id = generate_string()
		self.k = kwargs.get('k', 1) # Implies rationality
		self.strategy_theta = kwargs.get('theta', 0) 
		self.strategy_phi = kwargs.get('phi', 0)
		self.strategy_rho = kwargs.get('rho', 0) 
		self.strategy_sigma = kwargs.get('sigma', 0)
		self.fitted_option = None
		self.fitted_k = None
		self.fitted_last_move = None
		self.fitted_strategy_theta = None
		self.fitted_strategy_phi = None
		self.fitted_strategy_rho = None
		self.fitted_strategy_sigma = None
		self.relaxed_parameters = []

	def strategise(self, **kwargs):
		self.strategy_theta = kwargs.get('theta', 0) 
		self.strategy_phi = kwargs.get('phi', 0)
		self.strategy_rho = kwargs.get('rho', 0) 
		self.strategy_sigma = kwargs.get('sigma', 0)

	def fit(self, **kwargs):
		self.fitted_option = kwargs.get('option', None)
		self.relaxed_parameters = kwargs.get('relax', [])
		return self











