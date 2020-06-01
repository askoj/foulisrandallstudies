# Init

from qogpsi.misc import *
import math as mt
import re
import collections
import itertools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
from qogpsi.networkx_mod.networkx_mod.networkx_mod import *
import picos as pic
import copy
import copy
from itertools import chain
import os
import inspect

import qogpsi.networkx_mod as nx
import matplotlib.colors
from qogpsi.stylization.themes import display_theme
from matplotlib.axes._axes import _log as matplotlib_axes_logger
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
matplotlib_axes_logger.setLevel('ERROR')
def nearest_square(limit):
	answer = 0
	while (answer+1)**2 < limit:
		answer += 1
	return (answer+1)**2


class random_variable(object):
	def __init__(self, arg_name, **kwargs):
		self.id = generate_string()
		self.name = arg_name
		self.states = kwargs.get('states',[])
		self.label = kwargs.get('label',self.name)

class state(object):
	def __init__(self, arg_name, **kwargs):
		self.id = generate_string()
		self.name = arg_name
		self.label = kwargs.get('label',self.name)


def flatten_lists_only(l):
	flat_list = []
	for sublist in l:
		if (type(sublist) == type([])):
			for item in sublist:
				flat_list.append(item)
		else:
			flat_list.append(sublist)
	return flat_list


def c( arg_x, arg_y, state_size, arg_x_n):
	x_aug = arg_x_n+1-arg_x
	y_aug = arg_y
	step = int(np.prod([x for x in state_size[:x_aug-1]]))
	return (mt.floor((y_aug-1)/step)%state_size[x_aug-1])

def procure_axis(axis_rvs):
	axis_rv_state_defs = []
	# These lists are deliberately of length 2 - for both the names and states of the RVs
	axis_annotations = [[],[]] 
	for rv_from in axis_rvs:
		starting = True
		for state_from_rv in rv_from.states:
			axis_rv_state_defs.append([rv_from, state_from_rv])
			if (starting):
				axis_annotations[0].append(rv_from.label)
			else:
				axis_annotations[0].append("")
			axis_annotations[1].append(state_from_rv.label)
			starting = False
	return axis_rv_state_defs, axis_annotations	
class distribution(object):
	def __init__(self, **kwargs):
		self.id = "*" # Need to edit this
		self.distribution_values =  kwargs.get('dv',None)
		self.x_rvs = kwargs.get('x',None)
		self.y_rvs = kwargs.get('y',None)
	def graph(self, **kwargs):
		#print("henlo")
		#print(self.x_rvs)
		highlight = kwargs.get('highlight', None)
		length_of_highlight = None
		if "hyperedges" in highlight:
			length_of_highlight = nearest_square(len(highlight["hyperedges"]))
			highlight = highlight["hyperedges"]
		else:
			length_of_highlight = nearest_square(len(highlight))
		
		scale_all = kwargs.get('scale', 1)

		x_axis_rv_state_defs, x_axis_annotations = procure_axis(self.x_rvs)
		y_axis_rv_state_defs, y_axis_annotations = procure_axis(self.y_rvs)
		#print(y_axis_annotations)
		#print(y_axis_rv_state_defs)

		cartesian_product = []
		cartesian_product_annotations = []
		for i in range(len(x_axis_rv_state_defs)):
			cartesian_product.append([])
			cartesian_product_annotations.append([])
			for j in range(len(y_axis_rv_state_defs)):
				cartesian_product[i].append([x_axis_rv_state_defs[i], y_axis_rv_state_defs[j]])
				#cartesian_product_annotations[i].append()
				for val in self.distribution_values:
					is_equal = True
					#print(val)
					for rv_w_s1 in val["combination"]:
						found = False
						for rv_w_s2 in cartesian_product[i][j]:
							if rv_w_s1 == rv_w_s2:
								found = True
						if (not found):
							is_equal = False
					if (is_equal):
						cartesian_product_annotations[i].append("%.2f" % val["p"])
		#print(cartesian_product_annotations)


		x_axis_length = len(x_axis_rv_state_defs) + 2 # We add 2 to compliment the length of the orth axis
		y_axis_length = len(y_axis_rv_state_defs) + 2
		#print(x_axis_annotations)
		table_values = y_axis_annotations
		table_values.extend(cartesian_product_annotations)

		for i in range(2):
			table_values[i].insert(0,"")
			table_values[i].insert(0,"")
		for i in range(2, len(table_values)):
			for j in range(0, 2):
				table_values[i].insert(j,x_axis_annotations[j][i-2])
		#print(table_values)
		#sys.exit()
		
		dt_type = "dark"
		def set_table_properties(table, iii):
			props = table.properties()
			#print(props)
			cells = props['children']
			k = 0
			for i in range(0, x_axis_length):
				for j in range(0, y_axis_length):
					cells[k].set_facecolor((0.0,0.0,0.00,0))

					if ((i >= 2) and (j >= 2)):
						this_coordinates_object = cartesian_product[i-2][j-2]
						#print(this_coordinates_object)

						#for edge in highlight:
						for coordinate in highlight[iii]:
							if (len(coordinate) == len(this_coordinates_object)):
								found_tally = 0
								for sub_coordinate_1 in coordinate:
									for sub_coordinate_2 in this_coordinates_object:
										#print("sub_coordinate_1",sub_coordinate_1)
										#print("sub_coordinate_2",sub_coordinate_2)
										#sys.exit()
										if (sub_coordinate_1[0][0] == sub_coordinate_2[0]) and (sub_coordinate_1[1][0] == sub_coordinate_2[1]):
											found_tally += 1
								if (found_tally == len(coordinate)):
									cells[k].set_facecolor(display_theme()[dt_type]["color_bank_alpha_25"][5])
										#print("yes")


					#sys.exit()
					#if (highlight["validated_hyperedges"][0]):
					#	cartes

					cells[k].set_edgecolor(display_theme()[dt_type]["line_color_alpha"])
					cells[k].set_height(1./len(table_values))
					cells[k].set_width(1./len(table_values[0]))

					if (j < 2):
						cells[k].set_width(0.4/len(table_values[0]))
						cells[k].get_text().set_rotation(90)
						cells[k].set_linewidth(0)
					if (i < 2):
						cells[k].set_height(0.4/len(table_values))
						cells[k].set_linewidth(0)
					#
					if ((i >= 2) and (j >= 2)):
						cells[k].set_text_props(c=display_theme()[dt_type]["line_color_alpha"],fontproperties=FontProperties(size='24',weight='normal',style='italic'))
					elif (((i < 2) or (j < 2)) 
						and (not ((i < 2) and (j < 2)))):
						if (j >= 2):
							if (i / 2.0 % 1 != 0):
								cells[k].set_text_props(c=display_theme()[dt_type]["line_color_alpha"],fontproperties=FontProperties(size='24',weight='normal',style='italic'))
							else:
								cells[k].set_text_props(c=display_theme()[dt_type]["line_color_alpha"],fontproperties=FontProperties(size='24',weight='bold'))
						else:
							if (j / 2.0 % 1 != 0):
								cells[k].set_text_props(c=display_theme()[dt_type]["line_color_alpha"],fontproperties=FontProperties(size='24',weight='normal',style='italic'))
							else:
								cells[k].set_text_props(c=display_theme()[dt_type]["line_color_alpha"],fontproperties=FontProperties(size='24',weight='bold'))
					k += 1
		highlight_sqrt = int(mt.sqrt(length_of_highlight))
		highlight_len_x = highlight_sqrt
		highlight_len_y = highlight_sqrt

		if (length_of_highlight-len(highlight)) >= highlight_len_y:
			highlight_len_y -= 1
		fig = plt.figure(constrained_layout=True)
		fig.set_size_inches(highlight_len_x*4.5*scale_all, highlight_len_y*4.2*scale_all)
		plt.axis('off')
		spec = gridspec.GridSpec(ncols=highlight_len_x, nrows=highlight_len_y, figure=fig)
		fig_ax = []
		for i in range(highlight_len_y):
			for j in range(highlight_len_x):
				fig_ax.append(fig.add_subplot(spec[i, j]))

		for i in range(len(fig_ax)):
			ax_table = fig_ax[i].table(cellText=table_values, loc='center', cellLoc='center')
			if len(highlight) > i:
				set_table_properties(ax_table, i)
				fig_ax[i].axis("off")
				fig_ax[i].set_facecolor((0.0,0.0,0.00,0))

		#fig = plt.figure( figsize=(len(table_values[0])*1*scale_all, len(table_values)*1*scale_all), dpi=100) # EDIT SIZE HERE
		
		#ax1_table = plt.table(cellText=table_values, loc='center', cellLoc='center')
		#set_table_properties(ax1_table)
		plt.show()
		
		'''

		for event_to_find in self.events_register:
			for i in range(len(temp_distribution_values)):
				equal = False
				if len(event_to_find["combination"]) == len(temp_distribution_values[i]["combination"]):
					# Lengths match, proceed...
					equal = True
					for rv_w_s1 in event_to_find["combination"]:
						found = False
						for rv_w_s2 in temp_distribution_values[i]["combination"]:
							if rv_w_s1 == rv_w_s2:
								found = True
						if (not found):
							equal = False
				if (equal):
					temp_distribution_values[i]["p"] = event_to_find["p"]

		'''
		#print(self.distribution_values)
		'''
		x_axis = self.x_rvs
		y_axis = self.y_rvs
		highlight = kwargs.get('highlight',None)
		def axis_length(said_axis,orth_axis):
			axis_length = len(said_axis[0].states)
			for i in range(1,len(said_axis)):
				axis_length = axis_length*len(said_axis[i].states)
			axis_length_reduced = len(orth_axis)*2
			axis_length += len(orth_axis)*2
			return axis_length, axis_length_reduced
		x_axis_length, x_axis_length_reduced = axis_length(x_axis,y_axis)
		y_axis_length, y_axis_length_reduced = axis_length(y_axis,x_axis)
		table_template = []
		def table_part_appendor( value_to_append, m, n, axis_length_reduced, axis, orth_axis, true_x):
			adding_dimension_reference = False
			adding_dimension = None
			axis_correspondent = None
			value_to_store = None
			value_to_store_axis_rep = None
			if (m >= len(orth_axis)*2):
				if (int(n/2) < len(axis)):
					if ((m == (len(orth_axis)*2)) and (n / 2.0 % 1 == 0)):
						# Set the titles of each variable in the axes
						value_to_append = axis[int(n/2)].label
					if (n / 2.0 % 1 != 0):
						adding_dimension_reference = True
						axis_correspondent = ((m-(axis_length_reduced))+1)
						if (axis == true_x):
							adding_dimension = "x"
						else:
							adding_dimension = "y"
						state_sizes = [len(rv.states) for rv in axis]
						state_to_retrieve = c((int(n/2)+1), axis_correspondent, state_sizes, len(axis))
						value_to_store = axis[int(n/2)].states[state_to_retrieve]
						value_to_store_axis_rep = axis[int(n/2)]
						value_to_append = value_to_store.label
			return value_to_append, adding_dimension_reference, adding_dimension, axis_correspondent, value_to_store, value_to_store_axis_rep

		cartesian_map = {}
		cartesian_random_var_layout = []
		for rv in self.distribution_values["random_variables"]:
			cartesian_random_var_layout.append(rv.name)
		table_verification_model = []
		for i in range(0, x_axis_length):
			table_template.append([])
			table_verification_model.append([])
			for j in range(0, y_axis_length):
				value_to_append_verification = []
				value_to_append = ""
				if (((i < len(y_axis)*2) or (j < len(x_axis)*2)) 
					and (not ((i < len(y_axis)*2) and (j < len(x_axis)*2)))):
					value_to_append, adding_val_x, dim_x, axis_x_correspondent, axis_x_value, axis_x_rep = table_part_appendor( value_to_append, j, i,
						y_axis_length_reduced, y_axis, x_axis, x_axis)
					value_to_append, adding_val_y, dim_y, axis_y_correspondent, axis_y_value, axis_y_rep = table_part_appendor( value_to_append, i, j,
						x_axis_length_reduced, x_axis, y_axis, x_axis)
					if (adding_val_x):
						cartesian_ref_x = ("%s%s" % (dim_x, axis_x_correspondent))
						if not cartesian_ref_x in cartesian_map:
							cartesian_map.update({cartesian_ref_x:{}})
						cartesian_map[cartesian_ref_x].update({axis_x_rep.name:axis_x_value})
					if (adding_val_y):
						cartesian_ref_y = ("%s%s" % (dim_y, axis_y_correspondent))
						if not cartesian_ref_y in cartesian_map:
							cartesian_map.update({cartesian_ref_y:{}})
						cartesian_map[cartesian_ref_y].update({axis_y_rep.name:axis_y_value})
				if (((i >= len(y_axis)*2) and (j >= len(x_axis)*2))):
					x_rep = cartesian_map["x%s" % (i-(len(y_axis)*2)+1)]
					y_rep = cartesian_map["y%s" % (j-(len(x_axis)*2)+1)]
					constructed_product_rep = [None for i in range(len(cartesian_random_var_layout))]
					for q in range(0, len(cartesian_random_var_layout)):
						rv = cartesian_random_var_layout[q]
						if rv in x_rep:
							constructed_product_rep[q] = x_rep[rv]
						if rv in y_rep:
							constructed_product_rep[q] = y_rep[rv]
					for d_val in self.distribution_values["outcomes"]:
						if (d_val[0] == constructed_product_rep):
							value_to_append = str(d_val[1])
							value_to_append_verification = d_val[0]
				table_template[i].append(value_to_append)
				table_verification_model[i].append(value_to_append_verification)
		table_values = table_template
		print(table_verification_model)
		dt_type = "dark"
		fig = plt.figure( figsize=(len(table_values[0])*1, len(table_values)*1), dpi=100)
		def set_table_properties( table):
			props = table.properties()
			cells = props['child_artists']
			k = 0
			for i in range(0, x_axis_length):
				for j in range(0, y_axis_length):
					cells[k].set_facecolor((0.0,0.0,0.00,0))

					# for setting the face color
					identity_ref_states = table_verification_model[i][j]
					if highlight != None:
						if len(identity_ref_states) > 0:
							identity_ref_object = []
							for i in range(0, len(identity_ref_states)):
								identity_ref_object.append([
									self.distribution_values["random_variables"][i],
									identity_ref_states[i]])
							print(identity_ref_object)

					cells[k].set_edgecolor(display_theme()[dt_type]["line_color_alpha"])
					cells[k].set_height(1./len(table_values))
					cells[k].set_width(1./len(table_values[0]))
					if (j < len(x_axis)*2):
						cells[k].set_width(0.4/len(table_values[0]))
						cells[k].get_text().set_rotation(90)
						cells[k].set_linewidth(0)
					if (i < len(y_axis)*2):
						cells[k].set_height(0.4/len(table_values))
						cells[k].set_linewidth(0)
					#
					if (((i >= len(y_axis)*2) and (j >= len(x_axis)*2))):
						cells[k].set_text_props(c=display_theme()[dt_type]["line_color_alpha"],fontproperties=FontProperties(size='24',weight='normal',style='italic'))
					elif (((i < len(y_axis)*2) or (j < len(x_axis)*2)) 
						and (not ((i < len(y_axis)*2) and (j < len(x_axis)*2)))):
						if (j >= len(x_axis)*2):
							if (i / 2.0 % 1 != 0):
								cells[k].set_text_props(c=display_theme()[dt_type]["line_color_alpha"],fontproperties=FontProperties(size='24',weight='normal',style='italic'))
							else:
								cells[k].set_text_props(c=display_theme()[dt_type]["line_color_alpha"],fontproperties=FontProperties(size='24',weight='bold'))
						else:
							if (j / 2.0 % 1 != 0):
								cells[k].set_text_props(c=display_theme()[dt_type]["line_color_alpha"],fontproperties=FontProperties(size='24',weight='normal',style='italic'))
							else:
								cells[k].set_text_props(c=display_theme()[dt_type]["line_color_alpha"],fontproperties=FontProperties(size='24',weight='bold'))
					k += 1
		plt.axis('off')
		ax1_table = plt.table(cellText=table_values, loc='center', cellLoc='center')
		set_table_properties(ax1_table)
		plt.show()
		'''
		
		

class probability_space(object):
	def __init__(self, **kwargs):
		self.id = "*" # Need to edit this
		# For when the x or y axis are unspecified
		dummy_state = state("")
		dummy_rv = random_variable("", states=[dummy_state])
		self.random_variables = kwargs.get('random_vars',[])
		self.random_variables_dict = {}
		for rv in self.random_variables:
			self.random_variables_dict[rv.name] = rv
		self.events_register = []
		self.nominated_x_rvs = None
		self.nominated_y_rvs = None
		self.states = kwargs.get('states',[])
		self.states_dict = {}
		for s in self.states:
			self.states_dict[s.name] = s

	def rv_from_id(id):
		for rv in self.random_variables:
			if rv.id == id:
				return rv
		return None

	def state_from_id(id):
		for st in self.states:
			if st.id == id:
				return st
		return None
	def event(self, condition, **kwargs):
		random_variable_names = condition.split("|")
		random_variable_identities = []
		for rv_name in random_variable_names:
			random_variable_identities.append(self.random_variables_dict[rv_name])
		sub_events = kwargs.get('p',None)
		if sub_events != None:
			for k, v in sub_events.items():
				state_names = k.split("|")
				state_identities = []
				for state_name in state_names:
					state_identities.append(self.states_dict[state_name])
				new_event = {"combination":[], "p":v}
				for i in range(len(random_variable_identities)):
					rv = random_variable_identities[i]
					st = state_identities[i]
					new_event["combination"].append([rv,st])
				self.events_register.append(new_event)
		return self.events_register
	'''
	def event(self, condition, sub_events=None, **kwargs):
		handled_probability = kwargs.get('p',None)
		if (handled_probability != None):
			return { "condition" : condition, "p" : handled_probability }
		else:
			if (sub_events != None) and (handled_probability == None):
				return_list = []
				sub_events_flattened = flatten_lists_only(sub_events)
				for sub_event in sub_events_flattened:
					gather_string = "%s|%s" % (condition, sub_event["condition"])
					removal_list = [m for m in enumerate(re.finditer(r"(?=\=).*?(\||$)", gather_string, re.MULTILINE), start=1)]
					for rl_str in removal_list:
						gather_string = gather_string.replace(rl_str[1][0], ',')
					gather_string_members = gather_string.split(",")
					gather_string_members_improved = []
					for i in range(0, len(gather_string_members)):
						if len(gather_string_members[i]) > 0:
							gather_string_members_improved.append("self.random_variables_dict['%s']," % gather_string_members[i])
					gather_string_members_improved_string = "".join(gather_string_members_improved)
					gather_list = eval(gather_string_members_improved_string)
					if set(gather_list) == set(self.random_variables):
						self.events_register.append({"%s|%s" % (condition, sub_event["condition"]): sub_event["p"]})
					else:
						return_list.append({ "condition" : "%s|%s" % (condition, sub_event["condition"]), "p" : sub_event["p"] })
				return return_list
	'''
	def distribution(self, **kwargs):
		self.nominated_x_rvs = kwargs.get('x',[])
		self.nominated_y_rvs = kwargs.get('y',[])
		temp_distribution_values = []
		for i in range(len(self.nominated_x_rvs)):
			for j in range(len(self.nominated_y_rvs)):
				for i_i in range(len(self.nominated_x_rvs[i].states)):
					for j_j in range(len(self.nominated_y_rvs[j].states)):
						temp_distribution_values.append({"combination":[
							[self.nominated_x_rvs[i], self.nominated_x_rvs[i].states[i_i]],
							[self.nominated_y_rvs[j], self.nominated_y_rvs[j].states[j_j]]
							], "p":0.0})

		# Matching process - match combinations with event probabilities

		for event_to_find in self.events_register:
			for i in range(len(temp_distribution_values)):
				equal = False
				if len(event_to_find["combination"]) == len(temp_distribution_values[i]["combination"]):
					# Lengths match, proceed...
					equal = True
					for rv_w_s1 in event_to_find["combination"]:
						found = False
						for rv_w_s2 in temp_distribution_values[i]["combination"]:
							if rv_w_s1 == rv_w_s2:
								found = True
						if (not found):
							equal = False
				if (equal):
					temp_distribution_values[i]["p"] = event_to_find["p"]

		return distribution(
			dv=temp_distribution_values,
			x=self.nominated_x_rvs,
			y=self.nominated_y_rvs
		)
	'''
	def distribution(self, **kwargs):
		self.nominated_x_rvs = kwargs.get('x',[])
		self.nominated_y_rvs = kwargs.get('y',[])
		distribution_values = {
			"random_variables" : self.random_variables, "outcomes" : None }
		prod_of_rvs = int(np.prod([len(rv.states) for rv in distribution_values["random_variables"]]))
		len_of_rvs = len(distribution_values["random_variables"])
		outcomes_list = [[None for i in range(len_of_rvs)] for i in range(prod_of_rvs)]
		for i in range(prod_of_rvs):
			for j in range(len_of_rvs):
				rv_in_q = distribution_values["random_variables"][j]
				states_length = len(rv_in_q.states)
				state_sizes = [len(rv.states) for rv in distribution_values["random_variables"]]
				outcomes_list[i][j] = rv_in_q.states[c(j+1, i+1, state_sizes, len_of_rvs)]
		outcomes_probs_list = []
		for i in range(prod_of_rvs):
			outcomes_probs_list.append([outcomes_list[i],0.0])
		distribution_values["outcomes"] = outcomes_probs_list
		
		# Add values
		ref_list = []
		for event_register_item in self.events_register:
			for k, v in event_register_item.items():
				keys = ([m for m in enumerate(re.finditer(r'((?<=^)|(?<=\|)).*?(?=\=)', k, re.MULTILINE), start=1)])
				vals = ([m for m in enumerate(re.finditer(r'(?<=\=).*?(?=\||$)', k, re.MULTILINE), start=1)])
				test_list = []
				test_list_k = []
				for va in vals:
					test_list.append(eval("self.states_dict['%s']" % va[1][0]))
				for ke in keys:
					test_list_k.append(eval("self.random_variables_dict['%s']," % ke[1][0]))
				ref_list.append([test_list_k, test_list, v])
		for h in range(0, len(ref_list)):
			for i in range(0, len(distribution_values["outcomes"])):
				similar = True
				for k in range(0, len(ref_list[h][0])):
					mapped_index = distribution_values["random_variables"].index(ref_list[h][0][k][0])
					if (not distribution_values["outcomes"][i][0][mapped_index] == ref_list[h][1][k]):
						similar = False
				if (similar):
					distribution_values["outcomes"][i][1] = ref_list[h][2]
		
		# Test integrity
		unique_scenarios = []
		for entry in distribution_values["outcomes"]:
			if not entry[0][:-2] in unique_scenarios:
				unique_scenarios.append(entry[0][:-2])
		unique_scenarios_tallies = [0 for x in range(len(unique_scenarios))]
		for entry in distribution_values["outcomes"]:
			this_index = unique_scenarios.index(entry[0][:-2])
			unique_scenarios_tallies[this_index] += entry[1]
		for val in unique_scenarios_tallies:
			if ((val > 1) or (val < 1)):
				raise ValueError("One or more of the joint distributions in the probability space is malformed.")
		
		return distribution(
			dv=distribution_values,
			x=self.nominated_x_rvs,
			y=self.nominated_y_rvs
		)
	'''
from random import choice, randrange
from string import ascii_lowercase

def generate_string(n=10):
	return "".join(choice(ascii_lowercase) for i in range(n))
def coordinate_to_tuple(coordinate):
	tup = []
	for part in coordinate:
		tup.append(tuple(part))
	return tuple(tup)
def prefix_from_coordinate(direct_map, coordinate):
	for k, v in direct_map.items():
		if (set(coordinate_to_tuple(v)) == set(coordinate_to_tuple(coordinate))):
			return k
def flat(val):
	ll = list(chain.from_iterable(val))
	while (type([]) in [type(l) for l in ll]):
		ll = list(chain.from_iterable(ll))
	return ll

def foulis_randall(this_axes, **kwargs):
	num_rvs_all_axes = [len(ax) for ax in this_axes]
	debug_mode = kwargs.get('debug', False)
	aliases = kwargs.get('aliases', None)
	# The universal_map is the data structure used to organise the G = (V, E)
	universal_map = {}
	flat_var_list = []
	for i in range(len(this_axes)): # axis
		universal_map.update({str(i):{}})
		flat_var_list.append([])
		for j in range(len(this_axes[i])): # rv
			universal_map[str(i)].update({str(j):{}})
			for k in range(len(this_axes[i][j].states)): # state
				universal_map[str(i)][str(j)].update({
					str(k):{ 
						"rv" : [this_axes[i][j]], 
						"state" : [this_axes[i][j].states[k]]
					}})
				flat_var_list[i].append([[this_axes[i][j]], [this_axes[i][j].states[k]]])
	#print("flat_var_list",flat_var_list)
	aliases_map = {}
	if aliases != None:
		k = 0
		for j in range(0,len(aliases)):
			for i in range(0,len(aliases[j])):
				aliases_map.update({ aliases[j][i] : flat_var_list[j][i] })
				k += 1
	#print("aliases_map", aliases_map)
	'''
		Calculate all FR associations
	'''

	def comb(set1, set2, k, first_iteration=False):
		if (first_iteration):
			set2 = set1
		mock_set = set2
		dest_set = {}
		if (first_iteration):
			mock_set = dest_set
		for a in set1:
			for b in set2:
				if (b != a) and (set1[a] not in set2[b]) and (set2[b] not in set1[a]):
					candid_val = [set1[a],set2[b]]
					if (candid_val not in mock_set.values()) and (candid_val[::-1] not in mock_set.values()):
						if (len(set(flat(candid_val))) == len(flat(candid_val))):
							dest_set.update({k:candid_val})
							k += 1
		return dest_set, k
	#
	a_set = {}
	indx = 0
	for party in universal_map:
		a_set.update({ indx : [str(indx)] })
		indx += 1
	#
	b_set = a_set
	k = len(a_set.items())
	end_set = copy.deepcopy(a_set)
	end_set_prev = []
	while (len(end_set) != len(end_set_prev)):
		end_set_prev = copy.deepcopy(end_set)
		b_set, k = comb(a_set, copy.deepcopy(b_set), k, first_iteration=(k==len(a_set.items())))
		end_set.update(b_set)
	end_set_list = []
	for k,v in end_set.items():
		if (not (type([]) in [type(i) for i in v])):
			end_set_list.append(universal_map[str(k)])
		else:
			end_set_list.append(None)
	for q in range(0, len(end_set_list)):
		if (end_set_list[q] == None):
			k = 0
			new_edges = {}
			for a in range(0, 2):
				b = [[0,1],[1,0]]
				Ea = end_set_list[list(end_set.values()).index(end_set[q][b[a][0]])]
				Eb = end_set_list[list(end_set.values()).index(end_set[q][b[a][1]])]
				for edge_k, edge_v in Ea.items():
					vertices_on_Ea = list(edge_v.keys())
					edges_of_Eb = list(Eb.keys())
					combinations_of_Eb = [list(result) for result in (itertools.product(''.join(edges_of_Eb), repeat=len(vertices_on_Ea)))]
					#print("edges_of_Eb",edges_of_Eb)
					for combination in combinations_of_Eb:
						one_edge = {}
						#print(combination)
						#print(vertices_on_Ea)
						j = 0
						for i in range(0, len(vertices_on_Ea)):
							vertex_on_Ea = Ea[edge_k][str(i)]
							edge_on_Eb = Eb[str(combination[i])]
							for vertex_on_Eb_k, vertex_on_Eb_v in edge_on_Eb.items():
								new_composite_vertex = copy.deepcopy(vertex_on_Ea)
								new_composite_vertex['rv'].extend(list(vertex_on_Eb_v['rv']))
								new_composite_vertex['state'].extend(vertex_on_Eb_v['state'])
								one_edge.update({ str(j) : new_composite_vertex })
								j += 1
						new_edges.update({str(k) : one_edge})
						k += 1
			end_set_list[q] = new_edges
	#print("end_set_list",(end_set_list))
	formal_fr_products = []
	#print("end_set",end_set)
	#print("a_set",flat(list(a_set.values())))
	orderings = []
	for k, val in end_set.items():
		if (set(flat(list(a_set.values()))) == set(flat(val))):
			orderings.append(val)
			formal_fr_products.append(end_set_list[list(end_set.values()).index(val)])
	#print("formal_fr_products",formal_fr_products)
	#	
	def get_universal_map_ref(rv, state):
		for axis_k, axis_v in universal_map.items():
			for edge_k, edge_v in universal_map[axis_k].items():
				for vertex_k, vertex_v in universal_map[axis_k][edge_k].items():
					if (rv.id == vertex_v['rv'][0].id and state.id == vertex_v['state'][0].id):
						return [axis_k, edge_k, vertex_k]
		'''
			for axis_k, axis_v in universal_map.items():
			for edge_k, edge_v in universal_map[axis_k].items():
				for vertex_k, vertex_v in universal_map[axis_k][edge_k].items():
					print("rv",vertex_v['rv'][0].id, "state",vertex_v['state'][0].id)
		return [rv.id,state.id]
		'''
		return None
		
	direct_map_new = {}
	refactored_aliased_hyperedges_new = []
	k = 0
	for fr_product in formal_fr_products:
		refactored_aliased_hyperedges_new.append(set())
		for hyperedge_k, hyperedge_v in fr_product.items():
			#refactored_aliased_hyperedges_new[0].append()
			hyperedge_string_set = set()
			for vertex_k, vertex_v in fr_product[hyperedge_k].items():
				name_set = set()
				list_of_refs = []
				for i in range(0, len(fr_product[hyperedge_k][vertex_k]['rv'])):
					list_of_refs.append(get_universal_map_ref(fr_product[hyperedge_k][vertex_k]['rv'][i],
							fr_product[hyperedge_k][vertex_k]['state'][i]))
					name_set.add(fr_product[hyperedge_k][vertex_k]['rv'][i].id+"_"+fr_product[hyperedge_k][vertex_k]['state'][i].id+".")
				alias = ''.join(sorted(name_set))[:-1]
				direct_map_new.update({alias:list_of_refs})
				hyperedge_string_set.add(alias)
			refactored_aliased_hyperedges_new[k].add(':'.join(sorted(hyperedge_string_set)))
		k += 1

	#print("refactored_aliased_hyperedges_new",refactored_aliased_hyperedges_new)

	frp_i = 0
	frp_hyp_i = 0
	frp_hyp_vtx_i = 0
	refactored_formal_hyperedges_new = []
	hyperedges_new = []
	for fr_product in refactored_aliased_hyperedges_new:
		refactored_formal_hyperedges_new.append([])
		hyperedges_new.append([])
		frp_hyp_i = 0
		for hyperedge in fr_product:
			frp_hyp_vtx_i = 0
			refactored_formal_hyperedges_new[frp_i].append([])
			hyperedges_new[frp_i].append([])
			vertex_set = hyperedge.split(":")
			for vertex in vertex_set:
				hyperedges_new[frp_i][frp_hyp_i].append([])
				refactored_formal_hyperedges_new[frp_i][frp_hyp_i].append(direct_map_new[vertex])
				for ref in direct_map_new[vertex]:
					hyperedges_new[frp_i][frp_hyp_i][frp_hyp_vtx_i].append([
						universal_map[ref[0]][ref[1]][ref[2]]['rv'],
						universal_map[ref[0]][ref[1]][ref[2]]['state']])
				frp_hyp_vtx_i += 1
			frp_hyp_i += 1
		frp_i += 1
	#print("refactored_formal_hyperedges_new",refactored_formal_hyperedges_new)
	#print("hyperedges_new",hyperedges_new)
	#print("direct_map_new",direct_map_new)

	#
	formal_fr_products_aliases = []
	alias_i = 0
	for fr_product in formal_fr_products:
		formal_fr_products_aliases.append([])
		edge_i = 0
		for k_edge, v_edge in fr_product.items():
			formal_fr_products_aliases[alias_i].append([])
			vertices_i = 0
			for k_set_of_vertices, v_set_of_vertices in v_edge.items():
				formal_fr_products_aliases[alias_i][edge_i].append([])
				for i in range(0, len(v_set_of_vertices['rv'])):
					alias_retrieved = None
					#print([[v_set_of_vertices['rv'][i].id],[v_set_of_vertices['state'][i].id]])

					for k,v in aliases_map.items():
						#print([v[0][0].id,v[1][0].id])
						if ([[v[0][0].id],[v[1][0].id]] == [[v_set_of_vertices['rv'][i].id],[v_set_of_vertices['state'][i].id]]):
							alias_retrieved = k
					#print(alias_retrieved)
					formal_fr_products_aliases[alias_i][edge_i][vertices_i].append(alias_retrieved)
				formal_fr_products_aliases[alias_i][edge_i][vertices_i] = set(formal_fr_products_aliases[alias_i][edge_i][vertices_i])
				vertices_i += 1
			edge_i += 1
		alias_i += 1
	#print("formal_fr_products_aliases",formal_fr_products_aliases)
	alias_sets_temp = []
	for i in range(len(hyperedges_new)):
		alias_sets_temp.append(set())
		for member in refactored_aliased_hyperedges_new[i]:
			alias_set = []
			for member_part in member.split(":"):
				cood_set = direct_map_new[member_part]
				alias_set_mp = []
				for cood in cood_set:
					for state in universal_map[cood[0]][cood[1]][cood[2]]['state']:
						alias_set_mp.append(state.label)
				alias_set.append("_".join(sorted(alias_set_mp)))
			alias_sets_temp[i].add(":".join(sorted(alias_set)))
	asn = []
	for pr in alias_sets_temp:
		member_reconstructed = []
		for member_flat in pr:
			groups = []
			for group in member_flat.split(":"):
				coods = []
				for cood in group.split("_"):
					coods.append(cood)
				groups.append(sorted(coods))
			member_reconstructed.append(sorted(groups))
		asn.append(sorted(member_reconstructed))
	#print("asn",asn)
	'''
	alias_sets = []
	for i in range(len(alias_sets_temp)):
		#alias_sets.append([])
		for members in alias_sets_temp[i]:
			members_new = []
			for member in members.split(":"):
				coords = []
				for cood in member.split("_"):
					coords.append(cood)
				members_new.append(coords)
		alias_sets.append(members_new)
	print("alias_sets",alias_sets)
	'''
	#print("alias_sets_temp",alias_sets_temp)


	rtn = []
	for i in range(0, len(hyperedges_new)):
		rtn.append({
			"orderings" : orderings,
			"direct_map" : direct_map_new,
			"universal_map" : universal_map,
			"hyperedges" : hyperedges_new[i],
			"refactored_formal_hyperedges" : refactored_formal_hyperedges_new[i],
			"refactored_aliased_hyperedges" : refactored_aliased_hyperedges_new[i],
			"formal_fr_products_aliases" : asn,
			"num_rvs_all_axes" : num_rvs_all_axes
			})
	return rtn
	
	return {
		"direct_map": direct_map_new,
		"universal_map": universal_map,
		"hyperedges" : hyperedges_new,
		"refactored_formal_hyperedges" : refactored_formal_hyperedges_new,
		"refactored_aliased_hyperedges" : refactored_aliased_hyperedges_new
	}
	
	# Supposed to end here
	#sys.exit()

	'''
	FR_associations = [[0,1]]
	for association in FR_associations:
		Ea = universal_map[str(association[0])]
		Eb = universal_map[str(association[1])]
		new_edges = {}
		k = 0
		for edge_k, edge_v in Ea.items():
			vertices_on_Ea = list(edge_v.keys())
			edges_of_Eb = list(Eb.keys())
			combinations_of_Eb = [list(result) for result in (itertools.product(''.join(edges_of_Eb), repeat=len(vertices_on_Ea)))]
			#print("edges_of_Eb",edges_of_Eb)
			for combination in combinations_of_Eb:
				one_edge = {}
				#print(combination)
				#print(vertices_on_Ea)
				j = 0
				for i in range(0, len(vertices_on_Ea)):
					vertex_on_Ea = Ea[edge_k][str(i)]
					edge_on_Eb = Eb[str(combination[i])]
					for vertex_on_Eb_k, vertex_on_Eb_v in edge_on_Eb.items():
						new_composite_vertex = copy.deepcopy(vertex_on_Ea)
						new_composite_vertex['rv'].extend(list(vertex_on_Eb_v['rv']))
						new_composite_vertex['state'].extend(vertex_on_Eb_v['state'])
						one_edge.update({ str(j) : new_composite_vertex })
						j += 1
					#print("vertex_on_Ea",)
					#print("edge_on_Eb", Eb[str(combination[i])])
					print("\n\n\n")
					#new_edges.append([vertices_on_Ea[i],combination[i]])
				new_edges.update({str(k) : one_edge})
				k += 1
		print("new_edges",new_edges)
		sys.exit()
			#print()
			#length_of_ordering = len(edge)
	'''
	
	#print(edge_paths_list)
	# Produce the axes for each of the mutual observers (a data structure that details the random variables in each axis)
	edge_map = {}
	for i in range(0, len(this_axes)):
		edge_map.update({"%s" % (i) : this_axes[i]})
	#print(edge_map)

	# Compute all non-repeating ordered combinations of all axes. This requires a traversal function
	edge_paths_list = []
	traverse_level = len(this_axes)
	def edge_paths(arg_hand_list, arg_traverse_level):
		for i in range(0, len(this_axes)):
			if (not i in arg_hand_list):
				list_to_put = [i]
				list_to_put.extend(arg_hand_list)
				if (arg_traverse_level > 1):
					edge_paths(list_to_put, arg_traverse_level-1)
				else:
					edge_paths_list.append(list_to_put)
	edge_paths([],traverse_level)
	# The edge paths
	#print(edge_paths_list)
	


	informal_hyperedges = []
	for configuration in edge_paths_list:
		if debug_mode:
			print(" ")
			print(" ")
			print("configuration",configuration)
			print(" ")
		start_axis_index = 0
		# Depending on the configuration, we present the first index as the reference point (step 1)
		hyper_edge_vertex_correspondents = []
		for random_variable_index, v in universal_map[str(configuration[start_axis_index])].items():
			hyper_edge_vertex_correspondent = []
			for state_index, v in universal_map[str(configuration[start_axis_index])][str(random_variable_index)].items():
				vertex_id = ([str(configuration[start_axis_index]), random_variable_index, state_index])
				hyper_edge_vertex_correspondent.append(vertex_id)
			hyper_edge_vertex_correspondents.append(hyper_edge_vertex_correspondent)
		#print("{hyperedges from this configuration's first index}")
		#print("hyperedges", hyper_edge_vertex_correspondents)
		# For every orthogonal axis
		#print("	","fb")
		#print("	",hyper_edge_vertex_correspondents)
		if debug_mode:
			print(hyper_edge_vertex_correspondents)
		orth_axis_combinations = []
		# every edge correspondent on the starting axis (step 2)
		for i in range(0, len(hyper_edge_vertex_correspondents)):
			# group combinations by correspondent
			orth_axis_combinations.append([])
			if debug_mode:
				print("computing for h-edge correspondent", i, "of starting axis")
			# for all orthogonal axes (step 3)
			for orth_axis in configuration[1:]:
				# For all RVs in the axis
				#
				if debug_mode:
					print("\treading orth axis", orth_axis)
				orth_axis_rv_combinations = []
				traverse_level = len(hyper_edge_vertex_correspondents[i])
				def traverse_orth_axis(arg_hand_list, arg_traverse_level):
					for h in range(0, len(universal_map[str(orth_axis)].items())):
						list_to_put = [h]
						list_to_put.extend(arg_hand_list)
						if (arg_traverse_level > 1):
							traverse_orth_axis(list_to_put, arg_traverse_level-1)
						else:
							orth_axis_rv_combinations.append(list_to_put)
				traverse_orth_axis([],traverse_level)
				#print("		fc - O%s"%(i))
				#print("		",orth_axis_rv_combinations)
				if debug_mode:
					print("\torth_axis_rv_combinations",orth_axis_rv_combinations)
				# Add this orth axis' combinations to the global list for this correspondent
				orth_axis_combinations[-1].append(orth_axis_rv_combinations)
			#print("")
			if debug_mode:
				print("\torth_axis_isolated_combinations_of_rvs",orth_axis_combinations)
			orth_axis_combined_combinations = []
			traverse_level = len(configuration[1:])
			def traverse_orth_axis_combinations(arg_hand_list, arg_traverse_level):
				current_axis = str(configuration[1:][len(configuration[1:])-arg_traverse_level])
				num_rvs_in_current_axis = len(universal_map[current_axis].items())
				num_vertices_in_hedge_correspondent = len(hyper_edge_vertex_correspondents[i])
				if debug_mode:
					print("\tnum_vertices_in_hedge_correspondent", num_vertices_in_hedge_correspondent)
					print("\tcombinations of rvs", num_rvs_in_current_axis**num_vertices_in_hedge_correspondent)
					print("\tnum_rvs_in_current_axis",num_rvs_in_current_axis)
					print("\tcurrent_axis",current_axis)
				for h in range(0, num_rvs_in_current_axis**num_vertices_in_hedge_correspondent):
						list_to_put = [h]
						list_to_put.extend(arg_hand_list)
						if (arg_traverse_level > 1):
							traverse_orth_axis_combinations(list_to_put, arg_traverse_level-1)
						else:
							list_to_put.reverse()
							orth_axis_combined_combinations.append(list_to_put)
			traverse_orth_axis_combinations([],traverse_level)
			#print("			fc-2")
			#print("			",orth_axis_combined_combinations)

			# Finished interpretation here...

			#print(orth_axis_combined_combinations)
			#sys.exit()
			if debug_mode:
				print("\torth_axis_combined_combinations",orth_axis_combined_combinations)
			# For each of the combined combinations, there is a multidimensional hyperedge that 
			# needs to be retrieved
			nominated_hedge_correspondent = hyper_edge_vertex_correspondents[i]
			if debug_mode:
				print("nominated_hedge_correspondent",nominated_hedge_correspondent)
			for combined_combination in orth_axis_combined_combinations:
				if debug_mode:
					print(" ")
					print("combined_combination",combined_combination)
					print(" ")
				#print(print(orth_axis_combinations[i]))
				for j in range(len(combined_combination)):
					#print(orth_axis_combinations[i][j])
					if debug_mode:
						print("axis", configuration[1:][j])
						print("rvs chosen in axis for hedge", orth_axis_combinations[i][j][combined_combination[j]])
					
			if debug_mode:
				print("break")   
				print("break")  
				print("break")  
			
			
			for combined_combination in orth_axis_combined_combinations:
				if debug_mode:
					print("\t Beginning hyperedge assembly")
				hyperedge_coordinates = []
				for nominated_hyperedge_vertex_i in range(0, len(nominated_hedge_correspondent)):
					traverse_level = len(combined_combination)
					if debug_mode:
						print("\trunning traversal of hyperedge correspondent ", nominated_hyperedge_vertex_i)
						print("\ttesting combined combination:",combined_combination)
					#
					traverse_combined_combination_list = []
					def traverse_combined_combination(arg_hand_list, arg_traverse_level):
						current_axis = str(configuration[1:][traverse_level-arg_traverse_level])
						j = traverse_level-arg_traverse_level
						rvs_of_axis = orth_axis_combinations[i][j][combined_combination[j]]
						selected_rv = rvs_of_axis[nominated_hyperedge_vertex_i]
						num_of_states_in_rv = len(universal_map[str(current_axis)][str(selected_rv)].items())
						if debug_mode:
							print("\t\tnum_of_states_in_rv", num_of_states_in_rv)
							print("\t\tselected_rv", selected_rv)
							print("\t\trvs of axis",rvs_of_axis)
							print("\t\tcurrent_axis", current_axis)
						#sys.exit()
						#print(combined_combination[nominated_hyperedge_vertex_i])
						#rv_of_axis_length = len(universal_map[str(current_axis)][str(combined_combination[nominated_hyperedge_vertex_i])].items())
						for h in range(0, num_of_states_in_rv):
								list_to_put = [[current_axis, str(selected_rv), str(h)]]
								list_to_put.extend(arg_hand_list)
								if (arg_traverse_level > 1):
									traverse_combined_combination(list_to_put, arg_traverse_level-1)
								else:
									list_to_put.reverse()
									traverse_combined_combination_list.append(list_to_put)
					traverse_combined_combination([],traverse_level)
					#print("				fc-3 - T")
					#print("				",traverse_combined_combination_list)
					if debug_mode:
						print("\trecalling that we are operating on hyperedge correspondent vertex: ", nominated_hyperedge_vertex_i)
						print("\tcombinations are as follows:", traverse_combined_combination_list)
					for coordinate_group in traverse_combined_combination_list:
						hyperedge_coordinates_candidate = [nominated_hedge_correspondent[nominated_hyperedge_vertex_i]]
						hyperedge_coordinates_candidate.extend(coordinate_group)
						hyperedge_coordinates.append(hyperedge_coordinates_candidate)
						if debug_mode:
							print("coordinate: ",nominated_hedge_correspondent[nominated_hyperedge_vertex_i],coordinate_group)
				if debug_mode:
					print("hyperedge_coordinates",hyperedge_coordinates)
				informal_hyperedges.append(hyperedge_coordinates) 
		#print("		fc - O")
		#print("		",orth_axis_combinations)
	support_map = []
	for axis, v1 in universal_map.items():
		support_map.append([])
		for rv, v2 in universal_map[axis].items():
			for state, v3 in universal_map[axis][rv].items():
				support_map[-1].append([generate_string(),[str(axis),str(rv),str(state)]])
	#print("support_map",support_map)
	axis_products = []
	traverse_level = len(support_map)
	def axis_product_traversal(arg_hand_list, arg_traverse_level):
		for i in range(0, len(support_map[len(support_map)-arg_traverse_level])):
			list_to_put = [[len(support_map)-arg_traverse_level, i]]
			list_to_put.extend(arg_hand_list)
			if (arg_traverse_level > 1):
				axis_product_traversal(list_to_put, arg_traverse_level-1)
			else:
				axis_products.append(list_to_put)
	axis_product_traversal([],traverse_level)
	#print("axis_products",axis_products)
	direct_map = {}
	for h in range(len(axis_products)):
		prefix = ""
		coordinates = []
		for i in range(0, len(axis_products[h])):
			prefix += support_map[axis_products[h][i][0]][axis_products[h][i][1]][0]
			coordinates.append(support_map[axis_products[h][i][0]][axis_products[h][i][1]][1])
		direct_map.update({prefix:coordinates})
	if debug_mode:
		print("")
	formal_hyperedges = []
	for ih in informal_hyperedges:
		if debug_mode:
			print(ih)
		formal_hyperedge = []
		for co in ih:
			formal_hyperedge.append(prefix_from_coordinate(direct_map, co))
		formal_hyperedges.append(formal_hyperedge)
	if debug_mode:
		print("formal_hyperedges",formal_hyperedges)
	
	refactored_aliased_hyperedges = []
	for ih in formal_hyperedges:
		if not set(ih) in refactored_aliased_hyperedges:
			refactored_aliased_hyperedges.append(set(ih))
	if debug_mode:
		print("refactored_aliased_hyperedges",refactored_aliased_hyperedges)
	refactored_formal_hyperedges = []
	for rh in refactored_aliased_hyperedges:
		hyperedge = []
		for co in rh:
			hyperedge.append(direct_map[co])
		refactored_formal_hyperedges.append(hyperedge)

	def direct_to_val(coordinate_alias):
		values = direct_map[coordinate_alias]
		identity_ref_object = []
		for value in values:
			constructed_info = universal_map[value[0]][value[1]][value[2]]
			identity_ref_object.append([constructed_info["rv"], constructed_info["state"]])
		return identity_ref_object

	validated_hyperedges = []
	for hyperedge in refactored_aliased_hyperedges:
		hyperedge_to_construct = []
		for coordinate_alias in hyperedge:
			hyperedge_to_construct.append(direct_to_val(coordinate_alias))
		validated_hyperedges.append(hyperedge_to_construct)

	return {
		"direct_map": direct_map,
		"universal_map": universal_map,
		"hyperedges" : validated_hyperedges,
		"refactored_formal_hyperedges" : refactored_formal_hyperedges,
		"refactored_aliased_hyperedges" : refactored_aliased_hyperedges
	}
	

def traverse(**kwargs):
	arg_level = kwargs.get('level', 0)
	arg_starting_level = kwargs.get('starting_level', arg_level)
	arg_major_list = kwargs.get('list', [])
	indicator = kwargs.get('indicator', [])
	extendor_list = kwargs.get('extendor_list', [])
	for i in range(0, indicator[arg_starting_level-arg_level]):
		l = [i]
		l.extend(extendor_list)
		if (arg_level > 1):
			traverse(level=arg_level-1, extendor_list=l, list=arg_major_list, indicator=indicator, starting_level=arg_starting_level)
		else:
			arg_major_list.append(l)

def axis_ref_to_index(universal_map, ref):
	index = 0
	for rv_step, v in universal_map[ref[0]].items():
		# Add the number of states to the index
		if int(rv_step) < int(ref[1]):
			index += len(v)
		elif (int(rv_step) == int(ref[1])):
			index += int(ref[2])
	return index

def indices_to_vector_index(rvs_and_states_in_axes, indices):
	vector_index = indices[0]
	for i in range(1, len(indices)):
		vector_index += rvs_and_states_in_axes[i-1]*indices[i]
	return vector_index

def convert_reading_to_flat(FR, reading):
	i = 0
	for ax in FR["universal_map"]:
		i = 0
		for rv in FR["universal_map"][ax]:
			for st in FR["universal_map"][ax][rv]:
				if ([ax,rv,st] == reading):
					return i+1
				i += 1

def coordinate(limits_on_axes, coords_on_axes):
	# The last axis does not typically have a limit imposed on it
	current_axis = len(coords_on_axes)
	total = 0
	while (current_axis > 0):
		subtractant_from_last_axis = [coords_on_axes[current_axis-1] - 1] # Subtract 1 from last axis
		subtractant_from_last_axis.extend(limits_on_axes[:current_axis-1])
		total += np.product(subtractant_from_last_axis) # Multiply by all
		current_axis -= 1 # Prepare for next axis inwards
	return (total+1) # Return the total

class non_orthogonality(object):
	def __init__(self, **kwargs):
		FR = kwargs.get('FR', None)
		num_of_axes = len(FR["universal_map"])
		#print("\n\n\n")
		lens_of_axes = []
		for axis, axis_v in FR["universal_map"].items():
			#print("axis",axis, axis_v)
			lens_of_axes.append(len(flat([[v for v in axis_v[k]] for k in axis_v])))
		lens_of_axes = lens_of_axes[::-1] # Reverse the axis lengths
		self.lens_of_axes = lens_of_axes
		NO_graph_refactored = {}
		back_map = {}
		for member in FR["direct_map"]:
			#print(member)
			orthogs = []
			for hyperedge in FR["refactored_aliased_hyperedges"]:
				if (member in hyperedge):
					#print("\t",hyperedge)
					orthogs.extend(hyperedge.split(":"))
			NO_graph_refactored.update({member:list(set(FR["direct_map"].keys()).difference(set(orthogs)))})

			coords_flat = []
			vals = []
			for coord in FR['direct_map'][member]:
				coords_flat.append(convert_reading_to_flat(FR, coord))
				vals.append(FR["universal_map"][coord[0]][coord[1]][coord[2]])
			back_map.update({coordinate(lens_of_axes, coords_flat):vals})
		#print("back_map",back_map)



		
		#print("lens_of_axes",lens_of_axes)
		#print("NO_graph_refactored",NO_graph_refactored)
		all_relations = {}
		for from_val, to_vals in NO_graph_refactored.items():
			collected_relation = []
			#print(FR['direct_map'][from_val], to_vals)
			from_coords_flat = []
			for coord in FR['direct_map'][from_val]:
				from_coords_flat.append(convert_reading_to_flat(FR, coord))
			#print("from_coords_flat",from_coords_flat)
			#print("from_coord", coordinate(lens_of_axes, from_coords_flat))
			collected_relation.append(coordinate(lens_of_axes, from_coords_flat))
			to_coords_coords = []
			for refac_hyperedge in to_vals:
				to_coords_flat = []
				for coord in FR['direct_map'][refac_hyperedge]:
					to_coords_flat.append(convert_reading_to_flat(FR, coord))
				to_coords_coords.append(coordinate(lens_of_axes, to_coords_flat)-1)
			#print("to_coords_coords",to_coords_coords)
			collected_relation.append(to_coords_coords)
			all_relations.update({collected_relation[0]-1:collected_relation[1]})
		#print("all_relations",all_relations)

		


		G = nx.networkx_mod.networkx_mod.Graph()

		num_of_vertices = len(back_map.values())
		for i in range(num_of_vertices):
			G.add_node(i)
		for i in range(num_of_vertices):
			for co in all_relations[i]:
				G.add_edge(i, co)

		self.num_of_axes = len(lens_of_axes)
		self.back_map = back_map
		self.num_of_vertices = num_of_vertices
		self.G = G

		self.G_a = []
		for h in range(0, 10):
			G_ = nx.networkx_mod.networkx_mod.Graph()
			for i in range(num_of_vertices):
				G_.add_node(i)
				self.G_a.append(G_)
				#G2.add_edge(i, no_graph[i][j])

			#print()
		'''
		self.cliques = None
		FR = kwargs.get('FR', None)
		num_of_axes = len(FR["universal_map"])
		rvs_and_states_in_axes = []
		rvs_and_states_refs = []
		for axis in FR["universal_map"]:
			num_of_rvs_and_states = 0
			rvs_and_states_refs.append([])
			for rv in FR["universal_map"][axis]:
				#rvs_and_states_refs[-1].append(FR["universal_map"][axis][rv].items())
				for k,v in FR["universal_map"][axis][rv].items():
					num_of_rvs_and_states += 1
					rvs_and_states_refs[-1].append(FR["universal_map"][axis][rv][k])
			rvs_and_states_in_axes.append(num_of_rvs_and_states) 
		print("rvs_and_states_in_axes",rvs_and_states_in_axes)
		#rvs_and_states_in_axes = rvs_and_states_in_axes[::-1]
		combinations_of_axes = []
		traverse(level=num_of_axes, list=combinations_of_axes, indicator=rvs_and_states_in_axes)
		for i in range(len(combinations_of_axes)):
			combinations_of_axes[i] = combinations_of_axes[i][::-1]
		back_map = []
		#print(combinations_of_axes)
		print("rvs_and_states_refs",rvs_and_states_refs)
		print("combinations_of_axes", combinations_of_axes)
		for i in range(len(combinations_of_axes)):
			combination = combinations_of_axes[i]
			back_map_part = []
			for given_axis in range(len(combination)):
				index_in_axis = combination[given_axis]
				back_map_part.append(rvs_and_states_refs[given_axis][index_in_axis])
			back_map.append(back_map_part)
		
		hyperedge_vector_indices = []
		for hyperedge in FR["refactored_formal_hyperedges"]:
			hyperedge_vector_indices.append([])
			for coordinates_list in hyperedge:
				indices = []
				for index in coordinates_list:
					indices.append(axis_ref_to_index(FR["universal_map"], index))
				hyperedge_vector_indices[-1].append(indices_to_vector_index(rvs_and_states_in_axes, indices))
		print("hyperedge_vector_indices",hyperedge_vector_indices)
		o_graph = []
		no_graph = []

		num_of_vertices = len(back_map)
		for i in range(num_of_vertices):
			o_graph.append([])
			no_graph.append([])
			for hvi_i in hyperedge_vector_indices:
				if i in hvi_i:
					o_graph[-1].extend(hvi_i)
			o_graph[-1] = set(o_graph[-1])
			for j in range(num_of_vertices):
				if (j != i) and (j not in o_graph[-1]):
					no_graph[-1].append(j)

		G = nx.networkx_mod.networkx_mod.Graph()
		for i in range(num_of_vertices):
			G.add_node(i)
		for i in range(num_of_vertices):
			for j in range(len(no_graph[i])):
				G.add_edge(i, no_graph[i][j])
				#G2.add_edge(i, no_graph[i][j])
		print("o_graph",o_graph)
		print("no_graph",no_graph)
		self.num_of_axes = num_of_axes
		self.back_map = back_map
		self.num_of_vertices = num_of_vertices
		self.rvs_and_states_in_axes = rvs_and_states_in_axes
		self.hyperedge_vector_indices = hyperedge_vector_indices
		self.G = G

		self.G_a = []
		for h in range(0, 10):
			G_ = nx.networkx_mod.networkx_mod.Graph()
			for i in range(num_of_vertices):
				G_.add_node(i)
				self.G_a.append(G_)
		'''
	def graph(self, dpi=80, scale=None):
		fig = plt.figure(dpi=dpi)
		figsize_scale = None
		if (scale == None):
			figsize_scale = (self.num_of_vertices/32.0)
		else:
			figsize_scale = scale
		fig.set_size_inches(self.num_of_vertices*figsize_scale,self.num_of_vertices*figsize_scale)
		plt.margins(0.1)
		dt_type = "dark"
		lbl = {}
		
		for k,v in self.back_map.items():
			label_string = ""
			for sup_member in v:
				for i in range(len(sup_member['rv'])):
					label_string += r"$\bf{%s}$ = " % (sup_member['rv'][i].label)
					label_string += r"$\it{%s}$" % (sup_member['state'][i].label)
					label_string += "\n"

			lbl.update({k-1:label_string[:-1]})
			'''
			label_string = ""
			label_ref = self.back_map[i]
			for j in range(len(label_ref)):
				label_string += r"$\bf{%s}$ = " % (label_ref[j]['rv'][0].label)
				label_string += r"$\it{%s}$" % (label_ref[j]['state'][0].label)
				if (j != (len(label_ref) - 1)):
					label_string += "\n"
			lbl.update({i:label_string})
			'''
		
		'''
		for h in range(0, 10):
			nx.draw(self.G_a[h], pos=nx.circular_layout(self.G_a[h]), with_labels=True, font_weight="normal", 
				node_size=250*self.num_of_vertices*figsize_scale, 
				node_color=(0,0,0,0),
				node_alpha=1.0,
				node_edge_alpha_formal=1.0,
				font_color=(1,1,1,1),
				font_size=display_theme()[dt_type]["graphs"]["font_size"]*self.num_of_vertices*figsize_scale,
				edge_color=display_theme()[dt_type]["color_bank"][8],
				font_x_offset=0, font_y_offset=-0.002*self.num_of_vertices*figsize_scale,
				edge_alpha=(h/10.0), 
				edge_width=display_theme()[dt_type]["graphs"]["edge_width"]*h, 
				edgecolors=matplotlib.colors.to_hex(display_theme()[dt_type]["line_color_alpha"]), 
				label_pos=0.5, verticalalignment='center', labels=lbl)
		
		'''
		'''
		nx.draw(self.G2, pos=nx.circular_layout(self.G2), with_labels=True, font_weight="normal", 
			node_size=250*self.num_of_vertices*figsize_scale, 
			node_color=matplotlib.colors.to_hex(display_theme()[dt_type]["color_bank_alpha_25"][5]),
			node_alpha=display_theme()[dt_type]["graphs"]["node_alpha"],
			node_edge_alpha_formal=display_theme()[dt_type]["graphs"]["node_border_alpha"],
			font_color=display_theme()[dt_type]["font_color"],
			font_size=display_theme()[dt_type]["graphs"]["font_size"]*self.num_of_vertices*figsize_scale,
			edge_color=matplotlib.colors.to_hex(display_theme()[dt_type]["line_color_alpha"]),
			font_x_offset=0, font_y_offset=-0.002*self.num_of_vertices*figsize_scale,
			edge_alpha=display_theme()[dt_type]["graphs"]["edge_alpha"], 
			edge_width=display_theme()[dt_type]["graphs"]["edge_width"], 
			edgecolors=matplotlib.colors.to_hex(display_theme()[dt_type]["line_color_alpha"]), 
			label_pos=0, verticalalignment='center', labels=lbl)
		
		'''
		
		#print(os.path.abspath(inspect.getfile(nx.draw)))
		nx.networkx_mod.networkx_mod.draw(self.G, pos=nx.networkx_mod.networkx_mod.circular_layout(self.G), with_labels=True, font_weight="normal", 
			node_size=250*self.num_of_vertices*figsize_scale, 
			node_color=matplotlib.colors.to_hex(display_theme()[dt_type]["color_bank_alpha_25"][5]),
			node_alpha=display_theme()[dt_type]["graphs"]["node_alpha"],
			node_edge_alpha_formal=display_theme()[dt_type]["graphs"]["node_border_alpha"],
			font_color=display_theme()[dt_type]["font_color"],
			font_size=display_theme()[dt_type]["graphs"]["font_size"]*self.num_of_vertices*figsize_scale,
			edge_color=matplotlib.colors.to_hex(display_theme()[dt_type]["line_color_alpha"]),
			font_x_offset=0,
			font_y_offset=-0.002*self.num_of_vertices*figsize_scale,
			edge_alpha=display_theme()[dt_type]["graphs"]["edge_alpha"], 
			edge_width=display_theme()[dt_type]["graphs"]["edge_width"], 
			edgecolors=matplotlib.colors.to_hex(display_theme()[dt_type]["line_color_alpha"]), 
			label_pos=0.5, verticalalignment='center', labels=lbl)
		
		fig.set_facecolor((0, 0, 0, 0))

	def is_classical(self, arg_probability_space):
		# Procure the flat outcomes of the probabilistic model
		flat_outcomes = [0.0 for i in range(self.num_of_vertices)]
		for i in range(self.num_of_vertices):
			# Test coordinate

			for combination in arg_probability_space.events_register:
				combination_alias = set()
				for rvs in combination['combination']:
					combination_alias.add(rvs[0].id+'_'+rvs[1].id)
				combination_alias_string = ":".join(sorted(combination_alias))
				for k,v in self.back_map.items():
					back_map_alias = set()
					for j in range(len(v)):
						for l in range(len(v[j]['rv'])):
							back_map_alias.add(v[j]['rv'][l].id+'_'+v[j]['state'][l].id)
					back_map_alias_string = ":".join(sorted(back_map_alias))
					#print(back_map_alias_string, combination_alias_string)
					if (back_map_alias_string == combination_alias_string):
						flat_outcomes[k-1] = combination['p']
		self.flat_outcomes = flat_outcomes
		#print("flat_outcomes", flat_outcomes)
		#sys.exit()
		# Determine the cliques
		#print(dir(nx))
		#print(dir(nx.networkx_mod.networkx_mod))
		#print(dir(nx.networkx_mod.networkx_mod.networkx_mod))
		cliques = list(nx.networkx_mod.networkx_mod.networkx_mod.algorithms.clique.enumerate_all_cliques(self.G))
		self.cliques = cliques
		# Create a probability to constraint map
		p_to_c_map = {}
		for i in range(self.num_of_vertices):
			p_index = "%s"%(i)
			p_to_c_map[p_index] = []
			for c in range(len(cliques)):
				if (i in cliques[c]):
					p_to_c_map[p_index].append(c)
		QLP  = pic.Problem()
		CL = [ QLP.add_variable('CL['+str(i)+']',1) for i in range(len(cliques)) ]
		PR = [ QLP.add_variable('PR['+str(i)+']',1) for i in range(self.num_of_vertices) ]

		for i in range(self.num_of_vertices):
			QLP.add_constraint( 1|PR[i] == flat_outcomes[i])

		k = 0
		for i in range(len(cliques)):
			QLP.add_constraint( 1|CL[i] >= 0)
			k += CL[i]

		QLP.add_constraint( 1|k == 1)

		cl_list = []
		for k, v in p_to_c_map.items():
			j = 0
			for l in v:
				j += CL[l]
			QLP.add_constraint( 1|PR[int(k)] <= j)
			cl_list.append(j)

		QLP.solve(verbose=0)
		#print(QLP)
		#print(CL[147].is_valued())#str(cl_list[0])))
		#print(eval('round(CL[147],3)'))
		classical = True
		for i in range(0, len(cliques)):
			if (not CL[i].is_valued()):
				classical = False
		return classical
	def minimal_k(self):
		cliques = list(nx.networkx_mod.networkx_mod.networkx_mod.algorithms.clique.enumerate_all_cliques(self.G))
		self.cliques = cliques
		# Create a probability to constraint map
		p_to_c_map = {}
		for i in range(self.num_of_vertices):
			p_index = "%s"%(i)
			p_to_c_map[p_index] = []
			for c in range(len(cliques)):
				if (i in cliques[c]):
					p_to_c_map[p_index].append(c)
		QLP  = pic.Problem()
		CL = [ QLP.add_variable('CL['+str(i)+']',1) for i in range(len(cliques)) ]
		PR = [ QLP.add_variable('PR['+str(i)+']',1) for i in range(self.num_of_vertices) ]
		for i in range(self.num_of_vertices):
			QLP.add_constraint( 1|PR[i] == self.flat_outcomes[i])
		K = QLP.add_variable('K',1)
		QLP.add_constraint( 1| K >= 0)
		Allc = 0
		for i in range(len(cliques)):
			QLP.add_constraint( 1|CL[i] >= 0)
			Allc += CL[i]
		QLP.add_constraint( 1| K == Allc)
		for k, v in p_to_c_map.items():
			j = 0
			for l in v:
				j += CL[l]
			QLP.add_constraint( 1|PR[int(k)] <= j)
		#
		QLP.set_objective('min',K)
		QLP.solve(verbose=0)
		return K
		#
	def determine_contextuality(self, **kwargs):
		FR = kwargs.get('FR', None)
		clique_list_4 = []
		for i in range(len(self.cliques)):
			if len(self.cliques[i]) == 4: # TODO generalise this
				clique_list_4.append(list(np.array(self.cliques[i])+1))
		combinations_ = list(itertools.product([0,1], repeat=len(clique_list_4)))
		isolated_backmap_refs_all = []
		for c in range(0, len(combinations_)):
			try:
				clique_list_comb = []
				for i in range(0,len(combinations_[c])):
					if combinations_[c][i] > 0:
						clique_list_comb.append(clique_list_4[i])
				if (False):
					if (c % 100 == 0):
						print(c)
				cond_a = list(set([item for sublist in clique_list_comb for item in sublist])) == [x for x in range(1,17)] # does it touch every vertex
				l = [item for sublist in clique_list_comb for item in sublist]
				d = {}
				for i in l: d[i] = i in d
				cond_e = len([k for k in d if not d[k]]) > 0 # Does it have 8 distinct vertices
				contexts = [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]] # would vary with type of system
				clique_list_comb_unique = []
				for elem in clique_list_comb:
					for vertex in elem:
						duplic = False
						for elem_b in clique_list_comb:
							if (vertex in elem_b and (elem_b is not elem)):
								duplic = True
						if not duplic:
							clique_list_comb_unique.append(vertex)
				clique_list_comb_counts = []
				for i in clique_list_comb:
					total = 0
					for num in clique_list_comb_unique:
						if (num in i):
							total += 1
					clique_list_comb_counts.append(total)
				cond_g = (len(set(clique_list_comb_counts)) == 1)
				num_uniques = set(clique_list_comb_counts)
				cond_i = (list(num_uniques)[0] != np.product(FR['num_rvs_all_axes']))
				if ((cond_a) and cond_e and cond_g and cond_i):
					QLP  = pic.Problem()
					p = [ QLP.add_variable('p['+str(i)+']',1) for i in range(16) ]
					CL = [ QLP.add_variable('CL['+str(i)+']',1) for i in range(len(clique_list_comb)) ]
					for i in range(len(p)):
						QLP.add_constraint( 1| p[i] >= 0)
					for i in range(len(CL)):
						summation = 0
						for j in clique_list_comb[i]:
							summation += p[j-1]
						QLP.add_constraint( 1| CL[i] == summation)
					QLP.solve(verbose=0)
					isolated_backmap_refs = []
					for k,v in dict(collections.Counter(flat(clique_list_comb))).items():
						if v == 1:
							isolated_backmap_refs.append(k)
					isolated_backmap_refs_all.append(isolated_backmap_refs)
					if (False):
						for cl in CL:
							print(round(cl,5))
						print("\n")
						for pp in p:
							print(round(pp,5))
						print(collections.Counter(flat(clique_list_comb)))
			except:
				pass
		remove_list = []
		correlation_pairs = []
		for member in isolated_backmap_refs_all:
			for member_b in isolated_backmap_refs_all:
				dummy_member = []
				dummy_member.extend(member)
				dummy_member.extend(member_b)
				if (len(set(dummy_member)) == self.num_of_vertices):
					if (not member in remove_list) and (not member_b in remove_list):
						remove_list.append(member)
						remove_list.append(member_b)
						correlation_pairs.append([member,member_b])
		#print("correlation_pairs",correlation_pairs)
		self.correlation_pairs = correlation_pairs
		correlations = []
		for member in correlation_pairs:
			correlation = 0
			correlation_a = 0
			correlation_b = 0
			member_refac_a = (np.array(member[0])-1)
			member_refac_b = (np.array(member[1])-1)
			for i in range(len(self.flat_outcomes)):
				if i in member_refac_a:
					correlation_a += self.flat_outcomes[i]
				if i in member_refac_b:
					correlation_b += self.flat_outcomes[i]
			if ([0,2] in member_refac_a):
				correlations.append(abs(correlation_b-correlation_a))
			else:
				correlations.append(abs(correlation_a-correlation_b))
		#print("correlations",correlations)

		self.signalling_box_values = correlations
		axes = []
		for axis in FR['universal_map']:
			rvs = []
			for k,v in FR['universal_map'][axis].items():
				#print(v)
				rv = []
				for kk,vv in v.items():
					rv_st_str = ""
					for i in range(len((vv['rv']))):
						rv_st_str += vv['rv'][i].id+"_"+vv['state'][i].id+"."
					rv.append(rv_st_str[:-1])
				rvs.append(rv)
			axes.append(rvs)
		contexts = []
		for element in itertools.product(*axes):
			contexts.append([elem for elem in itertools.product(*element)])
		contexts_hyperedges = []
		for c in contexts:
			this_context_set = set()
			for cood in c:
				this_context_set.add(".".join(sorted(cood)))
			contexts_hyperedges.append(":".join(sorted(this_context_set)))
		paired_hyperedges = set()
		for hyperedge_a in FR['refactored_aliased_hyperedges']:
			if (not hyperedge_a in contexts_hyperedges):
				for hyperedge_b in FR['refactored_aliased_hyperedges']:
					if (not hyperedge_b in contexts_hyperedges):
						together_coods = hyperedge_a.split(":")
						together_coods.extend(hyperedge_b.split(":"))
						together_coods = sorted(set(together_coods))
						if len(set(hyperedge_a.split(":")).intersection(set(hyperedge_b.split(":")))) == 0:
							for c in contexts_hyperedges:
								if set(c.split(":")).issubset(together_coods):
									paired_hyperedges.add("!".join(sorted(set([hyperedge_a,hyperedge_b]))))
		self.paired_hyperedges = list(paired_hyperedges)
		#print("paired_hyperedges",self.paired_hyperedges)
		all_hyperedge_coods = []
		for pair in self.paired_hyperedges:
			this_top_hyperedge_coods = []
			for i in range(2):
				this_hyperedge = pair.split("!")[i]
				this_hyperedge_coods = []
				for reading_top in (this_hyperedge.split(":")):
					val = FR['direct_map'][reading_top]
					coords = []
					for reading in val:
						coords.append(convert_reading_to_flat(FR, reading))
					this_coordinate = coordinate(self.lens_of_axes,coords)
					this_hyperedge_coods.append(this_coordinate)
				this_top_hyperedge_coods.append(this_hyperedge_coods)
			all_hyperedge_coods.append(this_top_hyperedge_coods)
		self.all_hyperedge_coods = all_hyperedge_coods
		#
		totals_marginals = []
		for pair in self.all_hyperedge_coods:
			total_amt_a = 0
			total_amt_b = 0
			for i in range(len(self.flat_outcomes)):
				if i+1 in pair[0]:
					total_amt_a += self.flat_outcomes[i]
				if i+1 in pair[1]:
					total_amt_b += self.flat_outcomes[i]
			totals_marginals.append(abs(total_amt_a - total_amt_b)/2) # Modified here for Aerts
		self.totals_marginals = totals_marginals
		self.marginal_violation = (sum(self.totals_marginals)/2.0)
		self.maximum_allowed_violation = max([(max(self.signalling_box_values)/2.0)-1, self.marginal_violation])












