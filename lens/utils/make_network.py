"""
Make formatted files for Gephi Network Visualization

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division

import os
import csv


def get_loose_nodes(stoichiometry):
	# get all needed exchanges
	nodes, edges = make_network(stoichiometry)

	sources = set()
	targets = set()
	for edge in edges:
		source = edge[0]
		target = edge[1]
		sources.add(source)
		targets.add(target)

	loose_nodes = sources.symmetric_difference(targets)

	return loose_nodes

def make_network(stoichiometry, node_type={}):
	'''Build the network'''

	nodes = {}
	edges = []

	for reaction_id, stoich in stoichiometry.iteritems():

		# add reaction to node list
		rxn_type = 'reaction' # TODO -- pass in type
		nodes[reaction_id] = {'label': reaction_id, 'type': rxn_type}
		if node_type.get(reaction_id):
			n_type = node_type.get(reaction_id)
			nodes[reaction_id] = {'label': reaction_id, 'type': n_type}

		# add molecules to node list, and connections to edge list
		for molecule_id, coeff in stoich.iteritems():

			mol_type = 'molecule' # TODO -- pass in type
			nodes[molecule_id] = {'label': molecule_id, 'type': mol_type}
			if node_type.get(molecule_id):
				n_type = node_type.get(molecule_id)
				nodes[molecule_id] = {'label': molecule_id, 'type': n_type}

			## add edge between reaction and molecule
			# a reactant
			if coeff < 0:
				edge = [molecule_id, reaction_id]
			# a product
			elif coeff > 0:
				edge = [reaction_id, molecule_id]
			else:
				print(reaction_id + ', ' + molecule_id + ': coeff = 0')

			edges.append(edge)

	return nodes, edges


def make_gephi_network(stoichiometry, plotOutDir, node_type={}):

	nodes, edges = make_network(stoichiometry, node_type)

	out_dir = os.path.join(plotOutDir, 'network')

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	nodes_out = os.path.join(out_dir, 'nodes.csv')

	edges_out = os.path.join(out_dir, 'edges.csv')

	## Save network to csv

	# nodes list
	with open(nodes_out, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')

		# write header
		writer.writerow(['Id', 'Label', 'Type'])

		for node, specs in nodes.iteritems():
			label = specs['label']
			type = specs['type']

			row = [node, label, type]
			writer.writerow(row)

	# edges list
	with open(edges_out, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')

		# write header
		writer.writerow(['Source', 'Target'])

		for edge in edges:
			source = edge[0]
			target = edge[1]

			row = [source, target]
			writer.writerow(row)
