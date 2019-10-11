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

def make_network(stoichiometry, info={}):
	'''
	Build the network
	TODO -- add edge weight to info
	'''
	node_types = info.get('node_types', {})
	node_sizes = info.get('node_sizes', {})

	nodes = {}
	edges = []

	for reaction_id, stoich in stoichiometry.iteritems():

		# add reaction to node list
		n_type = node_types.get(reaction_id, 'reaction')
		n_size = node_sizes.get(reaction_id, 1)
		nodes[reaction_id] = {
			'label': reaction_id,
			'type': n_type,
			'size': n_size}

		# add molecules to node list, and connections to edge list
		for molecule_id, coeff in stoich.iteritems():
			n_type = node_types.get(molecule_id, 'molecule')
			n_size = node_sizes.get(molecule_id, 1)
			nodes[molecule_id] = {
				'label': molecule_id,
				'type': n_type,
				'size': n_size}

			## add edge between reaction and molecule
			# a reactant
			if coeff < 0:
				edge = [molecule_id, reaction_id]
			# a product
			elif coeff > 0:
				edge = [reaction_id, molecule_id]
			else:
				print(reaction_id + ', ' + molecule_id + ': coeff = 0')
				break
			edges.append(edge)

	return nodes, edges

def save_network(stoichiometry, plotOutDir='out/network', info={}):
	'''
	Makes a gephi network with weighted edges.
	stoichiometry and weights need to have the same reaction ids
	info can contain two types of data in dictionaries, edge_weights and node_types
	info = {
		'node_sizes': node_sizes (dict),
		'node_types': node_types (dict)
	}
	weights = {
	'''
	nodes, edges = make_network(stoichiometry, info)

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
		writer.writerow(['Id', 'Label', 'Type', 'Size'])

		for node, specs in nodes.iteritems():
			label = specs['label']
			type = specs['type']
			size = specs['size']

			row = [node, label, type, size]
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
