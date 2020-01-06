# THIS FILE IS EXPERIMENTAL AND INCOMPLETE
# TODO(Ryan): finish the bigraph class

def extract_sites(place):
	sites = []
	for node, children in place.items():
		if isinstance(children, dict):
			sites += extract_sites(children)
		else:
			sites.append(node)
	return sites

def extract_ports(link, label):
	ports = []
	for edge in link.values():
		ports += edge.get(label, [])
	return ports

class Bigraph(object):
	'''
	Bigraph: superimposed trees and hypergraphs with composition and dynamic transition rules.

	A bigraph consists of a set of nodes and an independent nested place graph and a link hypergraph
	involving those nodes.

	The place graph (essentially a tree) will be of the form:

	  {'a': {'b': {0: 0}}, 'c': {'d': {}, 'e': {}}}, ....}

	The hypergraph will be a dict of edge names to connected nodes of the form:

	  {'x': {'nodes': ['a', 'b', 'c'], 'outer': ['A'], 'inner': []},
	   'y': {'nodes': ['b', 'c'], 'outer': ['B'], 'inner': ['C', 'D']},
	   ....}

	For the place graph, roots are represented as top level keys ('a' and 'c' above) and
	sites are represented as nodes whose value is a site label (instead of a further nested dict).

	For link graphs, the 'nodes' key contains all nodes within the bigraph that are connected
	by this edge, and outer and inner names are contained under the keys 'outer' and 'inner'
	respectively.
	'''

	def __init__(self, signature, place, link):
		self.signature = signature
		self.place = place
		self.link = link
		self.assign_ports()

	def derive_ports(self):
		return {
			'roots': self.place.keys(),
			'sites': extract_sites(self.place),
			'outer': extract_ports(self.link, 'outer'),
			'inner': extract_ports(self.link, 'inner')}

	def assign_ports(self):
		ports = self.derive_ports()
		self.roots = ports['roots']
		self.sites = ports['sites']
		self.outer = ports['outer']
		self.inner = ports['inner']

	def compose(self, other):
		return self

# Namespaces

# N - node
# C - control
# P - port
# E - (hyper)edge
# M - name
# I - inner name
# O - outer name
# S - site

signature = {
	'C:agent': 2,
	'C:building': 1,
	'C:computer': 2,
	'C:room': 0}

example_place = {
	'N:agent:a': {
		'control': 'C:agent',
		'ports': ['E:call:callA', 'E:login:1111'],
		'children': [],
		'sites': []},
	'N:agent:b': {
		'control': 'C:agent',
		'ports': ['E:call:callA', None],
		'children': [],
		'sites': []},
	'N:agent:c': {
		'control': 'C:agent',
		'ports': ['E:call:callA', 'E:login:2222'],
		'children': [],
		'sites': []},
	'N:agent:d': {
		'control': 'C:agent',
		'ports': ['E:call:callA', 'E:login:3333'],
		'children': [],
		'sites': []},
	'N:agent:e': {
		'control': 'C:agent',
		'ports': ['E:call:callA', 'E:login:4444'],
		'children': [],
		'sites': []},
	'N:building:1': {
		'control': 'C:building',
		'ports': ['E:network:network1'],
		'children': ['N:room:X1', 'N:agent:b', 'N:room:X2'],
		'sites': ['S:building:1']},
	'N:building:2': {
		'control': 'C:building',
		'ports': ['E:network:network2'],
		'children': ['N:room:X3', 'N:room:X4'],
		'sites': ['S:building:2']},
	'N:computer:X1': {
		'control': 'C:computer',
		'ports': ['E:login:1111', 'E:network:network1'],
		'children': [],
		'sites': []},
	'N:computer:X2': {
		'control': 'C:computer',
		'ports': ['E:login:2222', 'E:network:network1'],
		'children': [],
		'sites': []},
	'N:computer:X3': {
		'control': 'C:computer',
		'ports': ['E:login:3333', 'E:network:network2'],
		'children': [],
		'sites': []},
	'N:computer:X4': {
		'control': 'C:computer',
		'ports': ['E:login:4444', 'E:network:network2'],
		'children': [],
		'sites': []},
	'N:room:X1': {
		'control': 'C:room',
		'ports': [],
		'children': ['N:agent:a', 'N:computer:X1'],
		'sites': []},
	'N:room:X2': {
		'control': 'C:room',
		'ports': [],
		'children': ['N:agent:c', 'N:computer:X2'],
		'sites': []},
	'N:room:X3': {
		'control': 'C:room',
		'ports': [],
		'children': ['N:agent:d', 'N:computer:X3'],
		'sites': []},
	'N:room:X4': {
		'control': 'C:room',
		'ports': [],
		'children': ['N:agent:e', 'N:computer:X4'],
		'sites': []}}

example_place = {
	'a': {
		'b': {
			0: 0}},
	'c': {
		'd': {},
		'e': {
			'f': {1: 1}}}}

example_link = {
	'x': {
		'nodes': ['a', 'b', 'c'],
		'outer': ['A'],
		'inner': []},
	'y': {
		'nodes': ['b', 'd', 'f'],
		'outer': ['B'],
		'inner': ['C']},
	'z': {
		'nodes': ['e', 'g'],
		'outer': [],
		'inner': ['D', 'E']}}

example_bigraph = Bigraph([], example_place, example_link)
