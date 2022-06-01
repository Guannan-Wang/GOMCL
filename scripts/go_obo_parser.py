#!/usr/bin/env python

__author__    = "Guannan Wang"
__email__     = "wangguannan2014@gmail.com"
__version__   = "1"

import sys, re, argparse, itertools
from collections import defaultdict
try:
	import networkx as nx
except ImportError as ErrorMessage:
	print(str(ErrorMessage))
	sys.exit()
	
def obo_parser(OBOInput):
	fin_obo = open(OBOInput,"rU")
	lines_obo = fin_obo.readlines()
	obo_go_id = defaultdict(str)
	obo_go_alt_id_dict = defaultdict(list)
	obo_go_namespace_dict = defaultdict(str)
	obo_go_name_dict = defaultdict(str)
	obo_go_synonym_dict = defaultdict(list)
	obo_go_is_obsolete_dict = defaultdict(str)
	obo_go_replaced_by_dict = defaultdict(str)
	obo_go_is_a_dict = defaultdict(list)
	obo_go_relationship_dict = defaultdict(list)
	obo_go_def_dict = defaultdict(list)
	obo_comment_dict = defaultdict(list)
	for linenum_obo in range(len(lines_obo)):
		if lines_obo[linenum_obo].strip() == "[Term]":
#			obo_go_numID = int(lines_obo[linenum_obo + 1].strip().split("id: GO:")[1])
			obo_go_id = lines_obo[linenum_obo + 1].strip().split("id: ")[1].strip()
			obo_go_name_dict[obo_go_id]= lines_obo[linenum_obo + 2].strip().split("name: ")[1].strip()
			obo_go_namespace_dict[obo_go_id] = lines_obo[linenum_obo + 3].strip().split("namespace: ")[1].strip()
			obo_go_marker = 4
			while lines_obo[linenum_obo + obo_go_marker].strip() not in ["[Term]","[Typedef]"]:
				if lines_obo[linenum_obo + obo_go_marker].strip().startswith("alt_id: "):
					obo_go_alt_id_dict[obo_go_id].append(lines_obo[linenum_obo + obo_go_marker].strip().split("alt_id: ")[1].strip())
				elif lines_obo[linenum_obo + obo_go_marker].strip().startswith("synonym: "):
					obo_go_synonym_dict[obo_go_id].append(lines_obo[linenum_obo + obo_go_marker].strip().split("synonym: ")[1].strip())
				elif lines_obo[linenum_obo + obo_go_marker].strip().startswith("is_obsolete: "):
					obo_go_is_obsolete_dict[obo_go_id] = str(lines_obo[linenum_obo + obo_go_marker].strip().split("is_obsolete: ")[1].strip())
				elif lines_obo[linenum_obo + obo_go_marker].strip().startswith("replaced_by: "):
					obo_go_replaced_by_dict[obo_go_id] = str(lines_obo[linenum_obo + obo_go_marker].strip().split("replaced_by: ")[1].strip())
				elif lines_obo[linenum_obo + obo_go_marker].strip().startswith("is_a: "):
					obo_go_is_a_dict[obo_go_id].append(re.split("is_a: |!",lines_obo[linenum_obo + obo_go_marker].strip())[1].strip())
				elif lines_obo[linenum_obo + obo_go_marker].strip().startswith("relationship: part_of"):
					obo_go_relationship_dict[obo_go_id].append(re.split("relationship: part_of|!",lines_obo[linenum_obo + obo_go_marker].strip())[1].strip())
				elif lines_obo[linenum_obo + obo_go_marker].strip().startswith("def: "):
					obo_go_def_dict[obo_go_id].append(lines_obo[linenum_obo + obo_go_marker].strip().split("def: ")[1].strip())
				elif lines_obo[linenum_obo + obo_go_marker].strip().startswith("comment: "):
					obo_comment_dict[obo_go_id].append(lines_obo[linenum_obo + obo_go_marker].strip().split("comment: ")[1].strip())
				obo_go_marker += 1
#			print(str(obo_go_id) + "\t" + str(obo_go_marker))
	return obo_go_name_dict, obo_go_namespace_dict, obo_go_alt_id_dict, obo_go_is_obsolete_dict, obo_go_is_a_dict, obo_go_relationship_dict

def construct_go_hierarchy_digraph(OBOInput):
	go_name_dict, go_namespace_dict, alt_id_dict, go_is_obsolete_dict, go_is_a_dict, go_relationship_dict = obo_parser(OBOInput)
	go_hierarchy_digraph = nx.DiGraph()
	for go_id in go_name_dict:
		if go_is_obsolete_dict[go_id] != "true":
			go_hierarchy_digraph.add_node(go_id, Name = go_name_dict[go_id], Namespace = go_namespace_dict[go_id], AltID = alt_id_dict[go_id], IsA = go_is_a_dict[go_id], Relationship = go_relationship_dict[go_id])
			go_hierarchy_digraph.add_edges_from([(parent_go_id,go_id) for parent_go_id in itertools.chain(go_is_a_dict[go_id], go_relationship_dict[go_id])], weight = 1) ## This reads the GO terms with following relationships into the network: "is_a" and "part_of".
#			go_hierarchy_digraph.add_edges_from([(parent_go_id,go_id) for parent_go_id in go_is_a_dict[go_id]], weight = 1)
	return go_hierarchy_digraph


def get_depth_levels(OBOInput):
	"""
	Measure depth of each GO term to the root terms: GO:0008150 (biological_process),GO:0003674 (molecular_function), GO:0005575 (cellular_component)
	This part should be used with caution, because calculating longest_simple_path_length takes a lot of time for large dataset. Suggestion: 'construct_go_hierarchy_digraph', and calculate longest_simple_path_length when needed
	"""
	
	RootTerm_dict = {"biological_process":"GO:0008150", "molecular_function":"GO:0003674", "cellular_component":"GO:0005575"}
	
	## Create weighted GO network and the negative weighted version for searching for longest path between two nodes
	go_hierarchy_network = construct_go_hierarchy_digraph(OBOInput)
	go_hierarchy_network_neg = nx.DiGraph(go_hierarchy_network)
	for u, v in go_hierarchy_network_neg.edges():
		go_hierarchy_network_neg[u][v]['weight'] *= -1
	
	go_depth_dict = dict()
	go_level_dict = dict()
	for node_GOID in list(go_hierarchy_network.nodes(data = False)):  # sorted(list(go_hierarchy_network.nodes(data = False)), key = lambda k:int(k.split("GO:")[1])):	
		node_GO_numID = int(node_GOID.split("GO:")[1])
#		nx.shortest_path(go_hierarchy_network, source = RootTerm_dict[go_hierarchy_network.nodes[node_GOID]["Namespace"]], target = node_GOID)
#		nx.all_simple_paths(go_hierarchy_network, source = RootTerm_dict[go_hierarchy_network.nodes[node_GOID]["Namespace"]], target = node_GOID) ## A simple path is a path with no repeated nodes. This returns a generator.	
#		len(max(nx.all_simple_paths(go_hierarchy_network, source = RootTerm_dict[go_hierarchy_network.nodes[node_GOID]["Namespace"]], target = node_GOID), key = len)) ## The use of "all_simple_paths" significantly reduce running effciency.		
#		nx.bellman_ford_path(go_hierarchy_network_neg, source = RootTerm_dict[go_hierarchy_network.nodes[node_GOID]["Namespace"]], target = node_GOID) ## Returns the shortest path from source to target in a weighted graph G.
		go_level_dict[node_GO_numID] = nx.shortest_path_length(go_hierarchy_network,source = RootTerm_dict[go_hierarchy_network.nodes[node_GOID]["Namespace"]],target = node_GOID)
		go_depth_dict[node_GO_numID] = nx.bellman_ford_path_length(go_hierarchy_network_neg, source = RootTerm_dict[go_hierarchy_network.nodes[node_GOID]["Namespace"]], target = node_GOID)
#		print(str(node_GOID) + "\t" + str(go_level_dict[node_GO_numID]) + "\t" + str(go_depth_dict[node_GO_numID]))
	return go_depth_dict, go_level_dict


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "go_obo_parser.py parses raw obo file to a tabulated table", formatter_class = argparse.RawTextHelpFormatter)
	parser.add_argument('OBOInput', metavar = '-OBO', help = 'obo file should be provided, e.g. go-basic.obo', action = 'store', nargs = None, const = None, default = None, type = None, choices = None)
	args = parser.parse_args()
	
	go_hierarchy_network = construct_go_hierarchy_digraph(args.OBOInput)
	for node_go_id in list(go_hierarchy_network.nodes(data = False)):
		print("Processing " + str(node_go_id))
