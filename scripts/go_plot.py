#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__    = "Guannan Wang"
__email__     = "wangguannan2014@gmail.com"
__version__   = "1"

import argparse, sys, re, os, io, operator 
from decimal import Decimal ## This is for float comparison when accuracy is important!
from collections import OrderedDict 
import numpy as np
import pandas as pd
import math
import random
from go_obo_parser import *
from go_clustering import *
from funs import *

try:
	import networkx as nx
except ImportError as ErrorMessage:
	print(str(ErrorMessage))
	sys.exit()

try:	
	import matplotlib.pyplot as plt
	from matplotlib.ticker import FuncFormatter
except ImportError as ErrorMessage:
	print(str(ErrorMessage))
	sys.exit()	
	
try:
	import seaborn as sns
except ImportError as ErrorMessage:
	print(str(ErrorMessage))
	sys.exit()


#####################################################################################################################################################################

def sim_plot(enGOclstr, SI, Ct, Inf):
	"""
	plot a similarity map for clustered GO
	enGOclstr	formatted GO file
	SI	similarity index between GOs, from ["JC","OC"]
	Ct	clustering threshold for the overlapping ratio between two GOs, any value between 0 and 1
	Inf	inflation value, main handle for cluster granularity, usually chosen somewhere in the range [1.2-5.0]
	"""
	enGOclstr_info_dict, gosim_dict = go_compare(enGOclstr, 2, 0, SI) ## This is for clustered GO file where the second column is the GO ID and the last column is the genes, pay attention to the use of "0" here.
	clstredGOlist = rowtolist(enGOclstr, 2, "\t", "Y")
	uniclstrids = unielement(enGOclstr, 1, "\t", "Y")
	clstrlist = rowtolist(enGOclstr, 1, "\t", "Y")
	
	with open(os.path.splitext(enGOclstr)[0] + ".simat", "w") as fin_simindex:
		fin_simindex.write("" + "\t" + "\t".join(map(str,clstredGOlist)) + "\n")
		for clstred_querygo in clstredGOlist:
			clstred_querygo_sim_list = [gosim_dict[clstred_querygo][clstred_subjectgo] for clstred_subjectgo in clstredGOlist]
			fin_simindex.write(str(clstred_querygo) + "\t" + "\t".join(map(str,clstred_querygo_sim_list)) + "\n")
	fin_simindex.close()
	
	try:
		plt.switch_backend('agg')
		plt.rcParams["figure.figsize"] = [9.5,8] ## width, height
		plt.rcParams['axes.xmargin'] = 0
		plt.rcParams['axes.ymargin'] = 0
		sns.set(style = "white", font_scale = 2)
		sorted_go_sim_matrix = pd.read_csv(os.path.splitext(enGOclstr)[0] + ".simat", sep = "\t", header = 0, index_col = 0)
		sorted_go_sim_matrix = np.tril(sorted_go_sim_matrix, k = 0) # with elements above the k-th diagonal zeroed.
		sorted_go_sim_matrix = pd.DataFrame(sorted_go_sim_matrix).mask(sorted_go_sim_matrix < float(Ct))
		sim_heatmap = sns.heatmap(sorted_go_sim_matrix, vmin = float(Ct), vmax = 1.0, cmap = "Reds", xticklabels = False, yticklabels = False, cbar = True, cbar_kws = {'label': 'Similarity', 'ticks': [float(Ct), 1.0]})
#		colbar = sim_heatmap.figure.colorbar(sim_heatmap.collections[0])
#		colbar.set_ticks([float(Ct), 1.0])
		linepos = 0
		for clstrid in uniclstrids:
			linepos += clstrlist.count(clstrid)
			if linepos < len(clstrlist):
#				plt.axvline(x = linepos, color = '#3D3D3D', linestyle = '--', linewidth = 1.2)
				plt.axhline(y = linepos, color = '#3D3D3D', linestyle = '--', linewidth = 1.2)
		plt.axis('off')
		plt.savefig(os.path.splitext(enGOclstr)[0] + "_sim.png", dpi = 600, transparent = False, bbox_inches = 'tight', pad_inches = 0)
		plt.close()
	except ValueError as ErrorMessage:
		print("Error: Unable to generate a GO similarity heatmap due to:\n%s\n" % ErrorMessage)
	os.remove(str(os.path.splitext(enGOclstr)[0] + ".simat"))
	return None
		


def sim_density(simlist, SI, increment = 0.02):
#def sim_density(simlist, SI, Ct, increment = 0.02):
	try:
		plt.switch_backend('agg')
		plt.rcParams["figure.figsize"] = [8, 6.6] ## width, height
		plt.rcParams['axes.xmargin'] = 0
		plt.rcParams['axes.ymargin'] = 0
		sns.set(style = "white", font_scale = 2)
		binlist = np.arange(0.0, 1.02, float(increment))
		sim_density_plot = sns.distplot(simlist, hist = True, kde = False, norm_hist = False, bins = binlist, color = 'blue', hist_kws = {"edgecolor": "black", "alpha": 0.5}, kde_kws = {"color": "white", "alpha": 0})
		plt.xlim(-0.02, 1.02)
		plt.xlabel('Similarity Index (' + str(SI) + ')', fontsize = 18)
#		plt.xlabel("Similarity Index ({0}), P(>={1}): {2:.2%}".format(SI, Ct, len([sim for sim in simlist if float(Ct) >= 0.5])/float(len(simlist))), fontsize = 18)
		plt.ylabel('Frequency', fontsize = 18)
		return sim_density_plot
	except ValueError as ErrorMessage:
		print("Error: Unable to generate a GO similarity density plot for resulting clusters due to:\n%s\n" % ErrorMessage)

#def sim_cumulative(simlist, SI, increment = 0.02):
def sim_cumulative(simlist, SI, Ct, increment = 0.02):
	try:
		plt.switch_backend('agg')
		plt.rcParams["figure.figsize"] = [8, 6.6] ## width, height
		plt.rcParams['axes.xmargin'] = 0
		plt.rcParams['axes.ymargin'] = 0
		plt.rcParams["axes.labelsize"] = 22
		sns.set(style = "white", font_scale = 2)
		binlist = np.arange(0.0, 1.02, float(increment))
		sim_cumulative_plot = sns.distplot(simlist, hist = True, kde = False, norm_hist = True, bins = binlist, color = 'blue', hist_kws = {"edgecolor": "black", "alpha": 0.5, "cumulative": True}, kde_kws = {"cumulative": True, "linestyle": "--", "color": "gray"})
		sim_cumulative_plot.text(0.25, 0.92, "P(>={0}): {1:.2%}".format(Ct, len([sim for sim in simlist if float(sim) >= float(Ct)])/float(len(simlist))), fontsize = 20, color = 'black', ha = 'center', va = 'bottom')
#		plt.yticks(sim_cumulative_plot.get_yticks(), sim_cumulative_plot.get_yticks() * 100)
		sim_cumulative_plot.yaxis.set_major_formatter(FuncFormatter('{0:.0%}'.format))
		plt.xlim(-0.02, 1.02)
#		plt.ylim(0, 1.02)
		plt.xlabel('Similarity Index (' + str(SI) + ')', fontsize = 18)
#		plt.xlabel("Similarity Index ({0}), P(>={1}): {2:.2%}".format(SI, Ct, len([sim for sim in simlist if float(sim) >= float(sim)])/float(len(simlist))), fontsize = 18)
		plt.ylabel('Cumulative Frequency', fontsize = 18)
		return sim_cumulative_plot
	except ValueError as ErrorMessage:
		print("Error: Unable to generate a GO similarity density plot for resulting clusters due to:\n%s\n" % ErrorMessage)

		
def sim_newtork(enGOclstr, SI, Ct, Inf, sig):
	"""
	plot a similarity-based network for clustered GO
	enGOclstr	formatted GO file
	SI	similarity index between GOs, from ["JC","OC"]
	Ct	clustering threshold for the overlapping ratio between two GOs, any value between 0 and 1
	Inf	inflation value, main handle for cluster granularity, usually chosen somewhere in the range [1.2-5.0]
	sig signifance level (p-value cutoff) used in the enrichment test, e.g. 0.05
	"""
	enGOclstr_info_dict, gosim_dict = go_compare(enGOclstr, 2, 0, SI)
	clstredGOlist = rowtolist(enGOclstr, 2, "\t", "Y")
	uniclstrids = unielement(enGOclstr, 1, "\t", "Y")
	go_sim_newtork = nx.Graph()
	colorlist = ["#FF4136","#0074D9","#9F54E8","#F1C61C","#A5014F","#005884","#54A883","#6F7300","#FF851B","#00FF00"]
	
	with open(enGOclstr, "rU") as fin_clstredGO:
		for line_clstredGO in fin_clstredGO:
			line_clstredGO = line_clstredGO.rstrip("\r\n")
			if line_clstredGO.split("\t")[0].isdigit():			
				clstredGO_Clstr = line_clstredGO.split("\t")[0]
				clstredGO_GO = line_clstredGO.split("\t")[1]
				clstredGO_Description = line_clstredGO.split("\t")[2]
				clstredGO_Type = line_clstredGO.split("\t")[3]
				clstredGO_Depth = line_clstredGO.split("\t")[4]
				clstredGO_pvalue = line_clstredGO.split("\t")[5] 
				clstredGO_corrpvalue = line_clstredGO.split("\t")[6]
				clstredGO_xcatstest = line_clstredGO.split("\t")[7] 
				clstredGO_ncatsref = line_clstredGO.split("\t")[8] if int(line_clstredGO.split("\t")[8]) <> 0 else 100
				clstredGO_Xtotaltest = line_clstredGO.split("\t")[9]
				clstredGO_Ntotalref = line_clstredGO.split("\t")[10]
				clstredGO_NodeColor = colorlist[int(clstredGO_Clstr)-1] if int(clstredGO_Clstr) <= 10 else "#D3D3D3"
				clstredGO_NodeColorAlpha = colalpha(clstredGO_NodeColor, 1 - ((float(clstredGO_corrpvalue) - 0)/(float(sig) - 0)))
				go_sim_newtork.add_node(clstredGO_GO, Clstr = clstredGO_Clstr, Description = clstredGO_Description, Type = clstredGO_Type, Depth = clstredGO_Depth, corrpvalue = clstredGO_corrpvalue, xcatstest = clstredGO_xcatstest, ncatsref = clstredGO_ncatsref, Xtotaltest = clstredGO_Xtotaltest, Ntotalref = clstredGO_Ntotalref, NodeColor = clstredGO_NodeColor, NodeColorAlpha = clstredGO_NodeColorAlpha)
			
	for key_querygo in gosim_dict:
		for key_subjectgo in gosim_dict:
			if key_querygo <> key_subjectgo:
				if (key_querygo, key_subjectgo) not in list(go_sim_newtork.edges()) and (key_subjectgo, key_querygo) not in list(go_sim_newtork.edges()):
#					if max(gosim_dict[key_querygo][key_subjectgo], gosim_dict[key_subjectgo][key_querygo]) >= float(Ct):
#						go_sim_newtork.add_edge(key_querygo, key_subjectgo, weight = max(gosim_dict[key_querygo][key_subjectgo], gosim_dict[key_subjectgo][key_querygo]))
					if float(gosim_dict[key_querygo][key_subjectgo]) >= float(Ct):
						go_sim_newtork.add_edge(key_querygo, key_subjectgo, weight = float(gosim_dict[key_querygo][key_subjectgo]))

	try:
		plt.switch_backend('agg')
		plt.rcParams["figure.figsize"] = [10,10] 
		nx.draw(go_sim_newtork, pos = nx.spring_layout(go_sim_newtork, dim = 2, center = [0.5,0.5], k = 1.4 * 1/math.sqrt(len(enGOclstr_info_dict.keys())), iterations = 25, scale = 0.5, weight = "weight"), node_size = [66 * math.sqrt(int(nx.get_node_attributes(go_sim_newtork,'xcatstest')[node])) for node in go_sim_newtork.nodes()], node_color = [nx.get_node_attributes(go_sim_newtork,'NodeColorAlpha')[node] for node in go_sim_newtork.nodes()], edge_color = "#CCCCCC", width = [0.4 * float(edgeweight) for edgeweight in nx.get_edge_attributes(go_sim_newtork, "weight").values()], with_labels = False, labels = OrderedDict([(node, nx.get_node_attributes(go_sim_newtork,'Description')[node]) for node in go_sim_newtork.nodes()]), font_size = 10, font_color = "gray", font_family = "sans") ## node_size = [12 * math.sqrt(int(nx.get_node_attributes(go_sim_newtork,'xcatstest')[node])) for node in go_sim_newtork.nodes()]. "sans" is the default family in ggplot2.
		plt.ylim([-0.02, 1.02])
		plt.xlim([-0.02, 1.02])
		plt.axis('off')	
		plt.savefig(os.path.splitext(enGOclstr)[0] + "_network.png", dpi = 600, transparent = True, bbox_inches = 'tight', pad_inches = 0)
		plt.close()
	except ValueError as ErrorMessage:
		print("Error: GOMCL is unable to generate a graphical GO cluster network due to:\n%s\n" % ErrorMessage)
	return go_sim_newtork
	


def construct_go_hierarchy_subgraph(OBOInput, nodelist, sig = 0.05, **kwargs):
	'''
	OBOInput	obo file should be provided, e.g. go-basic.obo
	nodelist	nodes that should be used to construct the subgraph
	sig	signifance level (p-value cutoff) used in the enrichment test, e.g. 0.05
	**kwargs 	keyword arguments, specifically for GO information dictionary
	'''
	go_hierarchy_digraph = construct_go_hierarchy_digraph(OBOInput)

#	for nodes_in_subgraph in nodelist:
#		if nodes_in_subgraph not in go_hierarchy_digraph.nodes():
#			nodelist.remove(nodes_in_subgraph)
#			nodes_in_origraph_list = [node for node, attrs in go_hierarchy_digraph.nodes(data = True) if str(nodes_in_subgraph) in str(attrs["AltID"])]
#			nodes_in_origraph = nodes_in_origraph_list[0] if len(nodes_in_origraph_list) >=1 else ""
#			nodelist.append(nodes_in_origraph)
	
	go_hierarchy_subgraph = nx.DiGraph(go_hierarchy_digraph.subgraph(nodelist)) ##SubGraph Views are readonly, and has to be converted to a DiGraph object.
	
	def hierarchy_pos(G, root, rootx, rooty, width = 1.0, vert_gap = 0.2, pos = None):
		'''
		If the graph is a tree this will return the positions to plot this in a  hierarchical layout.
		G: the graph
		root: the root node of current branch 
		rootx: horizontal location of root, has to be a float
		rooty: vertical location of root, has to be a float
		width: horizontal space allocated for this branch - avoids overlap with other branches, has to be a float
		vert_gap: gap between levels of hierarchy, has to be a float
		pos: a dict saying where all nodes go if they have been assigned
		parent: parent of this branch. - only affects it if non-directed
		'''
		if pos is None:
			pos = {root:(rootx, rooty)}
		else:
			if root not in pos:
				pos[root] = (rootx, rooty)
#		children = list(G.neighbors(root))
		children = [child for child in list(G.neighbors(root)) if child not in pos and len(max(nx.all_simple_paths(G, source = root, target = child), key = len)) == 2]
#		children = [child for child in list(G.neighbors(root)) if child not in pos]
		if children:
			dx = width/len(children) 
			nextx = rootx - width/2 - dx/2
			for child in children:
				nextx += dx
				pos = hierarchy_pos(G, child, width = dx, rootx = nextx, rooty = rooty - random.uniform(vert_gap - 0.1, vert_gap + 0.1), vert_gap = vert_gap, pos = pos)
		return pos
	
	## If two nodes are connected in the orignal network, but not directly (meaning there are intermediate nodes, regardless of how many), these two nodes won't be connected in the subgraph. Need to identify these nodes and connect them.
	for nodeA_subgraph in nodelist:
		for nodeB_subgraph in nodelist:
			if nodeA_subgraph <> nodeB_subgraph:
				try:
					if nx.has_path(go_hierarchy_subgraph, source = nodeA_subgraph, target = nodeB_subgraph) is False and nx.has_path(go_hierarchy_subgraph, source = nodeB_subgraph, target = nodeA_subgraph) is False: ## Check if two nodes are not connected, no matter how many steps away
						if nx.has_path(go_hierarchy_digraph, source = nodeA_subgraph, target = nodeB_subgraph):
							go_hierarchy_subgraph.add_edge(nodeA_subgraph, nodeB_subgraph, weight = 0)
						elif nx.has_path(go_hierarchy_digraph, source = nodeB_subgraph, target = nodeA_subgraph):
							go_hierarchy_subgraph.add_edge(nodeB_subgraph, nodeA_subgraph, weight = 0)
				except nx.exception.NodeNotFound:
					ori_nodeA_subgraph_list = [node for node, attrs in go_hierarchy_digraph.nodes(data = True) if str(nodeA_subgraph) in str(attrs["AltID"])]
					ori_nodeB_subgraph_list = [node for node, attrs in go_hierarchy_digraph.nodes(data = True) if str(nodeB_subgraph) in str(attrs["AltID"])]
					ori_nodeA_subgraph = ori_nodeA_subgraph_list[0] if len(ori_nodeA_subgraph_list) >= 1 else ""
					ori_nodeB_subgraph = ori_nodeB_subgraph_list[0] if len(ori_nodeB_subgraph_list) >= 1 else ""
					try:
						if nx.has_path(go_hierarchy_digraph, source = nodeA_subgraph, target = ori_nodeB_subgraph):
							go_hierarchy_subgraph.add_node(ori_nodeB_subgraph, Name = nx.get_node_attributes(go_hierarchy_digraph,'Name')[ori_nodeB_subgraph], Namespace = nx.get_node_attributes(go_hierarchy_digraph,'Namespace')[ori_nodeB_subgraph], AltID = nx.get_node_attributes(go_hierarchy_digraph,'AltID')[ori_nodeB_subgraph], IsA = nx.get_node_attributes(go_hierarchy_digraph,'IsA')[ori_nodeB_subgraph], Relationship = nx.get_node_attributes(go_hierarchy_digraph,'Relationship')[ori_nodeB_subgraph], weight = 0)
							go_hierarchy_subgraph.add_edge(nodeA_subgraph, ori_nodeB_subgraph, weight = 0)
						elif nx.has_path(go_hierarchy_digraph, source = ori_nodeB_subgraph, target = nodeA_subgraph): 
							go_hierarchy_subgraph.add_edge(ori_nodeB_subgraph, nodeA_subgraph, weight = 0)
							go_hierarchy_subgraph.add_node(ori_nodeB_subgraph, Name = nx.get_node_attributes(go_hierarchy_digraph,'Name')[ori_nodeB_subgraph], Namespace = nx.get_node_attributes(go_hierarchy_digraph,'Namespace')[ori_nodeB_subgraph], AltID = nx.get_node_attributes(go_hierarchy_digraph,'AltID')[ori_nodeB_subgraph], IsA = nx.get_node_attributes(go_hierarchy_digraph,'IsA')[ori_nodeB_subgraph], Relationship = nx.get_node_attributes(go_hierarchy_digraph,'Relationship')[ori_nodeB_subgraph], weight = 0)
					except nx.exception.NodeNotFound:
						try:
							if nx.has_path(go_hierarchy_digraph, source = ori_nodeA_subgraph, target = nodeB_subgraph):
								go_hierarchy_subgraph.add_node(ori_nodeA_subgraph, Name = nx.get_node_attributes(go_hierarchy_digraph,'Name')[ori_nodeA_subgraph], Namespace = nx.get_node_attributes(go_hierarchy_digraph,'Namespace')[ori_nodeA_subgraph], AltID = nx.get_node_attributes(go_hierarchy_digraph,'AltID')[ori_nodeA_subgraph], IsA = nx.get_node_attributes(go_hierarchy_digraph,'IsA')[ori_nodeA_subgraph], Relationship = nx.get_node_attributes(go_hierarchy_digraph,'Relationship')[ori_nodeA_subgraph], weight = 0)
								go_hierarchy_subgraph.add_edge(ori_nodeA_subgraph, nodeB_subgraph, weight = 0)
							elif nx.has_path(go_hierarchy_digraph, source = nodeB_subgraph, target = ori_nodeA_subgraph): 
								go_hierarchy_subgraph.add_node(ori_nodeA_subgraph, Name = nx.get_node_attributes(go_hierarchy_digraph,'Name')[ori_nodeA_subgraph], Namespace = nx.get_node_attributes(go_hierarchy_digraph,'Namespace')[ori_nodeA_subgraph], AltID = nx.get_node_attributes(go_hierarchy_digraph,'AltID')[ori_nodeA_subgraph], IsA = nx.get_node_attributes(go_hierarchy_digraph,'IsA')[ori_nodeA_subgraph], Relationship = nx.get_node_attributes(go_hierarchy_digraph,'Relationship')[ori_nodeA_subgraph], weight = 0)
								go_hierarchy_subgraph.add_edge(nodeB_subgraph, ori_nodeA_subgraph, weight = 0)
						except nx.exception.NodeNotFound:
							try:
								if nx.has_path(go_hierarchy_digraph, source = ori_nodeA_subgraph, target = ori_nodeB_subgraph):
									go_hierarchy_subgraph.add_node(ori_nodeA_subgraph, Name = nx.get_node_attributes(go_hierarchy_digraph,'Name')[ori_nodeA_subgraph], Namespace = nx.get_node_attributes(go_hierarchy_digraph,'Namespace')[ori_nodeA_subgraph], AltID = nx.get_node_attributes(go_hierarchy_digraph,'AltID')[ori_nodeA_subgraph], IsA = nx.get_node_attributes(go_hierarchy_digraph,'IsA')[ori_nodeA_subgraph], Relationship = nx.get_node_attributes(go_hierarchy_digraph,'Relationship')[ori_nodeA_subgraph], weight = 0)
									go_hierarchy_subgraph.add_node(ori_nodeB_subgraph, Name = nx.get_node_attributes(go_hierarchy_digraph,'Name')[ori_nodeB_subgraph], Namespace = nx.get_node_attributes(go_hierarchy_digraph,'Namespace')[ori_nodeB_subgraph], AltID = nx.get_node_attributes(go_hierarchy_digraph,'AltID')[ori_nodeB_subgraph], IsA = nx.get_node_attributes(go_hierarchy_digraph,'IsA')[ori_nodeB_subgraph], Relationship = nx.get_node_attributes(go_hierarchy_digraph,'Relationship')[ori_nodeB_subgraph], weight = 0)
									go_hierarchy_subgraph.add_edge(ori_nodeA_subgraph, ori_nodeB_subgraph, weight = 0)
								elif nx.has_path(go_hierarchy_digraph, source = ori_nodeB_subgraph, target = ori_nodeA_subgraph):
									go_hierarchy_subgraph.add_node(ori_nodeA_subgraph, Name = nx.get_node_attributes(go_hierarchy_digraph,'Name')[ori_nodeA_subgraph], Namespace = nx.get_node_attributes(go_hierarchy_digraph,'Namespace')[ori_nodeA_subgraph], AltID = nx.get_node_attributes(go_hierarchy_digraph,'AltID')[ori_nodeA_subgraph], IsA = nx.get_node_attributes(go_hierarchy_digraph,'IsA')[ori_nodeA_subgraph], Relationship = nx.get_node_attributes(go_hierarchy_digraph,'Relationship')[ori_nodeA_subgraph], weight = 0)
									go_hierarchy_subgraph.add_node(ori_nodeB_subgraph, Name = nx.get_node_attributes(go_hierarchy_digraph,'Name')[ori_nodeB_subgraph], Namespace = nx.get_node_attributes(go_hierarchy_digraph,'Namespace')[ori_nodeB_subgraph], AltID = nx.get_node_attributes(go_hierarchy_digraph,'AltID')[ori_nodeB_subgraph], IsA = nx.get_node_attributes(go_hierarchy_digraph,'IsA')[ori_nodeB_subgraph], Relationship = nx.get_node_attributes(go_hierarchy_digraph,'Relationship')[ori_nodeB_subgraph], weight = 0)
									go_hierarchy_subgraph.add_edge(ori_nodeB_subgraph, ori_nodeA_subgraph, weight = 0)
							except nx.exception.NodeNotFound:
								print("Either "+ str(nodeA_subgraph) + " or " + str(nodeB_subgraph) + " is not in GO hierarchy")
								pass
							pass
						pass
					pass
					
	edges_present = [(u,v) for (u,v,d) in go_hierarchy_subgraph.edges(data = True) if d['weight'] == 1] ## Edges present in the subgraph.
	edges_absent = [(u,v) for (u,v,d) in go_hierarchy_subgraph.edges(data = True) if d['weight'] == 0] ## Edges present in the orignal network, but absent from the subgraph.

	if "GOinfodict" in kwargs:
		for nodes_in_subgraph in go_hierarchy_subgraph.nodes():
			try:
				if nodes_in_subgraph in kwargs["GOinfodict"]:
					node_attributes = kwargs["GOinfodict"][nodes_in_subgraph].strip()
				else:
					for nodes_in_subgraph_altID in go_hierarchy_digraph.nodes(data = True)[nodes_in_subgraph]["AltID"]:
						if nodes_in_subgraph_altID in kwargs["GOinfodict"]:
							node_attributes = kwargs["GOinfodict"][nodes_in_subgraph_altID].strip()
			except Exception as exceptions:
				print("Couldn't generate a GO hierarchy graph for each cluster due to: \n%s" % exceptions)
			go_hierarchy_subgraph.nodes()[nodes_in_subgraph]["Depth"] = node_attributes.split("\t")[3]
			go_hierarchy_subgraph.nodes()[nodes_in_subgraph]["corrpvalue"] = node_attributes.split("\t")[5]
			go_hierarchy_subgraph.nodes()[nodes_in_subgraph]["xcatstest"] = node_attributes.split("\t")[6] 
			go_hierarchy_subgraph.nodes()[nodes_in_subgraph]["ncatsref"] = node_attributes.split("\t")[7] if int(node_attributes.split("\t")[7]) <> 0 else 100
			go_hierarchy_subgraph.nodes()[nodes_in_subgraph]["Xtotaltest"] = node_attributes.split("\t")[8]
			go_hierarchy_subgraph.nodes()[nodes_in_subgraph]["Ntotalref"] = node_attributes.split("\t")[9]
#			go_hierarchy_subgraph.nodes()[nodes_in_subgraph]["Nodecolor"] = "#FF8C00"
			go_hierarchy_subgraph.nodes()[nodes_in_subgraph]["NodeColorAlpha"] = colalpha("#FF8C00", 1 - ((float(node_attributes.split("\t")[5]) - 0)/(float(sig) - 0)))
	else:
		nx.set_node_attributes(go_hierarchy_subgraph, "", "Depth")
		nx.set_node_attributes(go_hierarchy_subgraph, 0, "corrpvalue")
		nx.set_node_attributes(go_hierarchy_subgraph, 100, "xcatstest")
		nx.set_node_attributes(go_hierarchy_subgraph, 100, "ncatsref")
		nx.set_node_attributes(go_hierarchy_subgraph, 100, "Xtotaltest")
		nx.set_node_attributes(go_hierarchy_subgraph, 100, "Ntotalref")
		nx.set_node_attributes(go_hierarchy_subgraph, "#FF8C00", "NodeColorAlpha")
	
	
	singletlist = []
	rootlist = []
	bottomchildlist = []
	try:
		if isinstance(go_hierarchy_subgraph, nx.DiGraph):
			singletlist.extend(list(nx.isolates(go_hierarchy_subgraph))) ## Capture all singlets.	
#			root = next(iter(nx.topological_sort(go_hierarchy_subgraph)))
#			topological_list = list(nx.topological_sort(go_hierarchy_subgraph))
			for node in list(go_hierarchy_subgraph.nodes()):
				if not list(go_hierarchy_subgraph.predecessors(node)): ## If empty, no predecessors. In Python, if a variable is a numeric zero or empty, or a None object then it is considered as False, otherwise True.
					if list(go_hierarchy_subgraph.neighbors(node)): ## If any successors. G.neighbors() in a DiGraph is equivalent to G.successors(). 
						rootlist.append(node)
				else:
					if not list(go_hierarchy_subgraph.neighbors(node)):
						bottomchildlist.append(node)
	except Exception as exceptions:
		print("Couldn't generate a hierarchical layout due to: \n%s\n" % exceptions)

	try:
		plt.switch_backend('agg')
#		plt.rcParams["figure.figsize"] = [10,10] 
#		Create hierarchical positions:

		if rootlist:
			pos = {}
			unassignedy = 2.0
			rootlist = sorted(rootlist, key = lambda rootnode : (int(go_hierarchy_subgraph.nodes()[rootnode]["xcatstest"]), len(list(go_hierarchy_subgraph.neighbors(rootnode)))), reverse = True)
			for root_num in range(len(rootlist)):
				rootx = float(root_num) * len(rootlist) * 30
				rooty = random.uniform(0.8, 1.6) 
				pos = hierarchy_pos(go_hierarchy_subgraph, rootlist[root_num], width = float((len(rootlist) - 1) * len(rootlist) * 10), rootx = rootx, rooty = rooty, vert_gap = 0.8, pos = pos)
			if singletlist: ## If not empty
				singletlist = sorted(singletlist, key = lambda singletnode : (len(list(go_hierarchy_subgraph.neighbors(singletnode))),int(go_hierarchy_subgraph.nodes()[singletnode]["xcatstest"])), reverse = True)
				for singlet_num in range(len(singletlist)):
					pos[singletlist[singlet_num]] = (singlet_num * ((len(rootlist) - 1) * len(rootlist)/len(singletlist)) * 30, random.uniform(1.8, 2.2))
			for node in list(go_hierarchy_subgraph.nodes()):
				if node not in pos:
					unassignedy -= 0.1
					pos[node] = (0.1, unassignedy)
#			pos = {u:(r * math.cos(theta),r * math.sin(theta)) for u, (theta, r) in pos.items()}
		else:
			pos = nx.spring_layout(go_hierarchy_subgraph, dim = 2, center = [0.5,0.5], k = 1.4 * 1/math.sqrt(len(nodelist)), iterations = 25, scale = 0.5)
		
		hierarchy_subgraph = nx.draw_networkx_nodes(go_hierarchy_subgraph, pos, node_size = [50 * math.sqrt(int(nx.get_node_attributes(go_hierarchy_subgraph,'xcatstest')[node])) for node in go_hierarchy_subgraph.nodes()], node_color = [nx.get_node_attributes(go_hierarchy_subgraph,'NodeColorAlpha')[node] for node in go_hierarchy_subgraph.nodes()], alpha = 1.0, edgecolors = "#696969", linewidths = 0.8) # node_size = [50 * (int(nx.get_node_attributes(go_hierarchy_subgraph,'xcatstest')[node]) ** (1.0/3.0)) for node in go_hierarchy_subgraph.nodes()]
#		nx.draw_networkx_labels(go_hierarchy_subgraph, pos, labels = OrderedDict([(node, nx.get_node_attributes(go_hierarchy_subgraph,'Name')[node]) for node in go_hierarchy_subgraph.nodes() if node in rootlist or node in singletlist]), font_size = 6.5, font_weight = 'medium', font_color = "black", font_family = "sans-serif")
		nx.draw_networkx_edges(go_hierarchy_subgraph, pos, edgelist = edges_present, arrows = True, width = 0.6, edge_color = '#000000', alpha = 0.6) # #808080 ##arrowsize = 8, arrowstyle = '->',
		nx.draw_networkx_edges(go_hierarchy_subgraph, pos, edgelist = edges_absent, arrows = True, width = 0.6, edge_color = '#D3D3D3', alpha = 0.8, style = 'dashed') ## style = "dashed" doesn't seem to work. ##arrowsize = 8, arrowstyle = '->', 
		xmin, xmax, ymin, ymax = plt.axis()
		plt.xlim([xmin - (xmax - xmin) * 0.1, xmax + (xmax - xmin) * 0.1])
		plt.ylim([ymin - (ymax - ymin) * 0.1, ymax + (ymax - ymin) * 0.1])
		plt.axis('off')	
	except ValueError as ErrorMessage:
		print("Error: GOMCL is unable to generate a graphical GO cluster network due to:\n%s\n" % ErrorMessage)
		
	return go_hierarchy_subgraph, hierarchy_subgraph



synopsis = "\n\
#############################################################################################################################################\n\
#go_clustering.py clusters GO terms using MCL based on overlapping ratios, OC (Overlap coefficient) or JC (Jaccard coefficient).            #\n\
#Use examples:                                                                                                                              #\n\
#   go_clustering.py -Ct 0.5 -I 1.5 EnrichedGO.txt                                                                                          #\n\
#   go_clustering.py -SI JC -Ct 0.5 -I 1.5 EnrichedGO.txt                                                                                   #\n\
#                                                                                                                                           #\n\
#############################################################################################################################################"


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = synopsis, formatter_class = argparse.RawTextHelpFormatter)
	parser.add_argument("enGOclstr", metavar = "-I_ClstredGO", help = "Formatted enriched GO input file", action = "store", nargs = None, const = None, default = None, type = None, choices = None) ## See below.
	parser.add_argument("-SI", metavar = None, dest = "simindex", help = "Method to calculate similarity between GO terms, OC (Overlap coefficient) or JC (Jaccard coefficient) (default: %(default)s)", action = "store", nargs = None, const = None, default = "OC", type = str, choices = ["OC", "JC"]) 
	parser.add_argument("-Ct", metavar = None, dest = "cutoff", help = "Clustering threshold for the overlapping ratio between two GO terms, any value between 0 and 1 (default: %(default)s)", action = "store", nargs = None, const = None, default = "0.5", type = float, choices = None) 
	parser.add_argument("-I", metavar = None, dest = "inflation", help = "Inflation value, main handle for cluster granularity, usually chosen somewhere in the range [1.2-5.0] (default: %(default)s)", action = "store", nargs = None, const = None, default = "2.0", type = float, choices = None) ## -I 5.0 will tend to result in fine-grained clusterings, and -I 1.2 will tend to result in very coarse grained clusterings.
	parser.add_argument("-SL", metavar = None, dest = "sig", help = "Signifance level (p-value cutoff) used in the enrichment test (default: %(default)s)", action = "store", nargs = None, const = None, default = "0.05", type = float, choices = None)
	parser.add_argument("-hm",  dest = None, help = "Only needed when a similarity heatmap is desired", action = "store_true", default = None ) ## Argument present --> true, not present --> false.
	parser.add_argument("-nw",  dest = None, help = "Only needed when a similarity-based network is desired", action = "store_true", default = None ) ## Argument present --> true, not present --> false.
	parser.add_argument("-d", dest = "dswitch", help = "Only needed if depth for input GO terms is desired", action = "store_true", default = None)
	args = parser.parse_args() 
	
	if args.hm:
		sim_plot(args.enGOclstr, args.simindex, args.cutoff, args.inflation)
	if args.nw:
		go_sim_newtork = sim_newtork(args.enGOclstr, args.simindex, args.cutoff, args.inflation, args.sig)