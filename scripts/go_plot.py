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
from go_clustering import *
from funs import *

try:
	import networkx as nx
except ImportError as ErrorMessage:
	print(str(ErrorMessage))
	sys.exit()

try:	
	import matplotlib.pyplot as plt
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
		sorted_go_sim_matrix = sorted_go_sim_matrix.mask(sorted_go_sim_matrix < float(Ct))
		sns.heatmap(sorted_go_sim_matrix, vmin = float(Ct), vmax = 1.0, cmap = "Reds", xticklabels = False, yticklabels = False, cbar_kws = {'label': 'Similarity'})
		linepos = 0
		for clstrid in uniclstrids:
			linepos += clstrlist.count(clstrid)
			if linepos < len(clstrlist):
##				plt.axvline(x = linepos, color = '#3D3D3D', linestyle = '--', linewidth = 1.2)
				plt.axhline(y = linepos - 1, color = '#3D3D3D', linestyle = '--', linewidth = 1.2)
		plt.axis('off')
		plt.savefig(os.path.splitext(enGOclstr)[0] + "_sim.png", dpi = 600, transparent = True)
		plt.close()
	except ValueError as ErrorMessage:
		print("Error: Unable to generate a GO similarity heatmap due to:\n%s\n" % ErrorMessage)
	os.remove(str(os.path.splitext(enGOclstr)[0] + ".simat"))
	return None
		
		
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
				clstredGO_ncatsref = line_clstredGO.split("\t")[8]
				clstredGO_Xtotaltest = line_clstredGO.split("\t")[9]
				clstredGO_Ntotalref = line_clstredGO.split("\t")[10]
				clstredGO_NodeColor = colorlist[int(clstredGO_Clstr)-1] if int(clstredGO_Clstr) <= 10 else "#D3D3D3"
				clstredGO_NodeColorAlpha = colalpha(clstredGO_NodeColor, 1 - ((float(clstredGO_corrpvalue) - 0)/(float(sig) - 0)))
				go_sim_newtork.add_node(clstredGO_GO, Clstr = clstredGO_Clstr, Description = clstredGO_Description, Type = clstredGO_Type, Depth = clstredGO_Depth, corrpvalue = clstredGO_corrpvalue, xcatstest = clstredGO_xcatstest, ncatsref = clstredGO_ncatsref, Xtotaltest = clstredGO_Xtotaltest, Ntotalref = clstredGO_Ntotalref, NodeColor = clstredGO_NodeColor, NodeColorAlpha = clstredGO_NodeColorAlpha)
			
	for key_querygo in gosim_dict:
		for key_subjectgo in gosim_dict:
			if key_querygo <> key_subjectgo:
				if (key_querygo, key_subjectgo) not in list(go_sim_newtork.edges()) and (key_subjectgo, key_querygo) not in list(go_sim_newtork.edges()):
					if max(gosim_dict[key_querygo][key_subjectgo], gosim_dict[key_subjectgo][key_querygo]) >= float(Ct):
						go_sim_newtork.add_edge(key_querygo, key_subjectgo, weight = max(gosim_dict[key_querygo][key_subjectgo], gosim_dict[key_subjectgo][key_querygo]))
	
	try:
		plt.switch_backend('agg')
		plt.rcParams["figure.figsize"] = [10,10]
		nx.draw(go_sim_newtork, pos = nx.spring_layout(go_sim_newtork, dim = 2, center = [0.5,0.5], k = 1.4 * 1/math.sqrt(len(enGOclstr_info_dict.keys())), iterations = 25, scale = 0.5, weight = "weight"), node_size = [12 * math.sqrt(int(nx.get_node_attributes(go_sim_newtork,'ncatsref')[node])) for node in go_sim_newtork.nodes()], node_color = [nx.get_node_attributes(go_sim_newtork,'NodeColorAlpha')[node] for node in go_sim_newtork.nodes()], edge_color = "#CCCCCC", width = [0.4 * float(edgeweight) for edgeweight in nx.get_edge_attributes(go_sim_newtork, "weight").values()], with_labels = False, labels = OrderedDict([(node, nx.get_node_attributes(go_sim_newtork,'Description')[node]) for node in go_sim_newtork.nodes()]), font_size = 10, font_color = "gray", font_family = "sans") ## "sans" is the default family in ggplot2
		plt.ylim([-0.02, 1.02])
		plt.xlim([-0.02, 1.02])
		plt.axis('off')	
		plt.savefig(os.path.splitext(enGOclstr)[0] + "_network.png", dpi = 600)
		plt.close()
	except ValueError as ErrorMessage:
		print("Error: GOMCL is unable to generate a graphical GO cluster network due to:\n%s\n" % ErrorMessage)
	return go_sim_newtork
	


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