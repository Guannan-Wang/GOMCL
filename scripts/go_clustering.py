#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__    = "Guannan Wang"
__email__     = "wangguannan2014@gmail.com"
__version__   = "1"

import argparse, sys, re, os, io, operator, subprocess
from decimal import Decimal ## This is for float comparison when accuracy is important!
from collections import OrderedDict 
import numpy as np
import pandas as pd
import math
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

def go_compare(enGOfmtfltr, enGOcol, genecol, SI):
	"""
	pairwise comparison for all GOs in the formatted GO list
	enGOfmtfltr	formatted GO file
	enGOcol	column number where GO ID is, should be integer
	genecol	column number where genes in each GO ID are, should be integer
	SI	similarity index, from ["JC","OC"]
	"""
	fin_enGOfmtfltr = open(enGOfmtfltr, "rU")
	oriGOlist = rowtolist(enGOfmtfltr, enGOcol, "\t", "Y")
	enGOfmtfltr_info_dict = dict()
	gosim_dict = dict()
	for line_enGOfmtfltr in fin_enGOfmtfltr:
		line_enGOfmtfltr = line_enGOfmtfltr.rstrip("\r\n")
#		if line_enGOfmtfltr.startswith("GO:") and line_enGOfmtfltr.split("\t")[0].split("GO:")[1].isdigit():
		if line_enGOfmtfltr.split("\t")[int(enGOcol)-1].startswith("GO:") and line_enGOfmtfltr.split("\t")[int(enGOcol)-1].split("GO:")[1].isdigit():
			enGOfmtfltr_ID = line_enGOfmtfltr.split("\t")[int(enGOcol)-1]
			enGOfmtfltr_info_dict[enGOfmtfltr_ID] = line_enGOfmtfltr
			gosim_dict[enGOfmtfltr_ID] = dict()
	for element_queryGO in oriGOlist:
		for element_subjectGO in oriGOlist:
			if SI == "OC":
				## A∩B/A,A∩B/B 
				gosim_dict[element_queryGO][element_subjectGO] = float(len(intersect(enGOfmtfltr_info_dict[element_subjectGO].split("\t")[int(genecol)-1].split("|"),enGOfmtfltr_info_dict[element_queryGO].split("\t")[int(genecol)-1].split("|"))))/float(len(enGOfmtfltr_info_dict[element_queryGO].split("\t")[int(genecol)-1].split("|")))
			if SI == "JC":
				## A∩B/A⋃B 
				gosim_dict[element_queryGO][element_subjectGO] = float(len(intersect(enGOfmtfltr_info_dict[element_subjectGO].split("\t")[int(genecol)-1].split("|"),enGOfmtfltr_info_dict[element_queryGO].split("\t")[int(genecol)-1].split("|"))))/float(len(set(enGOfmtfltr_info_dict[element_queryGO].split("\t")[int(genecol)-1].split("|") + enGOfmtfltr_info_dict[element_subjectGO].split("\t")[int(genecol)-1].split("|"))))
	return enGOfmtfltr_info_dict, gosim_dict

def go_clustering(enGOfmtfltr, SI, Ct, Inf):
	## This part will need to be updated to use markov clustering in python.
	"""
	markove clustering of all GOs in the formatted GO list
	enGOfmtfltr	formatted GO file
	SI	similarity index between GOs, from ["JC","OC"]
	Ct	clustering threshold for the overlapping ratio between two GOs, any value between 0 and 1
	Inf	inflation value, main handle for cluster granularity, usually chosen somewhere in the range [1.2-5.0]
	will think about using mcl clustering python module (it requires python 3+).
	"""
	enGOfmtfltr_info_dict, gosim_dict = go_compare(enGOfmtfltr, 1, 0, SI) ## This is for formatted or filtered GO file where the first column is the GO ID and the last column is the genes, pay attention to the use of "0" here.
	oriGOlist = rowtolist(enGOfmtfltr, 1, "\t", "Y")
	
	fin_enGOcmped = open(os.path.splitext(enGOfmtfltr)[0] + "_Ct" + str(Ct) + "I" + str(Inf) + ".enGOcmped.temp", "w")
	ProcessedPairs = []
	for element_GO_A in oriGOlist:
		for element_GO_B in oriGOlist:
			GO_PairAB = str(element_GO_A) + "_" + str(element_GO_B)
			GO_PairBA = str(element_GO_B) + "_" + str(element_GO_A)
			if GO_PairAB not in ProcessedPairs:
				if max(gosim_dict[element_GO_A][element_GO_B],gosim_dict[element_GO_B][element_GO_A]) >= float(Ct):
					fin_enGOcmped.write(str(element_GO_A) + "\t" + str(element_GO_B) + "\t" + str(max(gosim_dict[element_GO_A][element_GO_B],gosim_dict[element_GO_B][element_GO_A])) + "\n")
					ProcessedPairs.extend([GO_PairAB,GO_PairBA])
	fin_enGOcmped.close()
	
	cmd_mcl = "mcl " + str(os.path.splitext(enGOfmtfltr)[0] + "_Ct" + str(Ct) + "I" + str(Inf) + ".enGOcmped.temp") + " --abc -I " + str(Inf) + " -o " + str(os.path.splitext(enGOfmtfltr)[0] + "_Ct" + str(Ct) + "I" + str(Inf) + ".enGOcmped.mcl.temp")
	os.system(cmd_mcl)
	
#	mcl_std = subprocess.Popen(["mcl", str(os.path.splitext(enGOfmtfltr)[0] + ".enGOcmped.temp"), "--abc", "-I", str(Inf), "-o", str(os.path.splitext(enGOfmtfltr)[0] + ".enGOcmped.mcl.temp")],  stdout = subprocess.PIPE,  stderr = subprocess.STDOUT)
#	stdout,stderr = mcl_std.communicate()
	
	fin_mcl = open(os.path.splitext(enGOfmtfltr)[0] + "_Ct" + str(Ct) + "I" + str(Inf) + ".enGOcmped.mcl.temp","rU")
	mclgroups = fin_mcl.readlines()
	go_mclgroup_dict = dict()
	gene_mclgroup_dict = dict()
	for groupnum in range(len(mclgroups)):
		mclgroup = mclgroups[groupnum].strip()
		mclgroupid = str("mcl" + str(groupnum + 1))
		go_mclgroup_dict[mclgroupid] = mclgroup
		gene_mclgroup_dict[mclgroupid] = []
		for go_in_mclgroup in mclgroup.split("\t"):
			gene_mclgroup_dict[mclgroupid].extend(enGOfmtfltr_info_dict[go_in_mclgroup].split("\t")[-1].split("|"))
	fin_mcl.close()
	
	os.remove(str(os.path.splitext(enGOfmtfltr)[0] + "_Ct" + str(Ct) + "I" + str(Inf) + ".enGOcmped.temp"))
	os.remove(str(os.path.splitext(enGOfmtfltr)[0] + "_Ct" + str(Ct) + "I" + str(Inf) + ".enGOcmped.mcl.temp"))
	return go_mclgroup_dict, gene_mclgroup_dict

def go_assign_cluster(enGOfmtfltr, SI, Ct, Inf):
	"""
	assign GOs to clusters
	enGOfmtfltr	formatted GO file
	SI	similarity index between GOs, from ["JC","OC"]
	Ct	clustering threshold for the overlapping ratio between two GOs, any value between 0 and 1
	Inf	inflation value, main handle for cluster granularity, usually chosen somewhere in the range [1.2-5.0]
	"""
#	enGOfmtfltr_info_dict, gosim_dict = go_compare(enGOfmtfltr, 1, 0, SI)
	go_mclgroup_dict, gene_mclgroup_dict = go_clustering(enGOfmtfltr, SI, Ct, Inf)
	
	go_clstr_dict = dict()
	gene_clstr_dict = dict()
	clstred_go_dict = dict()
	clstred_gene_dict = dict()
	clstrnum = 1
	for mclgroupid in sorted(gene_mclgroup_dict, key = lambda key_mclgroupid : len(set(gene_mclgroup_dict[key_mclgroupid])), reverse = True):
		clstred_go_dict[clstrnum] = go_mclgroup_dict[mclgroupid].split("\t")
		clstred_gene_dict[clstrnum] = gene_mclgroup_dict[mclgroupid]
		for go_in_mclgroup in go_mclgroup_dict[mclgroupid].split("\t"):
			go_clstr_dict[go_in_mclgroup] = clstrnum
		for gene_in_mclgroup in set(gene_mclgroup_dict[mclgroupid]):
			gene_clstr_dict[gene_in_mclgroup] = clstrnum
		clstrnum += 1
	return go_clstr_dict, gene_clstr_dict, clstred_go_dict, clstred_gene_dict
	


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
	parser.add_argument("enGOfmtfltr", metavar = "-I_EnrichedGO", help = "Formatted enriched GO input file", action = "store", nargs = None, const = None, default = None, type = None, choices = None) ## See below.
	parser.add_argument("-SI", metavar = None, dest = "simindex", help = "Method to calculate similarity between GO terms, OC (Overlap coefficient) or JC (Jaccard coefficient) (default: %(default)s)", action = "store", nargs = None, const = None, default = "OC", type = str, choices = ["OC", "JC"]) 
	parser.add_argument("-Ct", metavar = None, dest = "cutoff", help = "Clustering threshold for the overlapping ratio between two GO terms, any value between 0 and 1 (default: %(default)s)", action = "store", nargs = None, const = None, default = "0.5", type = float, choices = None) 
	parser.add_argument("-I", metavar = None, dest = "inflation", help = "Inflation value, main handle for cluster granularity, usually chosen somewhere in the range [1.2-5.0] (default: %(default)s)", action = "store", nargs = None, const = None, default = "2.0", type = float, choices = None) ## -I 5.0 will tend to result in fine-grained clusterings, and -I 1.2 will tend to result in very coarse grained clusterings. 
#	parser.add_argument("-N", metavar = None, dest = "Name", help = "Unique name append to cluster numbers, e.g. Name_1, this is especially important if visulization in Cytoscape is desired", action = "store", nargs = None, const = None, default = "", type = str, choices = None)
	args = parser.parse_args() 
	
	go_clstr_dict, gene_clstr_dict, clstred_go_dict, clstred_gene_dict = go_assign_cluster(args.enGOfmtfltr, args.simindex, args.cutoff, args.inflation)

	with open(os.path.splitext(args.enGOfmtfltr)[0] + ".clstrinfo", "w") as fout_clstrinfo:
		fout_clstrinfo.write("GO.Clstr" + "\t" + "# of GOs" + "\t" + "# of genes" + "\n")
		for clstrid in clstred_go_dict:
			fout_clstrinfo.write("Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_" + str(clstrid) + "\t" + str(len(set(clstred_go_dict[clstrid]))) + "\t" + str(len(set(clstred_gene_dict[clstrid]))) + "\n")	
