#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__    = "Guannan Wang"
__email__     = "wangguannan2014@gmail.com"
__version__   = "1"

import argparse, sys, re, os, io, operator 
from decimal import Decimal ## This is for float comparison when accuracy is important!

try:
	import networkx as nx
except ImportError as ErrorMessage:
	print(str(ErrorMessage))
	sys.exit()

from .go_obo_parser import obo_parser, construct_go_hierarchy_digraph
from .funs import *

synopsis = "\n\
#############################################################################################################################################\n\
#go_enrichment_result_formatter.py formats output from different GO enrichment analysis tools (e.g. BiNGO, agriGO, AmiGO, etc.) to a        #\n\
#standard format, and sort either by different columns.                                                                                     #\n\
#                                                                                                                                           #\n\
#############################################################################################################################################"


def goea_formatter(OBOInput, goeatool, enGOraw, dswitch = False):
	"""
	gene ontology enrichment analysis (goea) result formatter, supporting BiNGO, agriGO, AmiGO, PANTHER, GOrilla, gProfiler
	Required:
	OBOInput	obo file should be provided, e.g. go-basic.obo
	goeatool	the go enrichment tools used: [BiNGO, agriGO, AmiGO, PANTHER,GOrilla, gProfiler, Enrichr]
	enGOraw	the enrichment test reults from the go enrichment tool used
	dswitch switch for depth ON and OFF
	"""
	obo_go_name_dict, obo_go_namespace_dict, obo_go_alt_id_dict, obo_go_is_obsolete_dict, obo_go_is_a_dict, obo_go_relationship_dict = obo_parser(OBOInput)
	go_hierarchy_network = construct_go_hierarchy_digraph(OBOInput)
	
	RootTerm_dict = {"biological_process":"GO:0008150", "molecular_function":"GO:0003674", "cellular_component":"GO:0005575"}
	
	fin_enGOraw = open(enGOraw, "rU")
#	fout_enGOfmtted = open(os.path.splitext(enGOraw)[0] + ".enGOfmtted", "w")
#	fout_enGOfmtted.write("Full GO-ID" + "\t" + "Description" + "\t" + "Type" + "\t" + "Depth" + "\t" + "p-value" + "\t" + "adj p-value" + "\t" + "x.cats.test" + "\t" + "n.cats.ref" + "\t" + "X.total.test" + "\t" + "N.total.ref" + "\t" + "Genes in test set" + "\n")
	enGOfmtted_list = []
	numoflines = 0
	for line_enGOraw in fin_enGOraw:
		line_enGOraw = line_enGOraw.rstrip("\r\n")
		numoflines += 1
		if numoflines % 50 == 0:
			sys.stdout.write("Done formatting %d GOs\r" % numoflines)
			sys.stdout.flush()
		
		## Read column information from enrichment results
		if str(goeatool).lower() == "bingo":
			if len(line_enGOraw.split("\t")) > 0 and isint(line_enGOraw.split("\t")[0]):
				goea_go_id = "GO:%07d" % int(line_enGOraw.split("\t")[0])
				goea_go_description = str(line_enGOraw.split("\t")[7])
				goea_go_pvalue = str(line_enGOraw.split("\t")[1])
				goea_go_adj_pvalue = str(line_enGOraw.split("\t")[2])
				goea_go_cats_test = int(line_enGOraw.split("\t")[3])
				goea_go_cats_ref = int(line_enGOraw.split("\t")[4])
				goea_go_total_test = int(line_enGOraw.split("\t")[5])
				goea_go_total_ref = int(line_enGOraw.split("\t")[6])
				goea_go_geneset = str(line_enGOraw.split("\t")[8])
			else:
				continue
		elif str(goeatool).lower() == "gprofiler":
			if line_enGOraw.startswith("GO:"):
				goea_go_id = str(line_enGOraw.split("\t")[2])
				goea_go_description = str(line_enGOraw.split("\t")[1])
				goea_go_pvalue = "na"
				goea_go_adj_pvalue = str(line_enGOraw.split("\t")[3])
				goea_go_cats_test = int(line_enGOraw.split("\t")[7])
				goea_go_cats_ref = int(line_enGOraw.split("\t")[5])
				goea_go_total_test = int(line_enGOraw.split("\t")[6])
				goea_go_total_ref = int(line_enGOraw.split("\t")[8])
				goea_go_geneset = str(line_enGOraw.split("\t")[9].replace(",","|"))
			else:
				continue
		elif str(goeatool).lower() == "gorilla":
#			if isint(line_enGOraw.split("\t")[0].split("GO:")[-1]):
			if line_enGOraw.startswith("GO:"):
				goea_go_id = str(line_enGOraw.split("\t")[0])
				goea_go_description = str(line_enGOraw.split("\t")[1])
				goea_go_pvalue = str(line_enGOraw.split("\t")[2])
				goea_go_adj_pvalue = str(line_enGOraw.split("\t")[3])
				goea_go_cats_test = int(line_enGOraw.split("\t")[8])
				goea_go_cats_ref = int(line_enGOraw.split("\t")[6])
				goea_go_total_test = int(line_enGOraw.split("\t")[7])
				goea_go_total_ref = int(line_enGOraw.split("\t")[5])
				goea_go_geneset = str("|".join([element.split(",")[-1].strip() for element in line_enGOraw.split("\t")[9].strip("[]").split(" - ")[:-1]]))
			else:
				continue
		elif str(goeatool).lower() == "agrigo":
#			if isint(line_enGOraw.split("\t")[0].split("GO:")[-1]):
			if line_enGOraw.startswith("GO:"):
				goea_go_id = str(line_enGOraw.split("\t")[0])
				goea_go_description = str(line_enGOraw.split("\t")[2])
				goea_go_pvalue = str(line_enGOraw.split("\t")[7])
				goea_go_adj_pvalue = str(line_enGOraw.split("\t")[8])
				goea_go_cats_test = int(line_enGOraw.split("\t")[3])
				goea_go_cats_ref = int(line_enGOraw.split("\t")[5])
				goea_go_total_test = int(line_enGOraw.split("\t")[4])
				goea_go_total_ref = int(line_enGOraw.split("\t")[6])
				goea_go_geneset = str("|".join(line_enGOraw.split("\t")[9].strip("// ").split(" // ")))
			else:
				continue
		elif str(goeatool).lower() == "goatools":
			if len(line_enGOraw.split("\t")) > 0 and line_enGOraw.split("\t")[0].strip(".").startswith("GO:"):
				goea_go_id = line_enGOraw.split("\t")[0].strip(".")
				goea_go_description = str(line_enGOraw.split("\t")[3])
				goea_go_pvalue = str(line_enGOraw.split("\t")[6])
				goea_go_adj_pvalue = str(line_enGOraw.split("\t")[9])
				goea_go_cats_test = int(line_enGOraw.split("\t")[4].split("/")[0])
				goea_go_cats_ref = int(line_enGOraw.split("\t")[5].split("/")[0])
				goea_go_total_test = int(line_enGOraw.split("\t")[4].split("/")[1])
				goea_go_total_ref = int(line_enGOraw.split("\t")[5].split("/")[1])
				goea_go_geneset = str("|".join(line_enGOraw.split("\t")[10].split(", ")))
			else:
				continue
		elif str(goeatool).lower() == "generic":
			if len(line_enGOraw.split("\t")) > 1 and line_enGOraw.startswith("GO:"):
				goea_go_id = str(line_enGOraw.split("\t")[0])
				goea_go_description = "na"
				goea_go_pvalue = 0
				goea_go_adj_pvalue = 0
				goea_go_cats_test = 0
				goea_go_cats_ref = 0
				goea_go_total_test = 0
				goea_go_total_ref = 0
				goea_go_geneset = str(line_enGOraw.split("\t")[1])
			else:
				continue
		else:
			print("The provided tool is currently not supported!") 
			break
			
		## Determine the level and depth of each GO term in the input
		if goea_go_id in obo_go_name_dict:
			if obo_go_is_obsolete_dict[goea_go_id] != "true":
				goea_go_type = "BP" if obo_go_namespace_dict[goea_go_id] == "biological_process" else "MF" if obo_go_namespace_dict[goea_go_id] == "molecular_function" else "CC" if obo_go_namespace_dict[goea_go_id] == "cellular_component" else ""
				if dswitch is True:
					if goea_go_id not in RootTerm_dict.values():
						goea_go_id_level = "D" + "%02d" % int(nx.shortest_path_length(go_hierarchy_network,source = RootTerm_dict[go_hierarchy_network.nodes[goea_go_id]["Namespace"]],target = goea_go_id) - 1)
						goea_go_id_depth = "D" + "%02d" % int(len(max(nx.all_simple_paths(go_hierarchy_network, source = RootTerm_dict[go_hierarchy_network.nodes[goea_go_id]["Namespace"]], target = goea_go_id), key = len)) - 1)
					else:
						goea_go_id_level = "L00"
						goea_go_id_depth = "D00"
					enGOfmtted_list.append(str(goea_go_id) + "\t" + str(goea_go_description) + "\t" + str(goea_go_type) + "\t" + str(goea_go_id_depth) + "\t" + str(goea_go_pvalue) + "\t" + str(goea_go_adj_pvalue) + "\t" + str(goea_go_cats_test) + "\t" + str(goea_go_cats_ref) + "\t" + str(goea_go_total_test) + "\t" + str(goea_go_total_ref) + "\t" + str(goea_go_geneset))
				else:
					enGOfmtted_list.append(str(goea_go_id) + "\t" + str(goea_go_description) + "\t" + str(goea_go_type) + "\t" + "na" + "\t" + str(goea_go_pvalue) + "\t" + str(goea_go_adj_pvalue) + "\t" + str(goea_go_cats_test) + "\t" + str(goea_go_cats_ref) + "\t" + str(goea_go_total_test) + "\t" + str(goea_go_total_ref) + "\t" + str(goea_go_geneset))
			else:
				print(str(goea_go_id) + " is labeled as \"obselete\" in the obo annotation file, will be skipped!")
		else:
			goea_ori_go_id_list = [key_go_id for key_go_id in obo_go_alt_id_dict if goea_go_id in obo_go_alt_id_dict[key_go_id]]
			if len(goea_ori_go_id_list) > 0:
				goea_ori_go_id = goea_ori_go_id_list[0]
				if obo_go_is_obsolete_dict[goea_ori_go_id] != "true":
					goea_go_type = "BP" if obo_go_namespace_dict[goea_ori_go_id] == "biological_process" else "MF" if obo_go_namespace_dict[goea_ori_go_id] == "molecular_function" else "CC" if obo_go_namespace_dict[goea_ori_go_id] == "cellular_component" else ""
					if dswitch is True:
						if goea_ori_go_id not in RootTerm_dict.values():
							goea_go_id_level = "D" + "%02d" % int(nx.shortest_path_length(go_hierarchy_network,source = RootTerm_dict[go_hierarchy_network.nodes[goea_ori_go_id]["Namespace"]],target = goea_ori_go_id) - 1)
							goea_go_id_depth = "D" + "%02d" % int(len(max(nx.all_simple_paths(go_hierarchy_network, source = RootTerm_dict[go_hierarchy_network.nodes[goea_ori_go_id]["Namespace"]], target = goea_ori_go_id), key = len)) - 1)
						else:
							goea_go_id_level = "L00"
							goea_go_id_depth = "D00"
						enGOfmtted_list.append(str(goea_go_id) + "\t" + str(goea_go_description) + "\t" + str(goea_go_type) + "\t" + str(goea_go_id_depth) + "\t" + str(goea_go_pvalue) + "\t" + str(goea_go_adj_pvalue) + "\t" + str(goea_go_cats_test) + "\t" + str(goea_go_cats_ref) + "\t" + str(goea_go_total_test) + "\t" + str(goea_go_total_ref) + "\t" + str(goea_go_geneset))
					else:
						enGOfmtted_list.append(str(goea_go_id) + "\t" + str(goea_go_description) + "\t" + str(goea_go_type) + "\t" + "na" + "\t" + str(goea_go_pvalue) + "\t" + str(goea_go_adj_pvalue) + "\t" + str(goea_go_cats_test) + "\t" + str(goea_go_cats_ref) + "\t" + str(goea_go_total_test) + "\t" + str(goea_go_total_ref) + "\t" + str(goea_go_geneset))
				else:
					print(str(goea_go_id) + " is labeled as \"obselete\" in the obo annotation file, will be skipped!")
			else:
				print(str(goea_go_id) + " is not found in the obo annotation file, will be skipped!")
				continue
	
	if dswitch is True:		
		enGOfmtted_list = sorted(enGOfmtted_list, key = lambda element: int(element.split("\t")[3].split("D")[1]))
	else:
		enGOfmtted_list = sorted(enGOfmtted_list, key = lambda element: int(element.split("\t")[6]), reverse = True)

	return enGOfmtted_list

def goea_filter(OBOInput, goeatool, enGOraw, gosize, gotype, dswitch = False):
	"""
	gene ontology enrichment analysis (goea) result formatter, supporting BiNGO, agriGO, AmiGO, PANTHER, GOrilla, gProfiler, Enrichr
	Required:
	OBOInput	obo file should be provided, e.g. go-basic.obo
	goeatool	the go enrichment tools used: [BiNGO, agriGO, AmiGO, PANTHER,GOrilla, gProfiler, Enrichr]
	enGOraw	the enrichment test reults from the go enrichment tool used
	gosize	threshold for the size of GO terms, only GO terms below this threshold will be printed out
	gotype  "BP","CC","MF" which category to use
	dswitch switch for depth ON and OFF
	"""
	enGOfmtted_list = goea_formatter(OBOInput, goeatool, enGOraw, dswitch)
	enGOfltrd_list = []
	for element_enGOfmtted in enGOfmtted_list:
		go_type = str(element_enGOfmtted.split("\t")[2])
		go_size_ref = int(element_enGOfmtted.split("\t")[7]) # number of genes in a go term from the reference annotation.
		if go_type in gotype and go_size_ref <= int(gosize):
			enGOfltrd_list.append(element_enGOfmtted)
	return enGOfltrd_list
		

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = synopsis, formatter_class = argparse.RawTextHelpFormatter)
	parser.add_argument("OBOInput", metavar = "-OBOInput", help = "OBO file from Gene Ontology", action = "store", nargs = None, const = None, default = None, type = None, choices = None) ## Cannot specify 'dest' for positional arguments
	parser.add_argument("enGO", metavar = "-enGO", help = "Enriched GO input file may be from different GO enrichment analysis tools (e.g. BiNGO, agriGO, AmiGO, etc..)", action = "store", nargs = None, const = None, default = None, type = None, choices = None) ## See below.
	parser.add_argument("-got", metavar = None, help = "GO enrichment tools used for enrichment test (default: %(default)s)", action = "store", nargs = None, const = None, default = "BiNGO", type = str, choices = ["BiNGO", "agriGO", "AmiGO", "PANTHER", "GOrilla", "gProfiler", "GOATOOLS"]) ## See below.
	parser.add_argument("-gosize", metavar = None, dest = None, help = "Threshold for the size of GO terms, only GO terms below this threshold will be printed out", action = "store", nargs = None, const = None, default = None, type = int, choices = None)
	parser.add_argument("-gotype", metavar = None, dest = None, help = "Type of GO terms, only GO terms in this or these categories will be printed out", action = "store", nargs = "*", const = None, default = "BP", type = None, choices = None)
	parser.add_argument("-d", dest = "dswitch", help = "Calculate depth for input GO terms", action = "store_true", default = None) ## Argument present --> true, not present --> false. The followings are not compatible with "store_true": metavar = None, nargs = None,const = None, type = None, choices = None
	args = parser.parse_args() 
#	print(args.OBOInput, args.got, args.enGO, args.gosize)

	if args.gosize is None:
		print("Printing full GO list")
		with open(os.path.splitext(args.enGO)[0] + ".enGOfmtted", "w") as fout_enGOfmtted:
			fout_enGOfmtted.write("Full GO-ID" + "\t" + "Description" + "\t" + "Type" + "\t" + "Depth" + "\t" + "p-value" + "\t" + "adj p-value" + "\t" + "x.cats.test" + "\t" + "n.cats.ref" + "\t" + "X.total.test" + "\t" + "N.total.ref" + "\t" + "Genes in test set" + "\n")
			enGOfmtted_list = goea_formatter(args.OBOInput, args.got, args.enGO, args.dswitch)
			for element_enGOfmtted in enGOfmtted_list:
				fout_enGOfmtted.write(element_enGOfmtted + "\n")
		fout_enGOfmtted.close()
	else:
		print("Printing filtered GO list")
		with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + ".enGOfltrd", "w") as fout_enGOfltrd:
			fout_enGOfltrd.write("Full GO-ID" + "\t" + "Description" + "\t" + "Type" + "\t" + "Depth" + "\t" + "p-value" + "\t" + "adj p-value" + "\t" + "x.cats.test" + "\t" + "n.cats.ref" + "\t" + "X.total.test" + "\t" + "N.total.ref" + "\t" + "Genes in test set" + "\n")
			enGOfltrd_list = goea_filter(args.OBOInput, args.got, args.enGO, args.gosize, args.gotype, args.dswitch)
			for element_enGOfltrd in enGOfltrd_list:
				fout_enGOfltrd.write(element_enGOfltrd + "\n")
		fout_enGOfltrd.close()
		

