#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__    = "Guannan Wang"
__email__     = "wangguannan2014@gmail.com"
__version__   = "1"

import argparse, sys, re, os, io, operator 
from scripts.go_enrichment_result_formatter import goea_formatter, goea_filter
from scripts.go_clustering import go_compare, go_assign_cluster
from scripts.go_plot import sim_plot, sim_density, sim_newtork, construct_go_hierarchy_subgraph, sim_cumulative

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

synopsis = "\n\
#############################################################################################################################################\n\
#GOMCL.py clusters GO terms using MCL based on overlapping ratios, OC (Overlap coefficient) or JC (Jaccard coefficient).                    #\n\
#Use examples:                                                                                                                              #\n\
#   GOMCL.py -Ct 0.5 -I 1.5 OBOfile EnrichedGO.txt                                                                                          #\n\
#   GOMCL.py -SI JC -Ct 0.5 -I 1.5 OBOfile EnrichedGO.txt                                                                                   #\n\
#                                                                                                                                           #\n\
#############################################################################################################################################"

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = synopsis, formatter_class = argparse.RawTextHelpFormatter)
	parser.add_argument('OBOInput', metavar = '-OBO', help = 'obo file should be provided, e.g. go-basic.obo', action = 'store', nargs = None, const = None, default = None, type = None, choices = None)
	parser.add_argument("enGO", metavar = "-enGO", help = "Enriched GO input file may be from different GO enrichment analysis tools (e.g. BiNGO, agriGO, AmiGO, etc..)", action = "store", nargs = None, const = None, default = None, type = None, choices = None)
	parser.add_argument("-d", dest = "dswitch", help = "Only needed if depth for input GO terms is desired", action = "store_true", default = None)
	parser.add_argument("-got", metavar = None, help = "GO enrichment tools used for enrichment test (default: %(default)s)", action = "store", nargs = None, const = None, default = "BiNGO", type = str, choices = ["BiNGO", "agriGO", "AmiGO", "PANTHER", "GOrilla", "gProfiler", "GOATOOLS", "generic"])
	parser.add_argument("-gosize", metavar = None, dest = None, help = "Threshold for the size of GO terms, only GO terms below this threshold will be printed out (default: %(default)s)", action = "store", nargs = None, const = None, default = 3000, type = int, choices = None)
	parser.add_argument("-gotype", metavar = None, dest = None, help = "Type of GO terms, only GO terms in this or these categories will be printed out", action = "store", nargs = "*", const = None, default = "BP", type = None, choices = None)
	parser.add_argument("-SI", metavar = None, dest = "simindex", help = "Method to calculate similarity between GO terms, OC (Overlap coefficient) or JC (Jaccard coefficient) (default: %(default)s)", action = "store", nargs = None, const = None, default = "OC", type = str, choices = ["OC", "JC"]) 
	parser.add_argument("-Ct", metavar = None, dest = "cutoff", help = "Clustering threshold for the overlapping ratio between two GO terms, any value between 0 and 1 (default: %(default)s)", action = "store", nargs = None, const = None, default = 0.5, type = float, choices = None) 
	parser.add_argument("-I", metavar = None, dest = "inflation", help = "Inflation value, main handle for cluster granularity, usually chosen somewhere in the range [1.2-5.0] (default: %(default)s)", action = "store", nargs = None, const = None, default = 1.5, type = float, choices = None) 
	parser.add_argument("-Sig", metavar = None, dest = "sig", help = "Signifance level (p-value cutoff) used in the enrichment test (default: %(default)s)", action = "store", nargs = None, const = None, default = 0.05, type = float, choices = None)
	parser.add_argument("-ssd", metavar = None, dest = None, help = "Only needed if a similarity score distribution is desired for clusters with number of GOs larger than this threshold", action = "store", nargs = "?", const = 0, default = None, type = int, choices = None)
	parser.add_argument("-hg", metavar = None, dest = None, help = "Only needed if a hierarchy graph is desired for clusters with number of GOs larger than this threshold", action = "store", nargs = "?", const = 0, default = None, type = int, choices = None)
	parser.add_argument("-hgt", dest = None, help = "Only needed if a tabular output of the GO hierarchy is desired for the clusters specified by option -hg, should always be used with option -hg", action = "store_true", default = None)
	parser.add_argument("-hm",  dest = None, help = "Only needed if a similarity heatmap is desired", action = "store_true", default = None ) 
	parser.add_argument("-nw",  dest = None, help = "Only needed if a similarity-based network is desired", action = "store_true", default = None ) 
	args = parser.parse_args() 
	
	enGOfltrd_list = goea_filter(args.OBOInput, args.got, args.enGO, args.gosize, args.gotype, args.dswitch)
	with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".enGOfltrd.temp", "w") as fin_enGOfltrd_temp:
		fin_enGOfltrd_temp.write("Full GO-ID" + "\t" + "Description" + "\t" + "Type" + "\t" + "Depth" + "\t" + "p-value" + "\t" + "adj p-value" + "\t" + "x.cats.test" + "\t" + "n.cats.ref" + "\t" + "X.total.test" + "\t" + "N.total.ref" + "\t" + "Genes in test set" + "\n")
		for element_enGOfltrd in enGOfltrd_list:
			fin_enGOfltrd_temp.write(element_enGOfltrd + "\n")
	fin_enGOfltrd_temp.close()

	enGOfmtfltr_info_dict, gosim_dict = go_compare(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".enGOfltrd.temp", 1, 0, args.simindex)
	
	with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".simfltred", "w") as fout_simfltr:
		fout_simfltr.write("GOtermID-A"	+ "\t" + "GOtermID-B" + "\t" + "Similarity (" + str(args.simindex) + ")" + "\n")
		processed = []
		for key_querygo in gosim_dict:
			for key_subjectgo in gosim_dict:
				if str(key_querygo) <> str(key_subjectgo):
					pair_qs = str(key_querygo) + "_" + str(key_subjectgo)
					pair_sq = str(key_subjectgo) + "_" + str(key_querygo)
					if pair_qs not in processed and pair_sq not in processed:
						if max(gosim_dict[key_querygo][key_subjectgo], gosim_dict[key_subjectgo][key_querygo]) >= float(args.cutoff):
							fout_simfltr.write(str(key_querygo) + "\t" + str(key_subjectgo) + "\t" +  str(max(gosim_dict[key_querygo][key_subjectgo], gosim_dict[key_subjectgo][key_querygo])) + "\n")
							processed.extend([pair_qs,pair_sq])
				else:
					fout_simfltr.write(str(key_querygo) + "\t" + "" + "\t" + "" + "\n")

	go_clstr_dict, gene_clstr_dict, clstred_go_dict, clstred_gene_dict = go_assign_cluster(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".enGOfltrd.temp", args.simindex, args.cutoff, args.inflation)
	
	with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".clstr", "w") as fout_clstr:
		fout_clstr.write("Clstr" + "\t" + "Full GO-ID" + "\t" + "Description" + "\t" + "Type" + "\t" + "Depth" + "\t" + "p-value" + "\t" + "adj p-value" + "\t" + "x.cats.test" + "\t" + "n.cats.ref" + "\t" + "X.total.test" + "\t" + "N.total.ref" + "\t" + "Genes in test set" + "\n")
		if args.dswitch:
			try:
				for key_enGO in sorted(enGOfmtfltr_info_dict, key = lambda dict_key : (int(go_clstr_dict[dict_key]), int(enGOfmtfltr_info_dict[dict_key].split("\t")[3].split("D")[1]), -int(enGOfmtfltr_info_dict[dict_key].split("\t")[7]))):
					fout_clstr.write(str(go_clstr_dict[key_enGO]) + "\t" + "\t".join(map(str,enGOfmtfltr_info_dict[key_enGO].split("\t"))) + "\n")
			except TypeError:
				for key_enGO in sorted(enGOfmtfltr_info_dict, key = lambda dict_key : (int(go_clstr_dict[dict_key]), int(enGOfmtfltr_info_dict[dict_key].split("\t")[3].split("D")[1]))):
					fout_clstr.write(str(go_clstr_dict[key_enGO]) + "\t" + "\t".join(map(str,enGOfmtfltr_info_dict[key_enGO].split("\t"))) + "\n")
		else:
			for key_enGO in sorted(enGOfmtfltr_info_dict, key = lambda dict_key : (int(go_clstr_dict[dict_key]), -int(enGOfmtfltr_info_dict[dict_key].split("\t")[7]))):
					fout_clstr.write(str(go_clstr_dict[key_enGO]) + "\t" + "\t".join(map(str,enGOfmtfltr_info_dict[key_enGO].split("\t"))) + "\n")
	fout_clstr.close()
	
#	with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".clstrinfo", "w") as fout_clstrinfo:
#		fout_clstrinfo.write("GO.Clstr" + "\t" + "# of GOs" + "\t" + "# of genes" + "\n")
#		for clstrid in clstred_go_dict:
#			fout_clstrinfo.write(str(args.simindex) + "Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_" + str(clstrid) + "\t" + str(len(set(clstred_go_dict[clstrid]))) + "\t" + str(len(set(clstred_gene_dict[clstrid]))) + "\n")
#	fout_clstrinfo.close()	

	with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".clstrinfo", "w") as fout_clstrinfo:
		Accumgenelist = []
		if args.nw:
			go_sim_newtork = sim_newtork(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".clstr", args.simindex, args.cutoff, args.inflation, args.sig)
			fout_clstrinfo.write("GO.Clstr" + "\t" + "# of GOs" + "\t" + "# of genes" + "\t" + "Accum # of genes" + "\t" + "Largest size (Description) (x/n, p-value)" + "\t" + "Smallest p-value (Description) (x/n, p-value)" + "\t" + "Most connected (Description) (x/n, p-value)" + "\n")
			for clstrid in clstred_go_dict:
				pvaluesortedlist = sorted(clstred_go_dict[clstrid], key = lambda element_GO : float(enGOfmtfltr_info_dict[element_GO].split("\t")[5]))
				sizesortedlist = sorted(clstred_go_dict[clstrid], key = lambda element_GO : int(enGOfmtfltr_info_dict[element_GO].split("\t")[7]), reverse = True)
				#minpvalueGO = min(clstred_go_dict[clstrid], key = lambda element_GO : float(enGOfmtfltr_info_dict[element_GO].split("\t")[5]))
				#largestGO = max(clstred_go_dict[clstrid], key = lambda element_GO : int(enGOfmtfltr_info_dict[element_GO].split("\t")[7]))
				edgesortedlist = sorted(clstred_go_dict[clstrid], key = lambda element_GO : int(dict(go_sim_newtork.degree())[element_GO]), reverse = True)
				Accumgenelist.extend(clstred_gene_dict[clstrid])
				
				fout_clstrinfo.write(str(args.simindex) + "Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_" + str(clstrid) + "\t" + str(len(set(clstred_go_dict[clstrid]))) + "\t" + str(len(set(clstred_gene_dict[clstrid]))) + "\t" + str(len(set(Accumgenelist))) + "\t" + str(sizesortedlist[0]) + " (" + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[1]) + ")" + " (" + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[6]) + "/" + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[7]) + ", " + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[5]) + ")" + "\t" + str(pvaluesortedlist[0]) + " (" + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[1]) + ")" + " (" + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[6]) + "/" + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[7]) + ", " + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[5]) + ")" + "\t" + str(edgesortedlist[0]) + " (" + str(enGOfmtfltr_info_dict[edgesortedlist[0]].split("\t")[1]) + ")" + " (" + str(enGOfmtfltr_info_dict[edgesortedlist[0]].split("\t")[6]) + "/" + str(enGOfmtfltr_info_dict[edgesortedlist[0]].split("\t")[7]) + ", " + str(enGOfmtfltr_info_dict[edgesortedlist[0]].split("\t")[5]) + ")"	+ "\n")
#				fout_clstrinfo.write(str(args.simindex) + "Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_" + str(clstrid) + "\t" + str(len(set(clstred_go_dict[clstrid]))) + "\t" + str(len(set(clstred_gene_dict[clstrid]))) + "\t" + str(len(set(Accumgenelist))) + "\t" + str(sizesortedlist[0]) + " (" + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[1]) + ")" + " (" + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[6]) + "/" + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[7]) + ", " + "{:.0E}".format(float(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[5])) + ")" + "\t" + str(pvaluesortedlist[0]) + " (" + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[1]) + ")" + " (" + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[6]) + "/" + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[7]) + ", " + "{:.0E}".format(float(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[5])) + ")" + "\t" + str(edgesortedlist[0]) + " (" + str(enGOfmtfltr_info_dict[edgesortedlist[0]].split("\t")[1]) + ")" + " (" + str(enGOfmtfltr_info_dict[edgesortedlist[0]].split("\t")[6]) + "/" + str(enGOfmtfltr_info_dict[edgesortedlist[0]].split("\t")[7]) + ", " + "{:.0E}".format(float(enGOfmtfltr_info_dict[edgesortedlist[0]].split("\t")[5])) + ")"	+ "\n")				
		else:
			fout_clstrinfo.write("GO.Clstr" + "\t" + "# of GOs" + "\t" + "# of genes" + "\t" + "Accum # of genes" + "\t" + "Largest size (Description) (x/n, p-value)" + "\t" + "Smallest p-value (Description) (x/n, p-value)" + "\n")
			for clstrid in clstred_go_dict:
				pvaluesortedlist = sorted(clstred_go_dict[clstrid], key = lambda element_GO : float(enGOfmtfltr_info_dict[element_GO].split("\t")[5]))
				sizesortedlist = sorted(clstred_go_dict[clstrid], key = lambda element_GO : int(enGOfmtfltr_info_dict[element_GO].split("\t")[7]), reverse = True)
				Accumgenelist.extend(clstred_gene_dict[clstrid])
				
				fout_clstrinfo.write(str(args.simindex) + "Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_" + str(clstrid) + "\t" + str(len(set(clstred_go_dict[clstrid]))) + "\t" + str(len(set(clstred_gene_dict[clstrid]))) + "\t" + str(len(set(Accumgenelist))) + "\t" + str(sizesortedlist[0]) + " (" + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[1]) + ")" + " (" + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[6]) + "/" + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[7]) + ", " + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[5]) + ")" + "\t" + str(pvaluesortedlist[0]) + " (" + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[1]) + ")" + " (" + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[6]) + "/" + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[7]) + ", " + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[5]) + ")" + "\n")
#				fout_clstrinfo.write(str(args.simindex) + "Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_" + str(clstrid) + "\t" + str(len(set(clstred_go_dict[clstrid]))) + "\t" + str(len(set(clstred_gene_dict[clstrid]))) + "\t" + str(len(set(Accumgenelist))) + "\t" + str(sizesortedlist[0]) + " (" + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[1]) + ")" + " (" + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[6]) + "/" + str(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[7]) + ", " + "{:.0E}".format(float(enGOfmtfltr_info_dict[sizesortedlist[0]].split("\t")[5])) + ")" + "\t" + str(pvaluesortedlist[0]) + " (" + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[1]) + ")" + " (" + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[6]) + "/" + str(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[7]) + ", " + "{:.0E}".format(float(enGOfmtfltr_info_dict[pvaluesortedlist[0]].split("\t")[5])) + ")" + "\n")
	fout_clstrinfo.close()
	
	
	if args.hm:
		sim_plot(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".clstr", args.simindex, args.cutoff, args.inflation)
	"""	
	if args.nw:
		sim_newtork(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".clstr", args.simindex, args.cutoff, args.inflation, args.sig)
	"""
	
		
	if args.ssd is not None:
		for clstrnum in clstred_go_dict:
			if len(clstred_go_dict[clstrnum]) >= max(int(args.ssd), 2):
				simlist = []
				processed = []
				for clstredGO_A in clstred_go_dict[clstrnum]:
					for clstredGO_B in clstred_go_dict[clstrnum]:
						if str(clstredGO_A) <> str(clstredGO_B):
							pairAB = str(clstredGO_A) + str(clstredGO_B)
							pairBA = str(clstredGO_B) + str(clstredGO_A)
							if pairAB not in processed and pairBA not in processed:
								simlist.append(gosim_dict[clstredGO_A][clstredGO_B])
								processed.extend([pairAB, pairBA])
				try:
					print("GOsize{0}_{1}_Ct{2}_I{3}_C{4}_P(>=0.5): {5:.2%}".format(args.gosize, args.simindex, args.cutoff, args.inflation, clstrnum, len([sim for sim in simlist if float(sim) >= float(args.cutoff)])/float(len(simlist))))
				except ZeroDivisionError:
					print("GOsize{0}_{1}_Ct{2}_I{3}_C{4}_P(>=0.5): Only one GO term is present in this cluster.".format(args.gosize, args.simindex, args.cutoff, args.inflation, clstrnum))
				sim_density_plot = sim_density(simlist, args.simindex, increment = 0.02)
				plt.savefig(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_C" + str(clstrnum) + "_simden.png", dpi = 600, transparent = False, bbox_inches = 'tight', pad_inches = 0)
				plt.close()
				
				sim_cumulative_plot = sim_cumulative(simlist, args.simindex, args.cutoff, increment = 0.02)
				plt.savefig(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_C" + str(clstrnum) + "_simcum.png", dpi = 600, transparent = False, bbox_inches = 'tight', pad_inches = 0)
				plt.close()
	if args.hg is not None:	
		for clstrnum in clstred_go_dict:
			if len(clstred_go_dict[clstrnum]) >= int(args.hg):
				bpgo_list = []
				mfgo_list = []
				ccgo_list = []
				for clstred_go in clstred_go_dict[clstrnum]:
					if enGOfmtfltr_info_dict[clstred_go].split("\t")[2] == "BP":
						bpgo_list.append(clstred_go)
					if enGOfmtfltr_info_dict[clstred_go].split("\t")[2] == "MF":
						mfgo_list.append(clstred_go)
					if enGOfmtfltr_info_dict[clstred_go].split("\t")[2] == "CC":
						ccgo_list.append(clstred_go)
				if bpgo_list:
					clstredgo_hierarchy_network, clstredgo_hierarchy_subgraph = construct_go_hierarchy_subgraph(args.OBOInput, bpgo_list, sig = args.sig, GOinfodict = enGOfmtfltr_info_dict)
					plt.savefig(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_C" + str(clstrnum) + "_BPhierplot.png", dpi = 600, transparent = False, bbox_inches = 'tight', pad_inches = 0)
					plt.close()
					if args.hgt:
						with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_C" + str(clstrnum) + ".BPhgt", "w") as fin_bpgo_hgt:
							fin_bpgo_hgt.write("GOtermA" + "\t" + "GOtermB" + "\t" + "weight" + "\n")
							for edge in clstredgo_hierarchy_network.edges(data = True):
								nodeA, nodeB, attributes = edge
								fin_bpgo_hgt.write("{0}\t{1}\t{2}\n".format(nodeA, nodeB, attributes["weight"]))
						fin_bpgo_hgt.close()
				if mfgo_list:
					clstredgo_hierarchy_network, clstredgo_hierarchy_subgraph = construct_go_hierarchy_subgraph(args.OBOInput, mfgo_list, sig = args.sig, GOinfodict = enGOfmtfltr_info_dict)
					plt.savefig(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_C" + str(clstrnum) + "_MFhierplot.png", dpi = 600, transparent = False, bbox_inches = 'tight', pad_inches = 0)
					plt.close()
					if args.hgt:
						with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_C" + str(clstrnum) + ".MFhgt", "w") as fin_mfgo_hgt:
							fin_mfgo_hgt.write("GOtermA" + "\t" + "GOtermB" + "\t" + "weight" + "\n")
							for edge in clstredgo_hierarchy_network.edges(data = True):
								nodeA, nodeB, attributes = edge
								fin_mfgo_hgt.write("{0}\t{1}\t{2}\n".format(nodeA, nodeB, attributes["weight"]))
						fin_mfgo_hgt.close()
				if ccgo_list:
					clstredgo_hierarchy_network, clstredgo_hierarchy_subgraph = construct_go_hierarchy_subgraph(args.OBOInput, ccgo_list, sig = args.sig, GOinfodict = enGOfmtfltr_info_dict)
					plt.savefig(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_C" + str(clstrnum) + "_CChierplot.png", dpi = 600, transparent = False, bbox_inches = 'tight', pad_inches = 0)
					plt.close()
					if args.hgt:
						with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_C" + str(clstrnum) + ".CChgt", "w") as fin_ccgo_hgt:
							fin_ccgo_hgt.write("GOtermA" + "\t" + "GOtermB" + "\t" + "weight" + "\n")
							for edge in clstredgo_hierarchy_network.edges(data = True):
								nodeA, nodeB, attributes = edge
								fin_ccgo_hgt.write("{0}\t{1}\t{2}\n".format(nodeA, nodeB, attributes["weight"]))
						fin_ccgo_hgt.close()
	os.remove(str(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_" + str(args.simindex) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".enGOfltrd.temp"))
	
	
	
	
	
	
	
	