#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__    = "Guannan Wang"
__email__     = "wangguannan2014@gmail.com"
__version__   = "1"

import argparse, sys, re, os, io, operator 
from scripts.go_enrichment_result_formatter import goea_formatter, goea_filter
from scripts.go_clustering import go_compare, go_assign_cluster
from scripts.go_plot import sim_plot, sim_newtork

synopsis = "\n\
#############################################################################################################################################\n\
#GOMCL.py clusters GO terms using MCL based on overlapping ratios, OC (Overlap coefficient) or JC (Jaccard coefficient).                    #\n\
#Use examples:                                                                                                                              #\n\
#   GOMCL.py -Ct 0.5 -I 1.5 EnrichedGO.txt                                                                                                  #\n\
#   GOMCL.py -SI JC -Ct 0.5 -I 1.5 EnrichedGO.txt                                                                                           #\n\
#                                                                                                                                           #\n\
#############################################################################################################################################"

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = synopsis, formatter_class = argparse.RawTextHelpFormatter)
	parser.add_argument('OBOInput', metavar = '-OBO', help = 'obo file should be provided, e.g. go-basic.obo', action = 'store', nargs = None, const = None, default = None, type = None, choices = None)
	parser.add_argument("enGO", metavar = "-enGO", help = "Enriched GO input file may be from different GO enrichment analysis tools (e.g. BiNGO, agriGO, AmiGO, etc..", action = "store", nargs = None, const = None, default = None, type = None, choices = None)
	parser.add_argument("-got", metavar = None, help = "GO enrichment tools used for enrichment test (default: %(default)s)", action = "store", nargs = None, const = None, default = "BiNGO", type = str, choices = ["BiNGO", "agriGO", "AmiGO", "PANTHER", "GOrilla", "gProfiler", "Enrichr"])
	parser.add_argument("-gosize", metavar = None, dest = None, help = "Threshold for the size of GO terms, only GO terms below this threshold will be printed out (default: %(default)s)", action = "store", nargs = None, const = None, default = "3000", type = int, choices = None)
	parser.add_argument("-SI", metavar = None, dest = "simindex", help = "Method to calculate similarity between GO terms, OC (Overlap coefficient) or JC (Jaccard coefficient) (default: %(default)s)", action = "store", nargs = None, const = None, default = "OC", type = str, choices = ["OC", "JC"]) 
	parser.add_argument("-Ct", metavar = None, dest = "cutoff", help = "Clustering threshold for the overlapping ratio between two GO terms, any value between 0 and 1 (default: %(default)s)", action = "store", nargs = None, const = None, default = "0.5", type = float, choices = None) 
	parser.add_argument("-I", metavar = None, dest = "inflation", help = "Inflation value, main handle for cluster granularity, usually chosen somewhere in the range [1.2-5.0] (default: %(default)s)", action = "store", nargs = None, const = None, default = "2.0", type = float, choices = None) 
	parser.add_argument("-SL", metavar = None, dest = "sig", help = "Signifance level (p-value cutoff) used in the enrichment test (default: %(default)s)", action = "store", nargs = None, const = None, default = "0.05", type = float, choices = None)
	parser.add_argument("-hm",  dest = None, help = "Only needed if a similarity heatmap is desired", action = "store_true", default = None ) 
	parser.add_argument("-nx",  dest = None, help = "Only needed if a similarity-based network is desired", action = "store_true", default = None ) 
	args = parser.parse_args() 
	
	enGOfltrd_list = goea_filter(args.OBOInput, args.got, args.enGO, args.gosize)
	with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + ".enGOfltrd", "w") as fin_enGOfltrd:
		fin_enGOfltrd.write("Full GO-ID" + "\t" + "Description" + "\t" + "Type" + "\t" + "Depth" + "\t" + "p-value" + "\t" + "adj p-value" + "\t" + "x.cats.test" + "\t" + "n.cats.ref" + "\t" + "X.total.test" + "\t" + "N.total.ref" + "\t" + "Genes in test set" + "\n")
		for element_enGOfltrd in enGOfltrd_list:
			fin_enGOfltrd.write(element_enGOfltrd + "\n")
	fin_enGOfltrd.close()

	enGOfmtfltr_info_dict, gosim_dict = go_compare(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + ".enGOfltrd", args.simindex)
	
	with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_Ct" + str(args.cutoff) + ".simfltred", "w") as fout_simfltr:
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

	
	go_clstr_dict, gene_clstr_dict, clstred_go_dict, clstred_gene_dict = go_assign_cluster(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + ".enGOfltrd", args.simindex, args.cutoff, args.inflation)
	with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".clstrinfo", "w") as fout_clstrinfo:
		fout_clstrinfo.write("GO.Clstr" + "\t" + "# of GOs" + "\t" + "# of genes" + "\n")
		for clstrid in clstred_go_dict:
			fout_clstrinfo.write("Ct" + str(args.cutoff) + "I" + str(args.inflation) + "_" + str(clstrid) + "\t" + str(len(set(clstred_go_dict[clstrid]))) + "\t" + str(len(set(clstred_gene_dict[clstrid]))) + "\n")
	fout_clstrinfo.close()
	
	with open(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + "_Ct" + str(args.cutoff) + "I" + str(args.inflation) + ".clstr", "w") as fout_clstr:
		fout_clstr.write("Full GO-ID" + "\t" + "Description" + "\t" + "Type" + "\t" + "Depth" + "\t" + "p-value" + "\t" + "adj p-value" + "\t" + "x.cats.test" + "\t" + "n.cats.ref" + "\t" + "X.total.test" + "\t" + "N.total.ref" + "\t" + "Genes in test set" + "\t" + "Clstr" + "\n")
		try:
			for key_enGO in sorted(enGOfmtfltr_info_dict, key = lambda dict_key : (int(go_clstr_dict[dict_key]), int(enGOfmtfltr_info_dict[dict_key].split("\t")[3].split("D")[1]), -int(enGOfmtfltr_info_dict[dict_key].split("\t")[7]))):
				fout_clstr.write("\t".join(map(str,enGOfmtfltr_info_dict[key_enGO].split("\t"))) + "\t" + str(go_clstr_dict[key_enGO]) + "\n")
		except TypeError:
			for key_enGO in sorted(enGOfmtfltr_info_dict, key = lambda dict_key : (int(go_clstr_dict[dict_key]),int(enGOfmtfltr_info_dict[dict_key].split("\t")[3].split("D")[1]))):
				fout_clstr.write("\t".join(map(str,enGOfmtfltr_info_dict[key_enGO].split("\t"))) + "\t" + str(go_clstr_dict[key_enGO]) + "\n")
	fout_clstr.close()
	
	if args.hm:
		sim_plot(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + ".enGOfltrd", args.simindex, args.cutoff, args.inflation)
	if args.nx:
		sim_newtork(os.path.splitext(args.enGO)[0] + "_GOsize" + str(args.gosize) + ".enGOfltrd", args.simindex, args.cutoff, args.inflation, args.sig)

	
	
	
	
	
	
	
	
	