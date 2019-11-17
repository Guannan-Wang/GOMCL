# GOMCL
## Overview
GOMCL is a tool to cluster and extract summarized associations of Gene Ontology based functions in omics data. It clusters GO terms using MCL based on overlapping ratios, OC (*Overlap coefficient*) or JC (*Jaccard coefficient*). The resulting clusters can be further analyzed and separated into sub-clusters using a second script, GOMCL-sub. This tool helps researchers to reduce time spent on manual curation of large lists of GO terms and minimize biases introduced by shared GO terms in data interpretation. 

---
## Getting Started
## Prerequisites - programs
* python 2.7
* python libraries: 
  - [networkx](https://networkx.github.io/)
  - [matplotlib](https://matplotlib.org/)
  - [seaborn](https://seaborn.pydata.org/)
  - [numpy](numpy)
  - [pandas](https://pandas.pydata.org/)
* [MCL](https://micans.org/mcl/) - other clustering methods will be added in the future
  - MCL can also be downloaded from https://www.micans.org/mcl/src/ or https://github.com/JohannesBuchner/mcl
  - MCL should be added to the running environment. The user can use "mcl -h" to check whether MCL is loaded into the running environment. 
## Prerequisites - input data
* An obo annotation file is required, and can be downloaded downloaded from http://purl.obolibrary.org/obo/go.obo or http://purl.obolibrary.org/obo/go/go-basic.obo 
  - More explanations for obo files can be found here: http://geneontology.org/docs/download-ontology/
## Installation
**1. Download zip file and install**
```
 wget https://github.com/Guannan-Wang/GOMCL/archive/master.zip
 unzip master.zip
 cd GOMCL-master
 chmod 755 *.py scripts/*.py
 export PATH=/path/to/GOMCL-master:$PATH 
```
After installation, please check if GOMCL is properly installed by simply typing the following command:
```
GOMCL.py -h # This will print all options for GOMCL
```
## Ready to run GOMCL
### Use examples:
#### GOMCL:
```
GOMCL.py OBOfile EnrichedGO -d -gosize 3500 -Ct 0.5 -I 1.5 
GOMCL.py OBOfile EnrichedGO -Ct 0.5 -I 1.5 -hm -nw -d -hg 0 -hgt -ssd 0 
```
#### GOMCL-sub:
```
GOMCL-sub.py OBOfile ClstrGO -C 1 -gosize 2000 -I 1.8 -ssd 0 -hg 0 -hgt -hm -nw # Cluster C1 will be further separated.
```
### Option expalanations:
This can be accessed by -h or --help.
#### required arguments:
```
  -OBO		obo file should be provided, e.g. go-basic.obo
  -enGO		Enriched GO input file may be from different GO enrichment analysis tools, currently supported GO enrichment tools are: BiNGO, agriGO, GOrilla, gProfiler
```
#### optional arguments:
```
  -d		Only needed if depth for input GO terms is desired 
  -got		GO enrichment tools used for enrichment test (default: BiNGO), 
  -gosize	Threshold for the size of GO terms, only GO terms below this threshold will be printed out (default: 3000)
  -gotype	Type of GO terms, only GO terms in this or these categories will be printed out 
  -SI		Method to calculate similarity between GO terms, OC (Overlap coefficient) or JC (Jaccard coefficient) (default: OC)
  -Ct		Clustering threshold for the overlapping ratio between two GO terms, any value between 0 and 1 (default: 0.5)
  -I		Inflation value, main handle for cluster granularity, usually chosen somewhere in the range [1.2-5.0] (default: 1.5)
  -Sig		Signifance level (p-value cutoff) used in the enrichment test (default: 0.05)
  -ssd		Only needed if a similarity score distribution is desired for clusters with number of GOs larger than this threshold
  -hg		Only needed if a hierarchy graph is desired for clusters with number of GOs larger than this threshold
  -hgt		Only needed if a tabular output of the GO hierarchy is desired for the clusters specified by option -hg, should always be used with option -hg
  -hm		Only needed if a similarity heatmap is desired
  -nw		Only needed if a similarity-based network is desired
```
### Note:
1. GOMCL is currently compatible with BiNGO, agriGO, GOrilla, g:Profiler and customized inputs. Support for other tools will be added upon request. 
2. Similarity between GO terms is calculated either as *Jaccard Coefficient* (JC) or *Overlap Coefficient* (OC), as described in [Merico et al., 2010](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0013984).
3. The use of -Ct and -I values heavily depends on the number of input GO terms and how similar they are. It is suggested to try different -Ct and -I values to select the best combination.
## Running the test
Will be available soon
