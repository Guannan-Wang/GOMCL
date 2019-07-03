# GOMCL
a Python tool for Gene Ontology gene sets clustering

## Getting Started
GOMCL.py clusters GO terms using MCL based on overlapping ratios, OC (Overlap coefficient) or JC (Jaccard coefficient). 

Use examples:

GOMCL.py -Ct 0.5 -I 1.5 OBOfile EnrichedGO

GOMCL.py -SI JC -Ct 0.5 -I 1.5 OBOfile EnrichedGO

GOMCL.py -Ct 0.5 -I 1.5 -nw -hm OBOfile EnrichedGO

Note:

OBO file can be downloaded from http://purl.obolibrary.org/obo/go.obo or http://purl.obolibrary.org/obo/go/go-basic.obo (more explanations for obo file can be found here: http://geneontology.org/docs/download-ontology/)

## Prerequisites
* python 2.7

* python libraries: networkx, matplotlib, seaborn, numpy, pandas

* [MCL](https://micans.org/mcl/index.html) (See installation instructions here: https://micans.org/mcl/)

## Installation
**1. Download zip file and install**
```
 wget https://github.com/Guannan-Wang/GOMCL/archive/master.zip
 unzip master.zip
 cd GOMCL-master
 chmod 755 *.py scripts/*.py
 export PATH=/path/to/GOMCL-master:$PATH 
```

## Running the test

Will be available soon
