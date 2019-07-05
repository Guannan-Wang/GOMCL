# GOMCL
a Python tool for Gene Ontology gene sets clustering

## Getting Started
GOMCL.py clusters GO terms using MCL based on overlapping ratios, OC (*Overlap coefficient*) or JC (*Jaccard coefficient*). 

#### Use examples:

GOMCL.py -Ct 0.5 -I 1.5 OBOfile EnrichedGO

GOMCL.py -SI JC -Ct 0.5 -I 1.5 OBOfile EnrichedGO

GOMCL.py -Ct 0.5 -I 1.5 -nw -hm OBOfile EnrichedGO

#### Note:

1. OBO file can be downloaded from http://purl.obolibrary.org/obo/go.obo or http://purl.obolibrary.org/obo/go/go-basic.obo (more explanations for obo files can be found here: http://geneontology.org/docs/download-ontology/)

2. GOMCL is currently compatible with BiNGO,agriGO,GOrilla,g:Profiler and customized inputs. Support for other tools will be added upon request. 

3. Similarity between GO terms is calculateed either as *Jaccard Coefficient* (JC) or *Overlap Coefficient* (OC), as described in [Merico et al., 2010](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0013984).

4. The use of -Ct and -I values heavily depends on the number of input GO terms and how similar they are. 

See -h or --help for more options

## Prerequisites
* python 2.7

* python libraries: networkx, matplotlib, seaborn, numpy, pandas

* [MCL](https://micans.org/mcl/) (It can be downloaded from https://www.micans.org/mcl/src/ or https://github.com/JohannesBuchner/mcl)

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
