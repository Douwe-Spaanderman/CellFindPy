# CellFindPy

CellFindPy is a program that uses Scanpy (https://github.com/theislab/scanpy) to cluster and further subcluster cells by single cell RNA sequencing 10x Data in an automated
fashion.

# What to expect

CellFindPy analysis 10x Data in the following steps.

- First it initializes an Anndata (https://github.com/theislab/anndata) object and filters out cells and genes as described in de Scanpy api. For example, you can set to 
filter out cells that have x% their gene expression assigned to mitochondrial genes.

- Secondly, it calculates pca, neighbor graph and find clusters using Scanpy (louvain clustering). 

- After clustering, marker genes are identified using both average differentiated expression and Mann–Whitney U test (to calculate p-values). 

- With this information a cluster test is conducted to identify if the set resolution gives distinct enough clusters. The thresholds for the cluster test is fully custumizable.
As an  example, here we used a threshold of at least 10 genes that are >log(2) differentially expressed between clusters.

- Lastly, looping through this system is done to find the highest possible resolution where the thresholds are met. Further subclustering and even subsubclustering can be done
using this same 5 step procedure 

Importantly, filters and threshold are fully custumizable in our scripts and can be interactively setup at the start of initializing using the _setup_.sh.

Alternatively, clustering can be initialized using the run_script.py and seperate functions can also be used by import Lib_.py. 

# Benchmark

To be determined

# How to use

### Install

First step of install is only required to be able to use both _setup_.sh and run_script.py. If you want to use functions separately than it is advised to use pip.
```
git clone https://github.com/Douwe-Spaanderman/CellFindPy
cd CellFindPy/
pip install .
```

### Run
Initializing using _setup_.sh should install all required packages used in CellFindPy. 

To run open terminal:
```
cd /CellFindPy/CellFindPy #locate _setup_.sh
bash _setup_.sh
```

Follow the instructions given by _setup_.sh to initialize CellFindPy.
Alternatively:
```
cd /CellFindPy/CellFindPy #locate run_script.py
python3 run_script.py -h
```

-h will show you all the possible arguments to set. Only the input folder (with the matrix.mtx, genes.tsv and barcodes.tsv), output folder and save name are required.

# License

To be determined

# Citation

To be determined

 
