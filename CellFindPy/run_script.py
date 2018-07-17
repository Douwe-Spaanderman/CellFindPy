import sys
import os
import time
import argparse
import scanpy.api as sc
import pandas as pd

import Lib_ as cf

def main(data_folder, output_name, output_folder, min_genes=200, min_cells=3, max_genes=7000, 
	mito_cutoff=False, n_pcs=10, initial_resolution=0.1, resolution_steps=0.05, threshold=10, subclustering=True, 
	subclustering_steps=0.2, save=True):
	#Check if Directory already exist and make new Directory
	if os.path.isdir('{}{}'.format(output_folder, output_name)) == True:
		print('WARNING')
		print('Output folder is already present')
		print('Continuing might overright files in folder')
		print('WARNING')

		input("Press Enter to continue...")
	else:
		os.mkdir('{}{}'.format(output_folder, output_name))

	#Setup adata
	adata = cf.preprocessing(data_folder, output_name, output_folder, min_genes, min_cells, max_genes, mito_cutoff)
	adata_raw = adata
	adata_raw.X = adata_raw.X.toarray()

	adata = cf.further_preprocessing(adata, n_pcs)

	#Get resolution
	df_initial_cluster, resolution = cf.find_resolution(adata, adata_raw, output_folder, output_name, initial_resolution, resolution_steps, threshold)

	#Calculate louvain using
	sc.tl.louvain(adata, random_state=100, resolution=resolution, flavor='vtraag')

	#Further subclustering
	if subclustering == True:
		adata = cf.subclustering(adata, adata_raw, output_folder, output_name, 
								initial_resolution, subclustering_steps, threshold, subclustering)
		sc.tl.tsne(adata, random_state=1, n_pcs=10, use_fast_tsne=True, learning_rate=200)
		sc.pl.tsne(adata, color='louvain')

		#FindMarkerGenes for each group
		df_all_clusters = cf.findmarker_genes(adata, adata_raw)
		if save == True:
			df_all_clusters.to_csv('{}{}/all_clusters_matrix.csv'.format(output_folder, output_name), sep=',')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="CellFindPy for RNA seq cell clustering")
	parser.add_argument("data_folder", help="path of files")
	parser.add_argument("output_folder", help="path to save file")
	parser.add_argument("output_name", help="output name")
	parser.add_argument("-min_g", "--min_genes", nargs="?", default="200", 
						help="Minimum amount of genes required for a cell to be valid")
	parser.add_argument("-min_c", "--min_cells", nargs="?", default="3", 
						help="Minimum amount of cells required for a gene to be valid")
	parser.add_argument("-max_g", "--max_genes", nargs="?", default="7000", 
						help="Maximum amount of genes permitted for a cell to be valid")
	parser.add_argument("-mito_c", "--mito_cutoff", nargs="?", default=False, 
						help="Percentage of genes permitted to be assigned to mitochondrial genes in a cell")
	parser.add_argument("-n_p", "--number_pcs", nargs="?", default="10", 
						help="Amount of principal component analysis to use for neighbor graph calculation")
	parser.add_argument("-i_r", "--initial_resolution", nargs="?", default="0.1", 
						help="The initial resolution to find clusters with")
	parser.add_argument("-r_s", "--resolution_steps", nargs="?", default="0.05", 
						help="Steps to itterate cluster finding with, a smaller resolution step helps more accurate results, however increases script run time")
	parser.add_argument("-t", "--threshold", nargs="?", default="10", 
						help="threshold of genes to be different between clusters")
	parser.add_argument("-sub", "--subclustering", nargs="?", default=True, 
						help="Option if you want to do subclustering")
	parser.add_argument("-sub_s", "--subclustering_steps", nargs="?", default="0.2", 
						help="steps similair to resolution_steps, but for subclustering")
	parser.add_argument("-sa", "--save", nargs="?", default="True", 
						help="Option to save the data, if return False only the tsne are going to be shown")
	args = parser.parse_args()

	data_folder = args.data_folder
	output_folder = args.output_folder
	output_name = args.output_name
	min_genes = float(args.min_genes)
	min_cells = float(args.min_cells)
	max_genes = float(args.max_genes)
	mito_cutoff = args.mito_cutoff
	n_pcs = int(args.number_pcs)
	initial_resolution = float(args.initial_resolution)
	resolution_steps = float(args.resolution_steps)
	threshold = float(args.threshold)
	subclustering = args.subclustering
	subclustering_steps = float(args.subclustering_steps)
	save = args.save

	print(sys.version)

	print('started')
	start = time.time()

	main(data_folder, 
		output_name, 
		output_folder, 
		min_genes, 
		min_cells, 
		max_genes, 
		mito_cutoff, 
		n_pcs, 
		initial_resolution, 
		resolution_steps, 
		threshold, 
		subclustering, 
		subclustering_steps, 
		save)

	end = time.time()
	m, s = divmod(end-start, 60)
	h, m = divmod(m, 60)
	print('script is completed')
	print("completed in %d:%02d:%02d" % (h, m, s))
	print('all output files are located at {}{}/'.format(output_folder, output_name))