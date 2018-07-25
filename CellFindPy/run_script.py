import sys
import os
import time
import argparse
import scanpy.api as sc
import pandas as pd

import Lib_ as cf

def main(
	data_folder,
	output_name,
	output_folder,
	data_type='mtx+tsv',
	min_genes=200,
	min_cells=3,
	max_genes=7000,
	mito_cutoff=False,
	n_pcs=10,
	initial_resolution=0.1,
	resolution_steps=0.05,
	threshold=10,
	subclustering=True,
	subclustering_steps=0.02,
	save=True
	):
	'''
	Main automated script to analyse single cell RNA sequencing data.

	CellFindPy uses Scanpy to loop through resolution and find suitable
	resolution for dataset.

	Input:
		data_folder = Path to data files
		output_name = Output name
		output_folder = Path to save folder
		data_type = Select data type, either mtx+tsv or h5ad (default mtx+tsv)
		min_genes = Minimum amount of genes required for a cell to be valid (default is set at 200)
		min_cells = Minimum amount of cells required for a gene to be valid (default is set at 3)
		max_genes = Maximum amount of genes permitted for a cell to be valid (default is set at 7000)
		mito_cutoff = Percentage of genes permitted to be assigned to mitochondrial genes in a cell 
					(default is set at False = 0)
		n_pcs = Amount of principal component analysis to use for neighbor graph calculation (default is set at 10)
		initial_resolution = The initial resolution to find clusters with (default is set at 0.1)
		resolution_steps = Steps to itterate cluster finding with, a smaller resolution step helps more accurate results, 
						however increases script run time (default is set at 0.05)
		threshold = Threshold of genes to be different between clusters (default is set at 10)
		subclustering = Option if you want to do subclustering (default is set at True = 1)
		subclustering_steps = Steps similair to resolution_steps, but for subclustering (default is set at 0.02)
		save = Option to save the data, if return False only the tsne are going to be shown
			(default set at True)

	Returns TSNE with clusters, subclusters and subsubclusters for the dataset with the best possible
	resolution. Also returns a matrix with the marker genes

	'''

	# Check if Directory already exists or make new Directory
	if os.path.isdir('{}/{}'.format(output_folder, output_name)) == True:
		print('WARNING')
		print('Output folder is already present')
		print('Continuing might overright files in folder')
		print('WARNING')

		input("Press Enter to continue...")

	else:
		os.mkdir('{}/{}'.format(output_folder, output_name))
		print('Created {}/{}'.format(output_folder, output_name))

	# Setup adata
	if data_type == 'mtx+tsv':
		adata = cf.preprocessing(data_folder, min_genes, min_cells, max_genes, mito_cutoff)
		adata_raw = adata
		sc.pp.log1p(adata_raw)
		adata_raw.X = adata_raw.X.toarray()

	elif data_type == 'h5ad':
		adata = sc.read('{}/merged_dataset.h5ad'.format(data_folder))
		adata_raw = adata
		sc.pp.log1p(adata_raw)

	else:
		raise ValueError("unrecognized data type is passed as argument")

	# Further Setup adata
	adata = cf.further_preprocessing(adata, n_pcs)

	# Get resolution
	df_initial_cluster, resolution = cf.find_resolution(adata, adata_raw, output_folder, output_name, initial_resolution, resolution_steps, threshold)

	# Calculate louvain using
	sc.tl.louvain(adata, random_state=100, resolution=resolution, flavor='vtraag')

	# Further subclustering
	if bool(subclustering) == True:
		adata = cf.subclustering(adata, adata_raw, output_folder, output_name, 
								initial_resolution, subclustering_steps, threshold, subclustering)
		sc.tl.tsne(adata, random_state=1, n_pcs=10, use_fast_tsne=True, learning_rate=200)

		#FindMarkerGenes for each group
		df_all_clusters = cf.findmarker_genes(adata, adata_raw)

		if bool(save) == True:
			df_all_clusters.to_csv('{}/{}/all_clusters_matrix.csv'.format(output_folder, output_name), sep=',')
			df_overall_data = cf.all_stats(adata, df_all_clusters)
			df_overall_data.to_csv('{}/{}/overall_matrix.csv'.format(output_folder, output_name), sep=',')
			sc.pl.tsne(adata, color='louvain', show=False, save=True)
			os.rename('./figures/tsne.pdf', '{}/{}/all_clusters_tsne.pdf'.format(output_folder, output_name))

		elif bool(save) == False:
			sc.pl.tsne(adata, color='louvain')

	elif bool(subclustering) == False:
		sc.tl.tsne(adata, random_state=1, n_pcs=10, use_fast_tsne=True, learning_rate=200)

		if bool(save) == True:
			df_initial_cluster.to_csv('{}/{}/initial_cluster_matrix.csv'.format(output_folder, output_name), sep=',')
			df_overall_data = cf.all_stats(adata, df_initial_cluster)
			df_overall_data.to_csv('{}/{}/overall_matrix.csv'.format(output_folder, output_name), sep=',')
			sc.pl.tsne(adata, color='louvain', show=False, save=True)
			os.rename('./figures/tsne.pdf', '{}/{}/initial_cluster_tsne.pdf'.format(output_folder, output_name))

		elif bool(save) == False:
			sc.pl.tsne(adata, color='louvain')

	#So scanpy save tsne is super stuppid as it let's you only save to ./figures/tnse.pdf so for every time
	#the tsne is plotted, os.rename is used to relocate it to the output folder with the wanted name
	#this last stap is to get rid of the ./figures directory
	os.rmdir('./figures')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="CellFindPy for RNA seq cell clustering"
		)
	parser.add_argument(
		"data_folder",
		help="Path to data files"
		)
	parser.add_argument(
		"output_folder",
		help="Path to save folder"
		)
	parser.add_argument(
		"output_name",
		help="Output name"
		)
	parser.add_argument(
		"data_type",
		choices=("mtx+tsv", "h5ad"),
		nargs="?",
		default="mtx+tsv",
		help="Select data type (default = mtx+tsv)"
		)
	parser.add_argument(
		"-min_g",
		"--min_genes",
		nargs="?",
		default="200",
		help="Minimum amount of genes required for a cell to be valid"
		)
	parser.add_argument(
		"-min_c",
		"--min_cells",
		nargs="?",
		default="3",
		help="Minimum amount of cells required for a gene to be valid"
		)
	parser.add_argument(
		"-max_g",
		"--max_genes",
		nargs="?",
		default="7000",
		help="Maximum amount of genes permitted for a cell to be valid"
		)
	parser.add_argument(
		"-mito_c",
		"--mito_cutoff",
		nargs="?",
		default="0",
		help="Percentage of genes permitted to be assigned to mitochondrial genes in a cell(0 = False, 1 = True)"
		)
	parser.add_argument(
		"-n_p",
		"--number_pcs",
		nargs="?",
		default="10",
		help="Amount of principal component analysis to use for neighbor graph calculation"
		)
	parser.add_argument(
		"-i_r",
		"--initial_resolution",
		nargs="?",
		default="0.1",
		help="The initial resolution to find clusters with"
		)
	parser.add_argument(
		"-r_s",
		"--resolution_steps",
		nargs="?",
		default="0.05",
		help="Steps to itterate cluster finding with, a smaller resolution step helps more accurate results, however increases script run time"
		)
	parser.add_argument(
		"-t",
		"--threshold",
		nargs="?",
		default="10",
		help="Threshold of genes to be different between clusters"
		)
	parser.add_argument(
		"-sub",
		"--subclustering",
		nargs="?",
		default="1",
		help="Option if you want to do subclustering(0 = False, 1 = True)"
		)
	parser.add_argument(
		"-sub_s",
		"--subclustering_steps",
		nargs="?",
		default="0.02",
		help="Steps similair to resolution_steps, but for subclustering"
		)
	parser.add_argument(
		"-sa",
		"--save",
		nargs="?",
		default="1",
		help="Option to save the data, if return False only the tsne are going to be shown(0 = False, 1 = True)"
		)
	args = parser.parse_args()

	data_folder = args.data_folder
	output_folder = args.output_folder
	output_name = args.output_name
	data_type = args.data_type
	min_genes = float(args.min_genes)
	min_cells = float(args.min_cells)
	max_genes = float(args.max_genes)
	mito_cutoff = float(args.mito_cutoff)
	n_pcs = int(args.number_pcs)
	initial_resolution = float(args.initial_resolution)
	resolution_steps = float(args.resolution_steps)
	threshold = float(args.threshold)
	subclustering = int(args.subclustering)
	subclustering_steps = float(args.subclustering_steps)
	save = int(args.save)

	#Some ways to clean up the arguments passed
	if data_folder.endswith('/'):
		data_folder = data_folder[:-1]
	if output_folder.endswith('/'):
		output_folder = output_folder[:-1]
	if output_name.endswith('/'):
		output_name = output_name[:-1]

	print(sys.version)

	print('started')
	start = time.time()

	main(
		data_folder,
		output_name,
		output_folder,
		data_type,
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
		save
		)

	end = time.time()
	m, s = divmod(end-start, 60)
	h, m = divmod(m, 60)
	print('script is completed')
	print("completed in %d:%02d:%02d" % (h, m, s))
	print('all output files are located at {}/{}/'.format(output_folder, output_name))