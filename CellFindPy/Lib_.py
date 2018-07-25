import sys
import itertools
import math
import os
import pandas as pd
import numpy as np
import scanpy.api as sc
import scipy.stats as ss

#R import for P-value adjust
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

sc.settings.set_figure_params(color_map='viridis')
stats = importr('stats')


def preprocessing(
	data_folder,
	min_genes=200,
	min_cells=3,
	max_genes=7000,
	mito_cutoff=False,
	normalize=True
	):
	"""
	Combined function for preprocessing using Scanpy. For a more complete documentation on preprocessing, please
	visit

	Input:
		data_folder = Path to data files
		min_genes = Minimum amount of genes required for a gene to be valid (default is set at 200)
		min_cells = Minimum amount of cells required for a gene to be valid (default is set at 3)
		max_genes = Maximum amount of genes permitted for a cell to be valid (default is set at 7000)
		mito_cutoff = Percentage of genes permitted to be assigned to  
		mitochondrial Genes in a cell (default is set at False=0)
		normalize = Normalize the Anndata object (default set at True)

	Returns AnnData type from matrix - genes - barcodes. Full documentation on AnnData can be found on github.

	"""

	#Read data and create initial AnnData Frame
	path = '{}/'.format(data_folder)
	adata = sc.read(path + 'matrix.mtx', cache=True).T  # transpose the data
	adata.var_names = pd.read_csv(path + 'genes.tsv', header=None, sep='\t')[1]
	adata.obs_names = pd.read_csv(path + 'barcodes.tsv', header=None)[0]

	adata.var_names_make_unique()

	#Filter data with min_genes per cell, max_genes per cell, min_cells per genes
	sc.pp.filter_cells(adata, min_genes=min_genes)
	sc.pp.filter_genes(adata, min_cells=min_cells)
	adata = adata[adata.obs['n_genes'] < max_genes, :]

	# add the total counts per cell as observations-annotation to adata
	adata.obs['n_counts'] = adata.X.sum(axis=1).A1

	#Create mito_genes and possible filter
	mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
	adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1

	if int(mito_cutoff) == False:
		pass
	else:
		adata = adata[adata.obs['percent_mito'] < float(mito_cutoff), :]

	#Normalize data option
	if normalize == True:
		sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

	return adata

def further_preprocessing(
	adata, 
	n_pcs=10
	):
	"""
	Further preprocessing consists of filtering genes based on dispersion,
	regress out on total counts per cell (n_counts), scale the data,
	doing principal component analysis and calculating neighbors.

	All these steps are done as described by Scanpy.

	Input:
		adata = Anndata derived from preprocessing 
		n_pcs = amount of principal component analysis used for neighbor
		calculation (default is set at 10)

	Returns preprocessed AnnData object including neigbor assignment

	"""

	# Select genes with high dispersion
	filter_result = sc.pp.filter_genes_dispersion(
    	adata.X, min_mean=0.0125, max_mean=3, min_disp=0.5)

	# Filter out these genes
	adata = adata[:, filter_result.gene_subset]

	# Regress out based on total read count per cell
	print(adata.X)
	sc.pp.regress_out(adata, ['n_counts'])
	print(adata.X)

	sc.pp.scale(adata, max_value=10)

	sc.tl.pca(adata, n_comps=40)

	# Neighbor assigment based on n_pcs number of principal component anaylsis
	sc.pp.neighbors(adata, n_pcs=n_pcs)

	return adata

def find_resolution(
	adata,
	adata_raw,
	output_folder,
	output_name,
	initial_resolution=0.1,
	resolution_steps=0.05,
	threshold=10,
	subclustering=False,
	subclustering_steps=0.02
	):
	"""
	Find Resolution uses looping through clustering to identify suitable resolution
	for adata. It starts using initial resolution for clustering -> identifies
	markers per cluster in that resolution and checks if resolution passes cluster
	check. If it does than it starts over with a higher resolution (resolution that
	just passed + resolution steps). It does this until cluster check passes False, 
	after it returns the resolution which was last suitable.

	Input:
		adata = AnnData object after further preprocessing
		adata_raw = AnnData object before further preprocessing
		output_folder = location to store output matrix
		output_name = directory name to store output matrix
		initial_resolution = The initial resolution to find clusters with (default = 0.1)
		resolution_steps = Steps to itterate cluster finding with, 
		a smaller resolution step helps more accurate results, 
		however increases script run time (default = 0.05)
		threshold = threshold of genes to be log(2) differentially expressed between clusters
		subclustering = option if you want to do subclustering(0 = False, 1 = True)
		subclustering_steps = steps similair to resolution_steps but for subclustering
	
	Returns last resolution where cluster check passed True and gives Dataframe of that
	resolution with marker genes

	"""

	n = 0
	do = True
	df_cluster = 0
	resolution_it = resolution_steps
	print('initial clustering')

	while do == True:
		#Save Dataframe from previous cluster as copy
		df_cluster_copy = df_cluster

		df_cluster = pd.DataFrame({'A' : []})
		resolution_it = resolution_steps*n
		resolution = initial_resolution + resolution_it
		print('resolution = {}'.format(resolution))

		adata_n = adata

		sc.tl.louvain(adata_n, random_state=100, resolution=resolution, flavor='vtraag')

		if len(adata_n.obs.louvain.unique()) == 1:
			print('{} gives only one cluster'.format(resolution))
			if resolution >= 0.7:
				if bool(subclustering) == False:
					print('still no initial clusters found at {}'.format(resolution))
					print('aborting script')
					sys.exit()
				else:
					print('no further subclusters found')
					resolution = initial_resolution
					do = False

		else:
			#First Make Dataframe Similair to Seurat FindMarker
			df_cluster = findmarker_genes(adata_n, adata_raw, df_cluster)
			do = cluster_check(df_cluster, resolution)

		n += 1

	#Possible Save for Dataframe
	if resolution >= initial_resolution+resolution_steps:
		if bool(subclustering) == False:
			resolution = resolution-resolution_steps
			print('resolution is set back to {}'.format(resolution))
		else:
			resolution = resolution-subclustering_steps
			print('resolution is set back to {}'.format(resolution))

	else:
		resolution = initial_resolution
		print('resolution is set back to {}'.format(resolution))

	return df_cluster_copy, resolution

def findmarker_genes(
	adata,
	adata_raw,
	df_cluster=pd.DataFrame({'A' : []}),
	pseudocount=1
	):
	"""
	Findmarker genes is inspired by Seurat's FindMarker genes, calculating the mean, 
	standard deviation, average difference expression of genes between clusters and
	p-value based on a mannwhitneyu statistical analysis. This analysis is conducted
	on the average difference expression.

	Input:
		adata = AnnData object after further preprocessing
		adata_raw = AnnData object before further preprocessing
		df_cluster = df_cluster to append data to (default = empty df)
		pseudocount = addition for average difference calculation to mimic Seurat
					(default = 1)

	Returns Dataframe with mean, standard devation, average difference expression
	and hochberg adjusted p-value for every cluster.

	"""
	nn = 0
	adata_c = adata_raw[adata.obs_names]
	adata_c.obs['louvain'] = adata.obs['louvain']

	#Change categorial first to numeric in order to sort correctly
	try:
		unique_clusters = sorted(pd.to_numeric(adata_c.obs.louvain.unique()))
	except ValueError:
		#For instance 2.0.0 can't be changed to numeric (ValueError), so this is a temporarily solution to get an unsorted big Matrix
		unique_clusters = adata_c.obs.louvain.unique()

	for i in itertools.repeat(None, len(unique_clusters)):
		nnn = 0
		wilcoxon_list = []
		adata_nn = 0
		adata_nn = adata_c[adata_c.obs['louvain'] == '{}'.format(unique_clusters[nn])]

		#Use all other cells for adata_c
		adata_cc = adata_c[list(set(list(adata_c.obs_names)) - set(list(adata_nn.obs_names))), :]

		#Min pct -> 0.25 to make calculating faster
		cluster_genes = sc.pp.filter_genes(adata_nn, min_cells=(0.25*len(adata_nn.obs_names)), copy=True)
		total_genes = sc.pp.filter_genes(adata_cc, min_cells=(0.25*len(adata_cc.obs_names)), copy=True)
		cluster_genes = list(cluster_genes.var_names)
		total_genes = list(total_genes.var_names)

		genes = list(set(cluster_genes + total_genes))

		adata_nn = adata_nn[:, genes]
		adata_ccc = adata_cc[:, genes]

		#Avg Dif
		#+1 is pseudocount to match seurat
		mean = np.log(np.mean(np.expm1(adata_nn.X), axis=0) + pseudocount)
		mean_total = np.log(np.mean(np.expm1(adata_ccc.X), axis=0) + pseudocount)
		avg_diff = mean - mean_total

		adata_nn.var['Avg_diff'] = avg_diff

		#Only take genes that show 0.25 dif (log_threshold)
		adata_nn = adata_nn[:, (adata_nn.var['Avg_diff'] >= 0.25) | (adata_nn.var['Avg_diff'] <= -0.25)]
		adata_ccc_wilcoxon = adata_ccc[:, list(adata_nn.var_names)]

		#Calculate values again for array creation (because genes got selected out)
		standard_deviation = np.std(adata_nn.X, axis=0)
		mean = np.log(np.mean(np.expm1(adata_nn.X), axis=0) + pseudocount)
		mean_total = np.log(np.mean(np.expm1(adata_ccc_wilcoxon.X), axis=0) + pseudocount)
		avg_diff = mean - mean_total

		#mannwhitneyu test is used because wilcoxon is a paired test which require equal N values (so equal amount of cells)
		for ii in itertools.repeat(None, len(adata_nn.var_names)):
			x = (adata_nn.X[:, [nnn]]).T
			y = (adata_ccc_wilcoxon.X[:, [nnn]]).T
			wilcoxon = ss.mannwhitneyu(x.ravel(), y.ravel())

			#Append all to list 
			wilcoxon_list.append(wilcoxon[1])

			nnn += 1

		#Seurat way to p_adjust by hochberg
		p_adjust = stats.p_adjust(FloatVector(wilcoxon_list), method = 'hochberg', n = len(adata_raw.var_names))
		#Convert wilcoxon list to array (faster to append to list first than appending to array)
		wilcoxon_array = np.asarray(p_adjust)

		sub_cluster_array = np.column_stack((mean, standard_deviation, avg_diff, wilcoxon_array))

		#Create DataFrame 
		if df_cluster.empty == True:
				df_cluster = pd.DataFrame(sub_cluster_array, 
											  index=adata_nn.var_names, 
											  columns=['mean_{}'.format(unique_clusters[nn]), 
											  'standard_deviation_{}'.format(unique_clusters[nn]), 
											  'Avg_diff_{}'.format(unique_clusters[nn]), 
											  'P_value_{}'.format(unique_clusters[nn])]
											  )

		else:
			df_2 = pd.DataFrame(sub_cluster_array, 
								index=adata_nn.var_names, 
								columns=['mean_{}'.format(unique_clusters[nn]), 
								'standard_deviation_{}'.format(unique_clusters[nn]), 
								'Avg_diff_{}'.format(unique_clusters[nn]), 
								'P_value_{}'.format(unique_clusters[nn])]
								)

			df_cluster = pd.concat([df_cluster, df_2], axis=1, sort=True)	

		nn += 1

	return df_cluster

def cluster_check(
	df_cluster,
	resolution,
	threshold=10
	):
	"""
	Checks in dataframe from findmarker genes if genes are > log(2) differentially
	expressed and if hochberg p-value is < 10e-5. Cluster check returns True if
	for both these cases X (default = 10) amount of genes are found.

	Note for developper -> make log(2) and 10e-5 custumizable

	Input:
		df_cluster = Dataframe from findmarker genes
		resolution = Resolution assigned to Dataframe
		threshold = Threshold of genes to be different between clusters (default = 10)

	Returns True or False on cluster check

	"""

	df_cluster_check = df_cluster[df_cluster.columns[df_cluster.columns.str.startswith('Avg') | df_cluster.columns.str.startswith('P')]]

	avg_diff_list = []
	p_value_list = []
	for columns in df_cluster_check:
		if columns.startswith('Avg'):
			avg_diff_list.append(sum(df_cluster_check[columns] > math.log(2)) > threshold)
		elif columns.startswith('P'):
			p_value_list.append(sum(df_cluster_check[columns] < 10e-5) > threshold)

	cluster_test = (sum(avg_diff_list)) >= (len(avg_diff_list) - 1) and (sum(p_value_list)) >= (len(p_value_list) - 1)

	if cluster_test == True:
		print('real clusters at {}'.format(resolution))
		do = True
	elif cluster_test == False:
		print('to far at {}'.format(resolution))
		do = False

	return do

def subclustering(
	adata,
	adata_raw,
	output_folder,
	output_name,
	initial_resolution=0.1, 
	subclustering_steps=0.02,
	threshold=10,
	subclustering=True,
	save=True
	):
	"""
	Subclustering data using initial clustering data, which is done
	similair to initial cluster looping. Option to save TSNE and
	matrix files for every single cluster and subcluster.

	Input:
		adata = AnnData object after initial clustering
		adata_raw = AnnData object before further preprocessing
		output_folder = Path to save folder
		output_name = Output name
		initial_resolution = The initial resolution to find clusters with 
						(default is set at 0.1)
		subclustering_steps = Steps similair to resolution_steps, 
						but for subclustering (default is set at 0.02)
		threshold = Threshold of genes to be different between clusters 
					(default is set at 10)
		subclustering = Option if you want to do subclustering 
					(default is set at True = 1)
		save = Option to save the data, if return False only the tsne are 
			going to be shown (default set at True)

	Returns new assigned louvain cluster to adata. Also able to save
	TSNE and matrix files for each cluster and subcluster 

	"""
	n = 0
	unique_clusters = adata.obs.louvain.unique()
	do = True
	louvain_dict = dict()
	print('running through clusters to find subclusters')

	# Loop through initial found clusters
	for i in itertools.repeat(None, len(unique_clusters)):
		print('clustering cluster {}'.format(n))
		adata_of_n = 0
		adata_of_n = adata[adata.obs['louvain'] == '{}'.format(n)]

		df_sub_cluster, resolution_sub = find_resolution(adata_of_n, adata_raw, output_folder, output_name, 
														initial_resolution, subclustering_steps, threshold, 
														subclustering, subclustering_steps)

		#Subclustering subclusters
		sc.tl.louvain(adata_of_n, random_state=100, resolution=resolution_sub, flavor='vtraag')
		unique_sub_clusters = adata_of_n.obs.louvain.unique()

		nn = 0

		if len(unique_sub_clusters) == 1:
			print('{} has no subclusters'.format(n))
			cluster_numbers = ['{}'.format(n)] * len(adata_of_n.obs_names)
			i_dict = dict(zip(adata_of_n.obs_names, cluster_numbers))
			louvain_dict = {**louvain_dict, **i_dict}

		else:
			print('{} has {} subclusters'.format(n, len(unique_sub_clusters)))

			if bool(save) == True:
				df_sub_cluster.to_csv('{}/{}/cluster_{}_matrix.csv'.format(output_folder, output_name, n), sep=',')
				sc.tl.tsne(adata_of_n, random_state=1, n_pcs=10, use_fast_tsne=True, learning_rate=200)
				sc.pl.tsne(adata_of_n, color='louvain', show=False, save=True)
				os.rename('./figures/tsne.pdf', '{}/{}/cluster_{}_tsne.pdf'.format(output_folder, output_name, n))

			# Loop through subclusters found in main cluster for subsubclustering
			for ii in itertools.repeat(None, len(unique_sub_clusters)):
				print('clustering sub_cluster {}.{}'.format(n, nn))
				adata_sub_of_n = 0
				adata_sub_of_n = adata_of_n[adata_of_n.obs['louvain'] == '{}'.format(nn)]

				df_sub_sub_cluster, resolution_sub_sub = find_resolution(adata_sub_of_n, adata_raw, output_folder, output_name, 
																		initial_resolution, subclustering_steps, threshold, 
																		subclustering, subclustering_steps)

				sc.tl.louvain(adata_sub_of_n, random_state=100, resolution=resolution_sub_sub, flavor='vtraag')
				unique_sub_sub_clusters = adata_sub_of_n.obs.louvain.unique()

				nnn = 0

				#Combine data by creating dictionary
				if len(unique_sub_sub_clusters) == 1:
					print('{}.{} has no subclusters'.format(n, nn))
					cluster_numbers = ['{}.{}'.format(n, nn)] * len(adata_sub_of_n.obs_names)
					i_dict = dict(zip(adata_sub_of_n.obs_names, cluster_numbers))
					louvain_dict = {**louvain_dict, **i_dict}
				else:
					print('{}.{} has {} subclusters'.format(n, nn, len(unique_sub_sub_clusters)))

					if bool(save) == True:
						df_sub_cluster.to_csv('{}/{}/cluster_{}_{}_matrix.csv'.format(output_folder, output_name, n, nn), sep=',')
						sc.tl.tsne(adata_sub_of_n, random_state=1, n_pcs=10, use_fast_tsne=True, learning_rate=200)
						sc.pl.tsne(adata_sub_of_n, color='louvain', show=False, save=True)
						os.rename('./figures/tsne.pdf', '{}/{}/cluster_{}_{}_tsne.pdf'.format(output_folder, output_name, n, nn))

					for iii in itertools.repeat(None, len(unique_sub_sub_clusters)):
						adata_sub_sub_of_n = 0
						adata_sub_sub_of_n = adata_sub_of_n[adata_sub_of_n.obs['louvain'] == '{}'.format(nnn)]
						cluster_numbers = ['{}.{}.{}'.format(n, nn, nnn)] * len(adata_sub_sub_of_n.obs_names) 
						i_dict = dict(zip(adata_sub_sub_of_n.obs_names, cluster_numbers))
						louvain_dict = {**louvain_dict, **i_dict}

						nnn += 1

				nn += 1

		n += 1

	# Remove old louvain data and add new one
	adata.obs = adata.obs.drop('louvain', axis=1)

	adata.obs = pd.concat([adata.obs, pd.DataFrame([louvain_dict], index=['louvain']).T], axis=1)

	return adata

def all_stats(
	adata,
	df
	):
	"""
	All stats creates an overall matrix from all the cluster together.
	Including number of cells in cluster and top 40 genes.

	Note for developer: still need average reads and genes

	Input:
		adata = AnnData object after clustering
		df = Dataframe from find markers

	Returns matrix with general data for every cluster

	"""

	#Loop through the clusters to get data for creating all_stats dataframe
	n = 0
	unique_clusters = adata.obs.louvain.unique()
	df_all_stats = pd.DataFrame({'A' : []})

	for i in itertools.repeat(None, len(unique_clusters)):

		cluster_name = unique_clusters[n]
		cluster_data = adata[adata.obs['louvain'] == '{}'.format(cluster_name)]
		number_cells = len(adata.obs_names)
		df_sorted = df.sort_values('Avg_diff_{}'.format(cluster_name))
		gene_list = df_sorted.index.values.tolist()[0:39]

		# Note should still add average reads and genes

		combined_list = [number_cells] + gene_list
		list_index = ['number_cells', 'top_1', 'top_2', 'top_3', 'top_4', 
					'top_5', 'top_6', 'top_7', 'top_8', 'top_9', 'top_10', 
					'top_11', 'top_12', 'top_13', 'top_14', 'top_15', 
					'top_16', 'top_17', 'top_18', 'top_19', 'top_20', 
					'top_21', 'top_22', 'top_23', 'top_24', 'top_25', 
					'top_26', 'top_27', 'top_28', 'top_29', 'top_30', 
					'top_31', 'top_32', 'top_33', 'top_34', 'top_35', 
					'top_36', 'top_37', 'top_38', 'top_39']
		combined_dict = dict(zip(list_index, combined_list))

		# Create DataFrame 
		if df_all_stats.empty == True:
			df_all_stats = pd.DataFrame.from_dict(combined_dict, orient='index', columns=[cluster_name])

		else:
			df_2 = pd.DataFrame.from_dict(combined_dict, orient='index', columns=[cluster_name])

			df_cluster = pd.concat([df_all_stats, df_2], axis=1, sort=True)	

		n += 1

	return df_all_stats

def merge_datasets(
	adata_1,
	adata_2,
	dataset_name_addition=False
	):
	'''
	Merge two datasets together.

	Input:
		adata_1 = AnnData dataset 1 after preprocessing
		adata_2 = AnnData dataset 2 after preprocessing
		dataset_name_addition = Option to add unique name to cells in dataset
								(default = False)

	Returns merged AnnData object

	'''
	if bool(dataset_name_addition) == True:
		# Add dataset_X to obs_names for when there are overlapping obs_names
		adata_1.obs_names = 'dataset_1_' + adata_1.obs_names
		adata_2.obs_names = 'dataset_2_' + adata_2.obs_names

		# Check if cell names are unique -> if not than stop script
	if set(adata_1.obs_names).isdisjoint(adata_2.obs_names) == False:
		sys.exit('overlapping names between datasets, aborting script. Datasets can be merged if dataset name addition is passed True')

	# Merge data provided by AnnData
	adata_merged = adata_1.concatenate(adata_2, index_unique=None, join='outer')

	# Fill var NaN values and condense columns into 1
	adata_merged.var = adata_merged.var.fillna(0)
	adata_merged.var['n_cells'] = adata_merged.var['n_cells-0'] + adata_merged.var['n_cells-1']

	# Remove additional not needed information
	del adata_merged.var['n_cells-0']
	del adata_merged.var['n_cells-1']
	del adata_merged.obs['batch']

	adata_merged.X = adata_merged.X.toarray()
	print(adata_merged.X)

	sc.pp.normalize_per_cell(adata_merged, counts_per_cell_after=1e4)

	return adata_merged