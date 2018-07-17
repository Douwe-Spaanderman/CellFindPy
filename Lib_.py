import sys
import itertools
import math
import pandas as pd
import numpy as np
import scanpy.api as sc
import scipy.stats as ss

#R import for P-value adjust
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

#sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3) for additional scanpy information
sc.settings.set_figure_params(color_map='viridis')
stats = importr('stats')


def preprocessing(data_folder, output_name, output_folder, min_genes=200, min_cells=3, max_genes=7000, mito_cutoff=False):
	"""
	Combined function for preprocessing using Scanpy. For a more complete documentation on preprocessing, please
	visit

	Input:
		data_folder = 
		output_name =
		output_folder =
		min_genes =
		min_cells =
		max_genes =
		mito_cutoff =

	Returns AnnData type

	
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

	if mito_cutoff == False:
		pass
	else:
		adata = adata[adata.obs['percent_mito'] < float(mito_cutoff), :]

	sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

	return adata

def further_preprocessing(adata, n_pcs=10):
	"""

	"""

	filter_result = sc.pp.filter_genes_dispersion(
    	adata.X, min_mean=0.0125, max_mean=3, min_disp=0.5)

	sc.pp.log1p(adata)

	adata = adata[:, filter_result.gene_subset]

	sc.pp.regress_out(adata, ['n_counts'])

	sc.pp.scale(adata, max_value=10)

	sc.tl.pca(adata, n_comps=40)

	#adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
	sc.pp.neighbors(adata, n_pcs=n_pcs)

	return adata

def find_resolution(adata, adata_raw, output_folder, output_name, initial_resolution=0.1, 
	resolution_steps=0.05, threshold=10, subclustering=False, subclustering_steps=0.02, save=True):
	"""

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
				if subclustering == False:
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
		if subclustering == False:
			resolution = resolution-resolution_steps
			print('resolution is set back to {}'.format(resolution))
		else:
			resolution = resolution-subclustering_steps
			print('resolution is set back to {}'.format(resolution))

	else:
		resolution = initial_resolution
		print('resolution is set back to {}'.format(resolution))

	if save == True:
		if subclustering == False:
			df_cluster_copy.to_csv('{}/{}/initial_cluster_matrix.csv'.format(output_folder, output_name), sep=',')

	return df_cluster_copy, resolution

def findmarker_genes(adata, adata_raw, df_cluster=pd.DataFrame({'A' : []}), pseudocount=1):
	"""

	"""
	nn = 0
	adata_c = adata_raw[adata.obs_names]
	adata_c.obs['louvain'] = adata.obs['louvain']

	#Change categorial first to numeric in order to sort correctly
	unique_clusters = sorted(pd.to_numeric(adata_c.obs.louvain.unique()))

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
		mean = np.log(np.mean(np.expm1(adata_nn.X), axis=0) + pseudocount)
		mean_total = np.log(np.mean(np.expm1(adata_ccc.X), axis=0) + pseudocount)
		avg_diff = mean - mean_total

		adata_nn.var['Avg_diff'] = avg_diff

		#Only take genes that show 0.25 dif (log_threshold)
		adata_nn = adata_nn[:, (adata_nn.var['Avg_diff'] >= 0.25) | (adata_nn.var['Avg_diff'] <= -0.25)]
		adata_ccc_wilcoxon = adata_ccc[:, list(adata_nn.var_names)]

		#Calculate values again for array creation (because genes got selected out)
		# +1 is pseudocount to match seurat
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

def cluster_check(df_cluster, resolution, threshold=10):
	"""

	"""
	#Now use Matrix to check if cluster is real
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

def subclustering(adata, adata_raw, output_folder, output_name, initial_resolution=0.1, 
	subclustering_steps=0.02, threshold=10, subclustering=True, save=True):
	"""

	"""
	n = 0
	unique_clusters = adata.obs.louvain.unique()
	do = True
	louvain_dict = dict()
	print('running through clusters to find subclusters')

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

			if save == True:
				df_sub_cluster.to_csv('{}/{}/cluster_{}_matrix.csv'.format(output_folder, output_name, n), sep=',')
				sc.pl.tsne(adata_of_n, color='louvain')

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

					if save == True:
						df_sub_cluster.to_csv('{}/{}/cluster_{}_{}_matrix.csv'.format(output_folder, output_name, n, nn), sep=',')
						sc.pl.tsne(adata_of_n, color='louvain')

					for iii in itertools.repeat(None, len(unique_sub_sub_clusters)):
						adata_sub_sub_of_n = 0
						adata_sub_sub_of_n = adata_sub_of_n[adata_sub_of_n.obs['louvain'] == '{}'.format(nnn)]
						cluster_numbers = ['{}.{}.{}'.format(n, nn, nnn)] * len(adata_sub_sub_of_n.obs_names) 
						i_dict = dict(zip(adata_sub_sub_of_n.obs_names, cluster_numbers))
						louvain_dict = {**louvain_dict, **i_dict}

						nnn += 1

				nn += 1

		n += 1

	#remove old louvain and add new one
	adata.obs = adata.obs.drop('louvain', axis=1)

	adata.obs = pd.concat([adata.obs, pd.DataFrame([louvain_dict], index=['louvain']).T], axis=1)

	return adata

