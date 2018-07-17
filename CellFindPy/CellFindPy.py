import sys

print(sys.version)

import argparse, time, itertools, math

import numpy as np
import pandas as pd
import scanpy.api as sc
import scipy.stats as ss
#import cellfindpy as cf

#R import for P-value adjust
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

stats = importr('stats')

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis')
sc.logging.print_versions()
results_file = '/Users/twardlab/Desktop/fhuCoch_15W_326P/hg19/pmbc3k.h5ad'
results_file_csv = '/Users/twardlab/Desktop/fhuCoch_15W_326P/hg19/Test.csv'
results_file_raw = '/Users/twardlab/Desktop/fhuCoch_15W_326P/hg19/pmbc3k_raw.h5ad'

path = '/Users/twardlab/Desktop/fhuCoch_15W_326P/hg19/'
adata = sc.read(path + 'matrix.mtx', cache=True).T  # transpose the data
adata.var_names = pd.read_csv(path + 'genes.tsv', header=None, sep='\t')[1]
adata.obs_names = pd.read_csv(path + 'barcodes.tsv', header=None)[0]

adata.var_names_make_unique()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary, as X is sparse - it transform to a dense array after summing
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

adata = adata[adata.obs['n_genes'] < 7000, :]
#adata = adata[adata.obs['percent_mito'] < 0.05, :]
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

filter_result = sc.pp.filter_genes_dispersion(
    adata.X, min_mean=0.0125, max_mean=3, min_disp=0.5)

#sc.pl.filter_genes_dispersion(filter_result)

adata_raw = adata
adata_raw.X = adata_raw.X.toarray()

sc.pp.log1p(adata)

adata = adata[:, filter_result.gene_subset]

sc.pp.regress_out(adata, ['n_counts'])

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, n_comps=40)

adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
sc.pp.neighbors(adata, n_pcs=10)

#Set items for res find
n = 0
do = True

#Find initial Resolution: AVG_DIFF NOT SAME!!! While DATA is the same || NAMES SUCK
#https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
while do == True:
	df_sub_cluster = pd.DataFrame({'A' : []})
	resolution_it = 0.125*n
	resolution = 0.125+resolution_it
	print(resolution)

	adata_n = adata

	sc.tl.louvain(adata_n, random_state=100, resolution=resolution, flavor='vtraag')

	if len(adata_n.obs.louvain.unique()) == 1:
		print('{} gives only one cluster'.format(resolution))
		if resolution >= 0.7:
			print('still no clusters found at {}'.format(resolution))
			print('aborting script')
			sys.exit()

	else:
		#test if threshold is met
		threshold = 10
		nn = 0
		pseudocount = 1

		adata_c = adata_raw[adata_n.obs_names]
		adata_c.obs['louvain'] = adata_n.obs['louvain']

		#First Make Matrix Similair to Seurat FindMarker
		unique_clusters = adata_c.obs.louvain.unique()

		for i in itertools.repeat(None, len(unique_clusters)):
			nnn = 0
			wilcoxon_list = []
			adata_nn = adata_c[adata_c.obs['louvain'] == '{}'.format(nn)]
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
			if df_sub_cluster.empty == True:
					df_sub_cluster = pd.DataFrame(sub_cluster_array, 
												  index=adata_nn.var_names, 
												  columns=['mean_{}'.format(nn), 
												  'standard_deviation_{}'.format(nn), 
												  'Avg_diff_{}'.format(nn), 
												  'P_value_{}'.format(nn)]
												  )

			else:
				df_2 = pd.DataFrame(sub_cluster_array, 
									index=adata_nn.var_names, 
									columns=['mean_{}'.format(nn), 
									'standard_deviation_{}'.format(nn), 
									'Avg_diff_{}'.format(nn), 
									'P_value_{}'.format(nn)]
									)

				df_sub_cluster = pd.concat([df_sub_cluster, df_2], axis=1, sort=True)

			nn += 1

		#Save Dataframe as copy
		df_sub_cluster_copy = df_sub_cluster

		#Now use Matrix to check if cluster is real
		df_cluster_check = df_sub_cluster[df_sub_cluster.columns[df_sub_cluster.columns.str.startswith('Avg') | df_sub_cluster.columns.str.startswith('P')]]

		avg_diff_list = []
		p_value_list = []
		for columns in df_cluster_check:
			if columns.startswith('Avg'):
				avg_diff_list.append(sum(df_cluster_check[columns] > math.log(2)) > threshold)
			elif columns.startswith('P'):
				p_value_list.append(sum(df_cluster_check[columns] < 10e-5) > threshold)

		cluster_test = (sum(avg_diff_list)) >= (len(avg_diff_list) - 1) and (sum(p_value_list)) >= (len(p_value_list) - 1)
		print(cluster_test)

		if cluster_test == True:
			print('real clusters at {}'.format(resolution))
			do = True
		elif cluster_test == False:
			print('to far at {}'.format(resolution))
			do = False

		#df_sub_cluster_copy.to_csv('/Users/twardlab/Desktop/fhuCoch_15W_326P/hg19_{}_test_matrix.csv'.format(resolution), sep=',')
		#sys.exit()
	n += 1

#Optional Save of the final cluster dataframe
df_sub_cluster_copy.to_csv('/Users/twardlab/Desktop/fhuCoch_15W_326P/hg19_{}_matrix.csv'.format(resolution), sep=',')

#Calculated first tsne and plot it
sc.tl.louvain(adata_n, random_state=100, resolution=(resolution-0.1), flavor='vtraag')
sc.tl.tsne(adata, random_state=1, n_pcs=10, use_fast_tsne=True, learning_rate=200)
sc.pl.tsne(adata, color='louvain')

#Subdivide clusters into clusters if genes vary by 10
unique_clusters = adata.obs.louvain.unique()

n = 1
for i in itertools.repeat(None, len(unique_clusters)):
	print(n)
	do = True

	adata_n = adata[adata.obs['louvain'] == '{}'.format(n)]

	while do == True:
		df_sub_cluster = pd.DataFrame({'A' : []})
		resolution_it = 0.125*nn
		resolution = 0.125+resolution_it

		sc.pp.neighbors(adata_n, n_pcs=10)
		sc.tl.louvain(adata_n, random_state=100, resolution=resolution, flavor='vtraag')

		sc.tl.tsne(adata_n, random_state=1, n_pcs=10, use_fast_tsne=True, learning_rate=200)

		if len(adata_n.obs.louvain.unique()) == 1:
			resolution_it += 0.1

		else:
			#Test if threshold is met
			threshold = 10

			#Get raw data with only cells in subcluster -> add subcluster
			adata_c = adata_raw[adata_n.obs_names]
			adata_c.obs['louvain'] = adata_n.obs['louvain']

			#Calculate mean/sd/avg_def/Pval for both groups and create dataframe
			unique_sub_clusters = adata_c.obs.louvain.unique()
			nn = 0

			#Check total sd and mean
			standard_deviation_total = np.mean(adata_c.X, axis=0)
			mean_total = np.mean(adata_c.X, axis=0)

			for ii in itertools.repeat(None, len(unique_sub_clusters)):

				adata_nn = adata_c[adata_c.obs['louvain'] == '{}'.format(nn)]
				#Use adata_nn.X tot calculate mean and sd for every gene
				standard_deviation = np.std(adata_nn.X, axis=0)
				mean = np.mean(adata_nn.X, axis=0)

				#Maybe faster in array instead of dataframe
				sub_cluster_array = np.column_stack((mean, standard_deviation))

				if df_sub_cluster.empty == True:
					df_sub_cluster = pd.DataFrame(sub_cluster_array, index=adata_nn.var_names, 
						columns=['mean_{}'.format(nn), 'standard_deviation_{}'.format(nn), 'Avg_diff_{}'.format(nn)])

				else:
					df_2 = pd.DataFrame(sub_cluster_array, index=adata_nn.var_names, 
						columns=['mean_{}'.format(nn), 'standard_deviation_{}'.format(nn), 'Avg_diff_{}'.format(nn)])

					df_sub_cluster = pd.concat([df_sub_cluster, df_2], axis=1)

				df_sub_cluster.to_csv('/Users/twardlab/Desktop/fhuCoch_15W_326P/hg19_{}_matrix.csv'.format(n), sep=',')

				nn += 1

			do = False

	else:
		print('{} is further clustered with resolution={}'.format(n, resolution))

	n += 1

	sc.pl.tsne(adata_n, color='louvain')

	sys.exit()




