import sys
import itertools
import os
import time
import argparse
import scanpy.api as sc
import pandas as pd

import Lib_ as cf

def main(dataset_1,
	dataset_2, 
	output_folder,
	dataset_name_addition=False,
	min_genes=200,
	min_cells=3,
	max_genes=7000,
	mito_cutoff=False,
	additional_datasets=False,
	normalize=False):
	'''
	Automated process to merge two datasets

	Note developer ability to merge more dataset than just two

	Input:
		dataset_1 = Path to dataset 1
		dataset_2 = Path to dataset 2
		output_folder = Path to save output file
		dataset_name_addition = Option to give datasets unique name (default = False)
		min_genes = Minimum amount of genes required for a cell to be valid (default is set at 200)
		min_cells = Minimum amount of cells required for a gene to be valid (default is set at 3)
		max_genes = Maximum amount of genes permitted for a cell to be valid (default is set at 7000)
		mito_cutoff = Percentage of genes permitted to be assigned to mitochondrial genes in a cell 
					(default is set at False = 0)
		additional_datasets = Path to all additional datasets as list
		normalize = Use of normalize in preprocessing (default = False)


	'''
	adata_1 = cf.preprocessing(dataset_1, min_genes, min_cells, max_genes, mito_cutoff, normalize)
	adata_2 = cf.preprocessing(dataset_2, min_genes, min_cells, max_genes, mito_cutoff, normalize)

	print(adata_1)
	if additional_datasets == False:
		merged = cf.merge_datasets(adata_1, adata_2, dataset_name_addition)
	elif isinstance(additional_datasets, (list,)):
		merged = cf.merge_datasets(adata_1, adata_2, dataset_name_addition, True)
		n = 0
		for i in itertools.repeat(None, len(additional_datasets)):
			dataset_a = additional_datasets[n]

			if dataset_a.endswith('/'):
				dataset_a = dataset_a[:-1]

			dataset_a = cf.preprocessing(dataset_a, min_genes, min_cells, max_genes, mito_cutoff, normalize)

			if n == (len(additional_datasets) - 1):
				merged = cf.merge_datasets(merged, dataset_a, dataset_name_addition, False, (int(n)+3))
			else:
				merged = cf.merge_datasets(merged, dataset_a, dataset_name_addition, True, (int(n)+3))

			n += 1
	else:
		raise ValueError("additional_datasets is neither a list or past False")

	merged.write('{}/merged_dataset.h5ad'.format(output_folder))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="Merge datasets"
		)
	parser.add_argument(
		"dataset_1",
		help="Path to dataset 1"
		)
	parser.add_argument(
		"dataset_2",
		help="Path to dataset 2"
		)
	parser.add_argument(
		"-a_d",
		"--additional_datasets",
		nargs="+",
		default=False,
		help="Path import multiple datasets (-a_d dataset3 dataset4)"
		)
	parser.add_argument(
		"output_folder",
		help="Path to save output file"
		)
	parser.add_argument(
		"-d_n_a",
		"--name_addition",
		nargs="?",
		default="0",
		help="Do you want to give a unique name to dataset (0 = False, 1 = True)"
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
	parser.add_argument("-mito_c",
		"--mito_cutoff",
		nargs="?",
		default="0",
		help="Percentage of genes permitted to be assigned to mitochondrial genes in a cell(0 = False, 1 = True)"
		)
	args = parser.parse_args()

	dataset_1 = args.dataset_1
	dataset_2 = args.dataset_2
	additional_datasets = args.additional_datasets
	output_folder = args.output_folder
	dataset_name_addition = int(args.name_addition)
	min_genes = float(args.min_genes)
	min_cells = float(args.min_cells)
	max_genes = float(args.max_genes)
	mito_cutoff = float(args.mito_cutoff)

	print(additional_datasets)

	if dataset_1.endswith('/'):
		dataset_1 = dataset_1[:-1]
	if dataset_2.endswith('/'):
		dataset_2 = dataset_2[:-1]
	if output_folder.endswith('/'):
		output_folder = output_folder[:-1]

	print(sys.version)

	print('started')
	start = time.time()

	main(dataset_1,
		dataset_2, 
		output_folder,
		dataset_name_addition,
		min_genes,
		min_cells,
		max_genes,
		mito_cutoff,
		additional_datasets
	)

	end = time.time()
	m, s = divmod(end-start, 60)
	h, m = divmod(m, 60)
	print('script is completed')
	print("completed in %d:%02d:%02d" % (h, m, s))
	print('merged dataset is located at {}/'.format(output_folder))