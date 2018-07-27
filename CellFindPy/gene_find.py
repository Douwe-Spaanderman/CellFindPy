import sys
import os
import argparse
import pandas as pd
import numpy as np
import scanpy.api as sc

def main(
	data_folder,
	gene1,
	gene2=False,
	gene3=False,
	gene4=False
	):
	'''
	Temporarily method for finding specific gene in dataset.

	NOTE! Combine this into run_script
	'''

	adata = sc.read('{}/adata_dataset.h5ad'.format(data_folder))

	if gene2 == False:
		sc.pl.tsne(adata, color=[gene1])
	elif gene3 == False:
		sc.pl.tsne(adata, color=[gene1, gene2])
	elif gene4 == False:
		sc.pl.tsne(adata, color=[gene1, gene2, gene3])
	else:
		sc.pl.tsne(adata, color=[gene1, gene2, gene3, gene4])

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="Temporarily method of finding specific genes"
		)
	parser.add_argument(
		"data_folder",
		help="Path to data files"
		)
	parser.add_argument(
		"gene_of_interest1",
		help="gene you want to investigate")
	parser.add_argument(
		"gene_of_interest2",
		nargs='?',
		default=False,
		help="gene you want to investigate")
	parser.add_argument(
		"gene_of_interest3",
		nargs='?',
		default=False,
		help="gene you want to investigate")
	parser.add_argument(
		"gene_of_interest4",
		nargs='?',
		default=False,
		help="gene you want to investigate")
	args = parser.parse_args()

	data_folder = args.data_folder
	gene1 = args.gene_of_interest1
	gene2 = args.gene_of_interest2
	gene3 = args.gene_of_interest3
	gene4 = args.gene_of_interest4

	if data_folder.endswith('/'):
		data_folder = data_folder[:-1]
	
	main(
		data_folder,
		gene1,
		gene2,
		gene3,
		gene4
		)

	print('tsne plotted')