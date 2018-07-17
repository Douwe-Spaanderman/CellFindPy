#!/bin/bash
printf '\e[8;50;100t'

echo "===================================================================================================="
echo "
				#######################
 				#Welcome at CellFindPy#																
 				#######################
"
echo "		CellFindPy uses Single Cell RNA sequencing to identify and 
		distinquisch cells based on clustering provided by Scanpy															
"
echo "		CellFindPy is a collaboration of Douwe Spaanderman, Kevin Yu, 
	Katharine Lee and Aaron Tward from the University of California: San Francisco
"

echo "----------------------------------------------------------------------------------------------------"

echo "Make sure to have this setup file located in the same folder as the run script 
	"
echo "To initialize CellFindPy, locate the folder where the RNA sequencing data files are located"
read -p 'File location: ' file_location

echo "Now please set the location where you want the output data to be stored aswell as giving it an 
unique name"
read -p 'Save location: ' store_folder
read -p 'Save name: ' save_name

read -r -p "Do you want to set multiple different parameters or use the standard parameters? [y/N]" response1
	case "$response1" in
    [yY][eE][sS]|[yY]) 
        echo "Minimum amount of genes required for a cell to be valid (default=200)"
        read -p 'Min_genes: ' min_g
        echo "Minimum amount of cells required for a gene to be valid (default=3)"
        read -p 'Min_cells: ' min_c
        echo "Maximum amount of genes permitted for a cell to be valid (default=7000)"
        read -p 'Max_genes: ' max_g
        echo "Percentage of genes permitted to be assigned to mitochondrial genes in a cell (default=False)"
        read -p 'Mito_cutoff: ' mito_c
        echo "Amount of principal component analysis to use for neighbor graph calculation (default=10)"
        read -p 'n_pcs: ' n_pcs
        echo "The initial resolution to find clusters with (default=0.1)"
        read -p 'initial_resolution: ' i_r
        echo "Steps to itterate cluster finding with, a smaller resolution step helps more accurate results, 
        however increases script run time (default=0.05)"
        read -p 'resolution_steps: ' r_s
        echo "threshold of genes to be different between clusters (default=10)"
        read -p 'gene_threshold: ' t
        echo "Option if you want to do subclustering (default=True)"
        read -p 'subclustering: ' sub
        echo "subclustering_steps. Similair to resolution_steps, but for subclustering (default=0.02)"
        read -p 'subclustering_steps: ' sub_s
        echo "Option to save the data (default=True)"
        read -p 'save: ' save
        ;;
    *)
        echo "default setting are selected:"
        echo "Min_genes = 200"
        echo "Min_cells = 3"
        echo "Max_genes = 7000"
        echo "Mito_cutoff = False"
        echo "Number_pcs = 10"
        echo "initial_resolution = 0.1"
        echo "resolution_steps = 0.05"
        echo "gene_threshold = 10"
        echo "subclustering = True"
        echo "subclustering_steps = 0.02"
        echo "Save = True"
        ;;
    esac

echo "----------------------------------------------------------------------------------------------------"
read -r -p "Do you want to initialize a virtual environment before initializing CellFindPy? [y/N]" response2
	case "$response2" in
    [yY][eE][sS]|[yY]) 
        echo "Please locate this virtual environment"
        echo "Example: /Users/UserName/virtualenvs/cellfindpy_environment "
        echo "! Note that you don't have to include /bin/activate" 
        read -p "virtualenv location: " virtualenv_location
        source "$virtualenv_location/bin/activate"
        ;;
    *)
		echo "Local environment will be used"
        ;;
    esac

read -r -p "Do you want to check if all the required packages are installed? [y/N] Note highly recommended on
first run. Also make sure to have requirements.txt one folder down if you want to install packages" response3
    case "$response3" in
    [yY][eE][sS]|[yY])
        pip install -r ../requirements.txt
        ;;
    *)
        echo "Packages weren't checked"
        ;;
    esac

echo "CellFindPy was optimized using MultiCore-TSNE from Ulyanov (2017), we highly recommend using
it aswell -> https://github.com/DmitryUlyanov/Multicore-TSNE"

echo "----------------------------------------------------------------------------------------------------"

echo "
		all parameters have been set -> initializing CellFindPy
	"

case "$response1" in
[yY][eE][sS]|[yY]) 
    python3 run_script.py -min_g $min_g -min_c $min_c -max_g $max_g -mito_c $mito_c -n_p $n_pcs -i_r $i_r -r_s $r_s -t $t -sub $sub -sub_s $sub_s -sa $save $file_location $store_folder $save_name
    ;;
*)
    python3 run_script.py $file_location $store_folder $save_name 
esac

echo "===================================================================================================="