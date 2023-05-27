import pandas as pd
import os
import subprocess
from Bio import Phylo
import sys

"""
Function to create necessary directatories

Parameters
____________

output_directory: defined place where faithpdsummaries will be saved

biom_dir: folder where sorting biom tables are created
fpd_dir: folder where faith pd outputs are stored
filter_tree_dir: folder where filtered newick trees are stored


Returns
____________
Three folders in which to save results of faith pd and Biom text files

"""

def create_directories(output_directory):
    
    biom_dir = os.path.join(output_directory, "Biom")
    fpd_dir = os.path.join(output_directory, "FPD")
    filter_tree_dir = os.path.join(output_directory, "Filter_tree")

    os.makedirs(biom_dir, exist_ok=True)
    os.makedirs(fpd_dir, exist_ok=True)
    os.makedirs(filter_tree_dir, exist_ok=True)


"""
Function to run faith pd

Parameters
____________

output_directory: defined place where faithpdsummaries will be saved

group: folder where sorting biom tables are created

input_directory: folder where faith pd outputs are stored

output_directory: output for faithpd summary

faithpd_test2.sh: command line script that runs faith pd

Returns
____________

A summary text file of the results

"""



def faith_pd(input_directory, group, output_directory):
    create_directories(output_directory)

    file_names = os.listdir(input_directory)

    """
    makes a list of file names at the input directory with the .treefile extension

    retrieves the basename by splitting on the first period

    creates a path to the output directories needed 
    - including the biom.txt file, the .biom file, the _pd.txt file, and the _filter.nwk file

    if in filtering too many sequences were removed then txt_pd_file_path wont exist

    """


    for file_name in file_names:
        
        if not file_name.endswith(".treefile"):
            continue

        file_path = os.path.join(input_directory, file_name)
        tree = Phylo.read(file_path, "newick")

        base_name = file_name.split(".", 1)[0]
        txt_file_path = os.path.join(output_directory, "Biom", base_name + ".txt")
        biom_file_path = os.path.join(output_directory, "Biom", base_name + ".biom")
        txt_pd_file_path = os.path.join(output_directory, "FPD", base_name + "_pd.txt")
        tree_path = os.path.join(output_directory, "Filter_tree", base_name + "_filter.nwk")

        """
        if txt_pd_file_path exists then a list of the tip labels is made

        from this the original biom table is read in as a dataframe, then subsetted by the strains in the filtered tree

        saves the new subsetted df to a csv

        additionally the tree is subsetted and written to a .nwk tree

        """



        if os.path.exists(txt_pd_file_path):
            continue

        tip_labels = [clade.name for clade in tree.get_terminals()]

        df = pd.read_csv(group, sep="\t", encoding='latin-1')

        df = df[["Sample", "Fecal"]]

        df = df[df.iloc[:, 1].isin([0, 1])]

        subset_df = df[df['Sample'].isin(tip_labels)]

        subset_tree = tree.common_ancestor(subset_df['Sample'].tolist())

        subset_df.to_csv(txt_file_path, sep='\t', index=False)

        sys.setrecursionlimit(10**6)
        Phylo.write(subset_tree, tree_path, format="newick")

        """

        the faithpd_test2.sh script is then called, creating the .biom file at its new directory

        then the tree for a given tree has faith pd run on it based on the .biom table

        results are saved to txt_pd_file_path

        """

        try:
            # Run the subprocess for faithpd_test2.sh
            subprocess.call(['Scripts/Bash/faithpd_test2.sh', txt_file_path, biom_file_path, txt_pd_file_path, tree_path])
        except subprocess.CalledProcessError as e:
            print(f"Error occurred during subprocess call for {file_name}: {e}")
            continue
