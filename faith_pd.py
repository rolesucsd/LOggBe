import pandas as pd
import os
import subprocess
from Bio import Phylo
import sys


def create_directories(output_directory):
    # Create the directories if they don't exist
    biom_dir = os.path.join(output_directory, "Biom")
    fpd_dir = os.path.join(output_directory, "FPD")
    filter_tree_dir = os.path.join(output_directory, "Filter_tree")

    os.makedirs(biom_dir, exist_ok=True)
    os.makedirs(fpd_dir, exist_ok=True)
    os.makedirs(filter_tree_dir, exist_ok=True)

def faith_pd(input_directory, group, output_directory):
    create_directories(output_directory)

    file_names = os.listdir(input_directory)

    for file_name in file_names:
        # Skip any non-IQTREE files
        if not file_name.endswith(".treefile"):
            continue

        # Read the tree from the IQTREE output file
        file_path = os.path.join(input_directory, file_name)
        tree = Phylo.read(file_path, "newick")

        # Create file paths for output files
        base_name = file_name.split(".", 1)[0]  # Extract base name up to the first "."
        txt_file_path = os.path.join(output_directory, "Biom", base_name + ".txt")
        biom_file_path = os.path.join(output_directory, "Biom", base_name + ".biom")
        txt_pd_file_path = os.path.join(output_directory, "FPD", base_name + "_pd.txt")
        tree_path = os.path.join(output_directory, "Filter_tree", base_name + "_filter.nwk")

        if os.path.exists(txt_pd_file_path):
            continue

        # Get all the tip labels from the tree
        tip_labels = [clade.name for clade in tree.get_terminals()]

        # Read the group biom file into a DataFrame
        df = pd.read_csv(group, sep="\t", encoding='latin-1')

        # Subset the dataframe to only include sample and fecal
        df = df[["Sample", "Fecal"]]

        # Filter the DataFrame to remove rows with values other than 0 or 1 in the second column
        df = df[df.iloc[:, 1].isin([0, 1])]

        # Subset the DataFrame based on tip labels
        subset_df = df[df['Sample'].isin(tip_labels)]

        # Subset the tree to only include tips present in the 'Sample' column
        subset_tree = tree.common_ancestor(subset_df['Sample'].tolist())

        # Write the subset DataFrame to a .txt file
        subset_df.to_csv(txt_file_path, sep='\t', index=False)

        # Write the subset tree to a new Newick file
        sys.setrecursionlimit(10**6)  # Increase recursion limit if needed
        Phylo.write(subset_tree, tree_path, format="newick")

        try:
            # Run the subprocess for faithpd_test2.sh
            subprocess.call(['Scripts/Bash/faithpd_test2.sh', txt_file_path, biom_file_path, txt_pd_file_path, tree_path])
        except subprocess.CalledProcessError as e:
            print(f"Error occurred during subprocess call for {file_name}: {e}")
            continue
