
#--------------------------------------------------------------
#Chu lab
#You already Know What It Is
#Not Entirely Sure What Goes Here
#
#Full Liscense is up your butt and around the corner
#--------------------------------------------------------------




import os
from ete3 import Tree


"""
Function to Calculate Branch Length

Parameters
____________

dir_path: path to tree files, trees should be in newick form

output: location where rates.txt file should written to

Returns
____________

avg_length: the average length of the branches of a phylogenetic tree, calculated via the phylogenetic distance of the trees

rates.txt: all genes and their respective average branch length written to a .txt file
"""


def calculate_branch_lengths(dir_path, output):
    """
    Getting the list of file names

    Create a dataframe for the branch lengths to be saved in 
    
    Loop over all the IQTREE outputs as .treefiles and skip anything else

    Get the gene base name by splitting by the first "."
    """

    file_names = os.listdir(dir_path)

    rates = []

    for file_name in file_names:
        if not file_name.endswith(".treefile"):
            continue

        basename = file_name.split(".", 1)[0]

        """
        tree: read in the tree based on dir_path

        total_length: set to 0 then add the length of every node to make total branch length

        num_branches: count the number of branhces / nodes by count the number of descendents

        avg_length: divide the total length by the number of branches to get the average branch length and write it to the rates vector
        """

        tree = Tree(os.path.join(dir_path, file_name))

        total_length = 0
        for node in tree.traverse():
            total_length += node.dist


        num_branches = len(tree.get_descendants())
        avg_length = total_length / num_branches
        rates.append((basename, avg_length))

    """
    rates: data frame with gene name and average branch length, sort the rates into descending order
    - Highest phylogenetic distance genes on top

    rates.txt: Write rates to a file based upon path at output
    """
    
    rates.sort(key=lambda x: x[1], reverse=True)

    with open(f"{output}/rates.txt", "w") as f:
        for rate in rates:
            f.write(f"{rate[0]}\t{rate[1]}\n")
