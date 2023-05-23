import os
from ete3 import Tree

def calculate_branch_lengths(dir_path, output):
    file_names = os.listdir(dir_path)
    # List of rates of evolution
    rates = []

    # Loop over all the IQTREE output files
    for file_name in file_names:
        # Skip any non-IQTREE files
        if not file_name.endswith(".treefile"):
            continue

        # Get the basename of the file without the extension
        basename = file_name.split(".", 1)[0]  # Extract base name up to the first "."

        # Load the tree from the IQTREE output file
        tree = Tree(os.path.join(dir_path, file_name))

        # Calculate the total branch length
        total_length = 0
        for node in tree.traverse():
            total_length += node.dist

        # Calculate the number of branches
        num_branches = len(tree.get_descendants())

        # Calculate the average branch length
        avg_length = total_length / num_branches

        # Add the rate to the list along with the basename
        rates.append((basename, avg_length))

    # Sort the rates in descending order based on the average branch length
    rates.sort(key=lambda x: x[1], reverse=True)

    # Write the rates to a file
    with open(f"{output}/rates.txt", "w") as f:
        for rate in rates:
            f.write(f"{rate[0]}\t{rate[1]}\n")
