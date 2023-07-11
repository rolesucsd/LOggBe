import os
from ete3 import Tree


def calculate_branch_lengths(dir_path: str, output: str) -> None:
    """
    Calculate branch lengths from IQTREE output files and write the rates to a file.

    Args:
        dir_path: Directory path containing the IQTREE output files.
        output: Output directory for writing the rates.
    """
    file_names = os.listdir(dir_path)
    rates = []

    for file_name in file_names:
        if not file_name.endswith(".treefile"):
            continue

        basename = file_name.split(".", 1)[0]
        tree = Tree(os.path.join(dir_path, file_name))

        total_length = 0
        for node in tree.traverse():
            total_length += node.dist

        num_branches = len(tree.get_descendants())
        avg_length = total_length / num_branches

        rates.append((basename, avg_length))

    rates.sort(key=lambda x: x[1], reverse=True)

    with open(os.path.join(output, "rates.txt"), "w") as f:
        for rate in rates:
            f.write(f"{rate[0]}\t{rate[1]}\n")