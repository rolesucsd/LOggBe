import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

from remove_doubles import iterate_files, get_basename_without_extension, filter_alignment
from evolution_rates import calculate_branch_lengths


def process_alignment_files(directory: str, output: str) -> None:
    """
    Process alignment files by filtering and saving the filtered alignment to the output directory.

    Args:
        directory: Directory containing the alignment files.
        output: Output directory for the filtered alignments.
    """
    files = iterate_files(directory)

    for file in files:
        basename = get_basename_without_extension(file)
        output_file = os.path.join(output, f"{basename}.aln_filtered.aln.fas")
        if os.path.exists(output_file):
            print(f"Skipping {file} - Output file already exists")
            continue

        filter_alignment(file, basename, output)


def create_trees(directory: str, num_jobs: int) -> None:
    """
    Create phylogenetic trees using alignment files.

    Args:
        directory: Directory containing the alignment files.
        num_jobs: Number of parallel jobs for creating trees.
    """
    files = iterate_files(directory)

    with ThreadPoolExecutor(max_workers=num_jobs) as executor:
        futures = []
        for file in files:
            output_file = file + ".treefile"
            if os.path.exists(output_file):
                print(f"Skipping {file} - Output file already exists")
                continue

            future = executor.submit(run_iqtree, file)
            futures.append(future)

        for future in futures:
            future.result()


def run_iqtree(file: str) -> None:
    """
    Run IQ-TREE to create phylogenetic trees.

    Args:
        file: Alignment file for creating the tree.
    """
    subprocess.call(['/panfs/roles/BF/Pangenome/iqtree.sh', file])


def main() -> None:
    """
    Main function to process alignment files and create phylogenetic trees.
    """
    parser = argparse.ArgumentParser(description='Process alignment files and create phylogenetic trees.')
    parser.add_argument('alignment_directory', help='Directory containing the alignment files')
    parser.add_argument('output_directory', help='Output directory for the filtered alignments')
    parser.add_argument('--jobs', '-j', type=int, default=10, help='Number of parallel jobs for creating trees')

    args = parser.parse_args()

    alignment_directory: str = args.alignment_directory
    output_directory: str = args.output_directory

    os.makedirs(output_directory, exist_ok=True)
    os.makedirs(os.path.join(output_directory, "Alignment"), exist_ok=True)
    os.makedirs(os.path.join(output_directory, "Evolution"), exist_ok=True)

    process_alignment_files(alignment_directory, os.path.join(output_directory, "Alignment"))
    create_trees(os.path.join(output_directory, "Alignment"), args.jobs)
    calculate_branch_lengths(os.path.join(output_directory, "Alignment"), os.path.join(output_directory, "Evolution"))


if __name__ == '__main__':
    main()
