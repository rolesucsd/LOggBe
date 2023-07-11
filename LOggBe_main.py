from remove_doubles import iterate_files, get_basename_without_extension, filter_alignment
from evolution_rates import calculate_branch_lengths
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

def process_alignment_files(directory, output):
    # Get a list of all alignment files in the directory
    files = iterate_files(directory)

    for file in files:
        # Get the basename of the file without the extension
        basename = get_basename_without_extension(file)

        # Check if the output file already exists
        output_file = os.path.join(output, f"{basename}.aln_filtered.aln.fas")
        if os.path.exists(output_file):
            print(f"Skipping {file} - Output file already exists")
            continue

        # Filter the alignment and save the filtered alignment to the output directory
        filter_alignment(file, basename, output)

def create_trees(directory, num_jobs):
    # Get a list of all alignment files in the directory
    files = iterate_files(directory)

    # Create a ThreadPoolExecutor with the specified number of jobs
    with ThreadPoolExecutor(max_workers=num_jobs) as executor:
        # Submit jobs to the executor
        futures = []
        for file in files:
            # Check if the output file already exists
            output_file = file + ".treefile"
            if os.path.exists(output_file):
                print(f"Skipping {file} - Output file already exists")
                continue

            # Submit the job to the executor
            future = executor.submit(run_iqtree, file)
            futures.append(future)

        # Wait for all jobs to complete
        for future in futures:
            future.result()

def run_iqtree(file):
    # Create phylogenetic trees
    subprocess.call(['/panfs/roles/BF/Pangenome/iqtree.sh', file])

# Specify the directory containing the alignment files
alignment_directory = "/panfs/roles/BF/Pangenome/Panaroo/aligned_gene_sequences"
group = "/panfs/roles/BF/Pangenome/group.txt"

# Specify the output directory for the filtered alignments
output_directory = "/panfs/roles/BF/Pangenome/Panaroo/Core"
os.makedirs(output_directory, exist_ok=True)
os.makedirs(output_directory+"/Alignment", exist_ok=True)
os.makedirs(output_directory+"/Evolution", exist_ok=True)

# Process the alignment files
#process_alignment_files(alignment_directory, output_directory+"/Alignment")

# Create phylogenetic trees
create_trees(output_directory+"/Alignment", 10)

# Calculate branch lengths from IQTREE output files
calculate_branch_lengths(output_directory+"/Alignment", output_directory+"/Evolution")
