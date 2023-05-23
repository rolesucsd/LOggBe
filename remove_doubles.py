import os
import random
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
import pandas as pd

def iterate_files(directory_name):
    # Check if the directory exists
    if not os.path.exists(directory_name):
        print(f"The directory '{directory_name}' does not exist.")
        return []

    # List to store the files
    file_list = []

    # Iterate over each file in the directory
    for filename in os.listdir(directory_name):
        file_path = os.path.join(directory_name, filename)
        
        # Check if the current item is a file
        if os.path.isfile(file_path):
            # Add the file to the list
            file_list.append(file_path)

    return file_list

def get_basename_without_extension(file_path):
    # Get the base name of the file
    base_name = os.path.basename(file_path)
    
    # Remove the extension from the base name
    basename_without_extension = os.path.splitext(base_name)[0]
    
    return basename_without_extension

def remove_duplicates_from_msa(alignment):
    # Create a dictionary to store the number of occurrences for each sequence name
    sequence_dict = {}

    # Count the occurrences of each sequence name
    for record in alignment:
        sequence_name = record.id
        if sequence_name in sequence_dict:
            sequence_dict[sequence_name].append(record)
        else:
            sequence_dict[sequence_name] = [record]

    # Randomly select one sequence for each name
    unique_records = [random.choice(sequences) for sequences in sequence_dict.values()]

    return unique_records


def filter_alignment(file_path, gene_name, output):
    alignment = AlignIO.read(file_path, "fasta")

    # Rename the sequences in the alignment
    for record in alignment:
        string = record.name
        new_name = string.split(";")
        record.id = new_name[0]
        record.name = record.id
        record.description = ""

    alignment_renamed = alignment

    threshold = 50

    # Filter sequences below the threshold
    filtered_alignment = [record for record in alignment_renamed if
                           ((record.seq.count('-') / len(record.seq)) * 100) <= threshold]

    filtered_alignment = remove_duplicates_from_msa(filtered_alignment)

    file_name_new = f"{output}/{gene_name}_{'problematic' if len(filtered_alignment) < 400 else 'filtered'}.aln.fas"

    # Save the filtered alignment
    with open(file_name_new, "w") as handle:
        SeqIO.write(filtered_alignment, handle, "fasta")