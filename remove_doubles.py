import os
import random
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
import pandas as pd


def iterate_files(directory_name: str) -> list:
    """
    Iterate over files in a directory and return a list of file paths.

    Args:
        directory_name: Name of the directory.

    Returns:
        List of file paths in the directory.
    """
    if not os.path.exists(directory_name):
        print(f"The directory '{directory_name}' does not exist.")
        return []

    file_list = []

    for filename in os.listdir(directory_name):
        file_path = os.path.join(directory_name, filename)

        if os.path.isfile(file_path):
            file_list.append(file_path)

    return file_list


def get_basename_without_extension(file_path: str) -> str:
    """
    Get the base name of a file without the extension.

    Args:
        file_path: Path of the file.

    Returns:
        Base name of the file without the extension.
    """
    base_name = os.path.basename(file_path)
    basename_without_extension = os.path.splitext(base_name)[0]

    return basename_without_extension


def remove_duplicates_from_msa(alignment):
    """
    Remove duplicate sequences from a multiple sequence alignment.

    Args:
        alignment: Bio.Align.MultipleSeqAlignment object.

    Returns:
        List of unique records.
    """
    sequence_dict = {}

    for record in alignment:
        sequence_name = record.id
        sequence_dict.setdefault(sequence_name, []).append(record)

    unique_records = [random.choice(sequences) for sequences in sequence_dict.values()]

    return unique_records


def filter_alignment(file_path: str, gene_name: str, output: str) -> None:
    """
    Filter a multiple sequence alignment and save the filtered alignment to a file.

    Args:
        file_path: Path of the alignment file.
        gene_name: Name of the gene.
        output: Output directory.
    """
    alignment = AlignIO.read(file_path, "fasta")

    for record in alignment:
        string = record.name
        new_name = string.split(";")
        record.id = new_name[0]
        record.name = record.id
        record.description = ""

    alignment_renamed = alignment

    threshold = 50

    filtered_alignment = [
        record
        for record in alignment_renamed
        if ((record.seq.count('-') / len(record.seq)) * 100) <= threshold
    ]

    filtered_alignment = remove_duplicates_from_msa(filtered_alignment)

    file_name_new = f"{output}/{gene_name}_{'problematic' if len(filtered_alignment) < 400 else 'filtered'}.aln.fas"

    with open(file_name_new, "w") as handle:
        SeqIO.write(filtered_alignment, handle, "fasta")
