import argparse
import os

from src import miniprot_pipeline as mp

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="This script automates gene prediction via miniprot.")

    # Required arguments
    parser.add_argument('-d', required=True, help="Output path")
    parser.add_argument('-r', required=True, help="Path to directory containing reference nucleotide sequences.")
    parser.add_argument('-f', required=True, help="Path to directory containing genome assemblies in fasta format")

    # Optional argument
    parser.add_argument('-i', required=False, help="Path to directory containing miniprot genome indexes")

    args = parser.parse_args()
    download_path = mp.make_unique_directory(os.path.join(args.d, "pipeline_output"))
    reference_sequences = args.r
    master_genome_fastas = args.f

    # miniprot indexing or not
    if not args.i:
        master_genome_indexes = mp.miniprot_index(master_genome_fastas, download_path)
    else:
        master_genome_indexes = args.i

    # reference sequences conversion from nucleotide to protein
    reference_sequences_prot = mp.reference_sequences_nuc_to_protein(reference_sequences, download_path)

    # miniprot alignment and prediction
    master_miniprot_alignments = mp.miniprot_align(reference_sequences_prot, master_genome_indexes, download_path)

    # agat sequence extraction
    master_agat_output = mp.agat_gff_cds_extraction(master_genome_fastas, master_miniprot_alignments, download_path)

    # combine predictions into single fasta files
    combined_predictions = mp.combine_agat_files_by_gene(master_agat_output, download_path)

    # add reference sequences to prediction fasta files
    mp.add_references_to_combined_files(reference_sequences, download_path)

