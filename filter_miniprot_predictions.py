import argparse
import os

from src import miniprot_pipeline as mp

if __name__ == "__main__":


    parser = argparse.ArgumentParser(description="This script automates gene prediction via miniprot.")

    # Required arguments
    parser.add_argument('-m', required=True, help="Filter mode: 'length' or 'score'")
    parser.add_argument('-d', required=True, help="Output path")
    parser.add_argument('-r', required=True, help="Path to directory containing reference nucleotide sequences.")
    parser.add_argument('-c', required=True, help="Path to directory containing gene predictions as combined fasta files.")

    # Optional arguments
    parser.add_argument('-u', required=False, help="Upper length filter threshold (if mode is 'length'")
    parser.add_argument('-l', required=False, help="Lower length filter threshold (if mode is 'length'")
    parser.add_argument('-s', required=False, help="Raw score threshold (if mode is 'score'")

    args = parser.parse_args()
    mode = args.m
    download_path = args.d
    reference_sequences = args.r
    combined_fastas = args.c

    # miniprot indexing or not
    if mode == "length":

        if not args.u:
            upper_length_filter = 1.2
        else:
            upper_length_filter = args.u

        if not args.l:
            lower_length_filter = 0.8
        else:
            lower_length_filter = args.l

        output = mp.seq_filter(mode, reference_sequences, combined_fastas, download_path,
                                              lower_length_bound=lower_length_filter,
                                              upper_length_bound=upper_length_filter)
        print("Filtered output is available at: " + output)

    elif mode == "score":

        if not args.s:
            score_filter = 0.5
        else:
            score_filter = args.s

        output = mp.seq_filter(mode, reference_sequences, combined_fastas, download_path, score_threshold=score_filter)

    else:
        print("Incorrect mode entered")
