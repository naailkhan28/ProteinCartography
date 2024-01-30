#!/usr/bin/env python
import argparse

import pandas as pd

# only import these functions when using import *
__all__ = ["aggregate_features"]


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", required=True, help="Path to input mmseqs2 .a3m alignment file"
    )
    parser.add_argument(
        "-l", "--list", required=True, help="Path to output .txt file where UniProt IDs will be saved"
    )
    parser.add_argument(
        "-t", "--hits", required=True, help="Path to output .csv file where hits will be saved"
    )
    parser.add_argument(
        "-e", "--max_evalue", required=True, help="Maximum E-Value cutoff to use to filter mmseqs2 hits"
    )

    args = parser.parse_args()

    return args


def process_mmseqs_hits(mmseqs_hits_file, max_evalue, output_list, output_hits):

    records = []

    with open(mmseqs_hits_file, "r") as infile:
        for line in infile:
            if line.startswith(">UniRef100"):
                line = line.split()
                records.append({
                    "UniRef ID": line[0].split("_")[1],
                    "Score": int(line[1]),
                    "Identity": float(line[2]),
                    "E-Value": float(line[3])
                })

    mmseqs_df = pd.DataFrame(records)
    mmseqs_df = mmseqs_df[mmseqs_df["E-Value"] <= max_evalue]
    mmseqs_df.to_csv(output_hits, index=False)
    #mmseqs_df.to_csv(f"{output_path}/{protid}_mmseqs_hits_processed.csv", index=False)
    mmseqs_uniprot_ids = mmseqs_df["UniRef ID"].unique()

    #with open(f"{output_path}/{protid}_mmseqs_uniprot_ids.txt", "w") as outfile:
    with open(output_list, "w") as outfile:
        outfile.write("\n".join(mmseqs_uniprot_ids))
        outfile.write("\n")


# run this if called from the interpreter
def main():
    args = parse_args()
    input_file = args.input
    listfile = args.list
    hitsfile = args.hits
    max_evalue = float(args.max_evalue)

    process_mmseqs_hits(input_file, max_evalue, listfile, hitsfile)


# check if called from interpreter
if __name__ == "__main__":
    main()
