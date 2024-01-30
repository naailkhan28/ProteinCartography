#!/usr/bin/env python
import argparse
import os

from mmseqs2_utils import run_mmseqs2
from tests import mocks

# if necessary, mock the `run_blast` method
# see comments in `tests.mocks` for more details
if os.environ.get("PROTEINCARTOGRAPHY_WAS_CALLED_BY_PYTEST") == "true":
    mocks.mock_run_blast()


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--query", required=True, help="path to the input peptide FASTA file.")
    parser.add_argument(
        "--folder",
        required=True,
        help="path to the file where the mmseqs2 output MSA will be saved.",
    )
    parser.add_argument(
        "--file",
        required=True,
        help="path to the file where the mmseqs2 output MSA will be saved.",
    )
    parser.add_argument(
        "--protid",
        required=True,
        help="path to the file where the mmseqs2 output MSA will be saved.",
    )
    
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    run_mmseqs2(args.query, args.folder, args.file, args.protid)

if __name__ == "__main__":
    main()
