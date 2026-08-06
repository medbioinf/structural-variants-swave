#!/usr/bin/env python3
"""
Based on the Swave software package (https://github.com/songbowang125/Swave).
Originally licensed under the GPL-3.0.

Publication:
Wang, S., Xu, T., Zhang, P. & Ye, K. Population-level structural variant
characterization using pangenome graphs. Nat Genet (2026).
https://doi.org/10.1038/s41588-026-02538-6

Modified and refactored for Nextflow integration.
Copyright (c) 2026 Jonah Kapski <Jonah.Kapski@edu.ruhr-uni-bochum.de>
"""

import argparse
import sys

from swave.split_alleles import split_and_interleave_fasta
from swave.version import __version__


def main():
    parser = argparse.ArgumentParser(description="Splits fasta into N length-balanced fastas using round-robin distribution.")
    
    parser.add_argument("--version", action="version", version=f"{__version__}")
    parser.add_argument("--fasta", required=True, help="Path to input fasta file.")
    parser.add_argument("--prefix", required=True, help="Prefix for output split files.")
    parser.add_argument("--seq_per_split", type=int, default=2500, help="Number of sequences per split (default: 2500).")

    options = parser.parse_args()

    split_and_interleave_fasta(options.fasta, options.prefix, options.seq_per_split)
    sys.exit(0)


if __name__ == "__main__":
    main()
