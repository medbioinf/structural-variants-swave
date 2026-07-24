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

import sys
import argparse

from swave.extract_alleles import (
    load_nodes_from_gfa_fasta,
    parse_minigraph_bed_to_snarls,
    extract_and_write_alleles_to_fasta
)


def main():
    parser = argparse.ArgumentParser(description="Extracts structural variant alleles from a Pangenome Graph and sample BED files.")
    
    parser.add_argument("--gfa_fasta", required=True, help="Path to the GFA FASTA file containing node sequences.")
    parser.add_argument("--bed", required=True, help="Path to the sample minigraph --call BED file.")
    parser.add_argument("--sample_id", required=True, help="Name/ID of the sample being processed.")
    parser.add_argument("--output", default="output.fa", help="Output FASTA file for extracted alleles (default: output.fa).")
    
    # Optional parameters
    parser.add_argument(
        "--spec_snarl",
        default=None,
        help="Specific snarl ID to process (e.g. '>s1>s3') for debugging/analysis. If not provided, all snarls will be extracted."
    )
    parser.add_argument(
        "--force_reverse",
        action="store_true",
        help="Enable original Swave inversion detection and rescue logic for reversed contigs."
    )
    parser.add_argument(
        "--remove_small", 
        action="store_true", 
        help="Filter out small snarls/variants below the minimum SV size threshold."
    )
    parser.add_argument(
        "--min_sv_size", 
        type=int, 
        default=50, 
        help="Minimum size (in base pairs) for a variant to be considered a structural variant (default: 50)."
    )
    
    options = parser.parse_args()
    
    nodes_dict, fasta_index = load_nodes_from_gfa_fasta(options.gfa_fasta)
    
    snarls_dict = {}
    parse_minigraph_bed_to_snarls(options.bed, options.sample_id, snarls_dict, nodes_dict, options)
    
    extract_and_write_alleles_to_fasta(snarls_dict, fasta_index, options.output)
    
    sys.exit(0)

if __name__ == "__main__":
    main()
