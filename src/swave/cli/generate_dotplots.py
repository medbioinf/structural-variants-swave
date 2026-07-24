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
import logging
import argparse

from src.swave.generate_dotplots_projections import process_sample_alleles_to_dotplots


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def main():
    parser = argparse.ArgumentParser(description="Generates dotplot objects (REF2REF, REF2ALT, ALT2ALT)"
                                     "from pangenome graph alleles and exports them into a compressed .pkl file.")
    
    parser.add_argument('--alt_fasta', required=True, help="Path to the extracted alleles fasta (ALT).")
    parser.add_argument('--ref_fasta', required=True, help="Path to the fasta file of the original reference genome (REF).")
    parser.add_argument('--gfa_fasta', required=True, help="Path to the pangenome gfa fasta file.")
    parser.add_argument('--pkl_out_prefix', required=True, help="Output prefix for the generated .pkl file.")
    parser.add_argument('--img_out_prefix', required=True, help="Output prefix for the generated dotplot images.")
    
    # Optional parameters
    parser.add_argument('--kmer_size', type=int, default=30, help="K-mer size (default: 30).")
    parser.add_argument('--max_sv_size', type=int, default=1000000, help="Maximum SV size to process (default: 1000000).")
    parser.add_argument('--spec_path', help="Process only a specific path or snarl ID.")
    parser.add_argument('--save_dotplot_images', action='store_true', help="Saves PNG dotplot images for each snarl.")
    parser.add_argument('--skip_forward', action='store_true', help="Skip forward k-mer matches in dotplot matrix generation.")
    parser.add_argument('--skip_reverse', action='store_true', help="Skip reverse k-mer matches in dotplot matrix generation.")
    
    options = parser.parse_args()
    
    logging.info(f"Starting dotplot matrix generation for {options.alt_fasta}...")
    
    process_sample_alleles_to_dotplots(
        alt_fasta_path=options.alt_fasta,
        ref_fasta_path=options.ref_fasta,
        pangenome_gfa_fasta_path=options.gfa_fasta,
        options=options
    )
    
    logging.info(f"Dotplot generation for {options.alt_fasta} completed.")
    
    sys.exit(0)

if __name__ == "__main__":
    main()
