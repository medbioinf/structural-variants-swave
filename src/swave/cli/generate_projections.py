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

from swave.generate_dotplots_projections import process_dotplot_bundles_to_projections


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def main():
    parser = argparse.ArgumentParser(description="Reads the dotplot bundles from a .pkl file and generates projections"
                                     "and segments them into matrices for each snarl, saving them into a compressed .pkl file.")
    
    parser.add_argument('--dotplot_bundles_pkl', required=True, help="Path to the dotplot bundles .pkl file.")
    parser.add_argument('--projections_out_prefix', required=True, help="Output prefix for the generated projections .pkl file.")
    parser.add_argument("--kmer_size", required=True, type=int, default=30, help="K-mer size used for matrix evaluations (default: 30).")
    
    # Optional parameters
    parser.add_argument(
        "--save_projections_images", 
        action="store_true", 
        help="If specified, saves the synthesized 6-subplot wave projection PNGs."
    )
    
    options = parser.parse_args()
    
    logging.info(f"Starting projection generation for {options.dotplot_bundles_pkl}...")
    
    process_dotplot_bundles_to_projections(
        dotplot_bundles_pickle_path=options.dotplot_bundles_pkl,
        projections_output_prefix=options.projections_out_prefix,
        options=options
    )
    
    logging.info(f"Projection generation for {options.dotplot_bundles_pkl} completed.")
    
    sys.exit(0)

if __name__ == "__main__":
    main()
