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
import os
import logging
import argparse

from swave.call_variants import process_predictions_to_tsv


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def main():
    parser = argparse.ArgumentParser(description="Reads predicted labels and dotplot bundles from a gzipped"
                                     ".pkl file and streams called variants into a structured TSV file.")
    
    parser.add_argument("--predictions_pkl", required=True, help="Path to the gzipped predictions .pkl file.")
    parser.add_argument("--projections_pkl", required=True, help="Path to the gzipped projections .pkl file.")
    parser.add_argument("--dotplots_pkl", required=True, help="Path to the gzipped dotplot bundles .pkl file.")
    parser.add_argument("--output_tsv", required=True, help="Path to the output called variants TSV file.")
    parser.add_argument("--sample_id", required=True, help="Sample ID to include in the TSV metadata column.")
    
    # Optional parameters
    parser.add_argument("--min_sv_size", type=int, default=50, help="Minimum SV length to be called.")
    parser.add_argument("--max_sv_size", type=int, default=1000000, help="Maximum SV length to be called.")
    parser.add_argument("--max_sv_comps", type=int, default=5, help="Maximum number of SV components for detailed output, otherwise, hyperCPX is used.")
    parser.add_argument("--dup_to_ins", action="store_true", default=False, help="Report duplications as insertions.")
    
    options = parser.parse_args()
    
    if not os.path.exists(options.predictions_pkl):
        logging.error(f"Predictions file not found: {options.predictions_pkl}")
        sys.exit(1)
    
    if not os.path.exists(options.projections_pkl):
        logging.error(f"Projections file not found: {options.projections_pkl}")
        sys.exit(1)
        
    if not os.path.exists(options.dotplots_pkl):
        logging.error(f"Dotplot bundles file not found: {options.dotplots_pkl}")
        sys.exit(1)
        
    logging.info(f"Processing variant calling for Sample: {options.sample_id}")
    
    process_predictions_to_tsv(
        predictions_pkl_path=options.predictions_pkl,
        projections_pkl_path=options.projections_pkl,
        dotplots_pkl_path=options.dotplots_pkl,
        output_tsv_path=options.output_tsv,
        options=options
    )

    sys.exit(0)

if __name__ == "__main__":
    main()
