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

from swave.write_vcf import process_tsv_files_to_vcf


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def main():
    parser = argparse.ArgumentParser(
        description="SWAVE_WRITE_VCF: Collects population-level variant calling TSVs "
                    "from samples and writes them into a finalized multi-sample "
                    "VCF as well as a split VCF."
    )
    
    parser.add_argument(
        "--tsv_files", required=True, nargs="+", help="Space-separated list of path(s) to the called variants TSV files."
    )
    parser.add_argument(
        "--ref_fasta", required=True, help="Path to the reference genome FASTA file (needed for VCF header contig lengths)."
    )
    parser.add_argument(
        "--output_vcf", required=True, help="Path to the output multi-sample VCF file. (The .split.vcf will be created automatically alongside)."
    )
    
    options = parser.parse_args()
    
    if not os.path.exists(options.ref_fasta):
        logging.error(f"Reference FASTA file not found: {options.ref_fasta}")
        sys.exit(1)
    
    valid_tsv_files = []
    for tsv_path in options.tsv_files:
        if os.path.exists(tsv_path):
            valid_tsv_files.append(tsv_path)
        else:
            logging.warning(f"Input TSV file not found, skipping: {tsv_path}")
            
    if not valid_tsv_files:
        logging.error("No valid input TSV files provided or found. Exiting.")
        sys.exit(1)
    
    logging.info(f"Starting SWAVE_WRITE_VCF for {len(valid_tsv_files)} TSV file(s).")
    
    process_tsv_files_to_vcf(tsv_files=valid_tsv_files, output_vcf_path=options.output_vcf, ref_asm_path=options.ref_fasta)
    
    logging.info("SWAVE_WRITE_VCF completed successfully.")
    sys.exit(0)

if __name__ == "__main__":
    main()
