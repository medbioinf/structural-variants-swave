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

from swave.extract_alleles import extract_alleles_for_sample
from swave.version import __version__


def main():
    parser = argparse.ArgumentParser(description="Extracts structural variant alleles for a sample or reference from a pangenome graph, using either a minigraph BED file or a vg-deconstruct VCF as the allele-path source.")
    
    parser.add_argument("--version", action="version", version=f"{__version__}")
    
    parser.add_argument("--graph_construction_tool", required=True, choices=["minigraph", "pggb", "cactus"], help="Pangenome graph construction tool that produced the input allele-path file.")
    parser.add_argument("--gfa_fasta", required=True, help="Path to the gfatools gfa2fa node FASTA file.")
    parser.add_argument("--sample_id", required=True, help="Name/ID of the sample being processed.")
    parser.add_argument("--output_dir", default=".", help="Directory to write the extracted allele FASTA file(s) into (default: current directory).")
    parser.add_argument("--bed", default=None, help="Path to the sample minigraph --call BED file. Required when --graph_construction_tool=minigraph.")
    parser.add_argument("--vcf", default=None, help="Path to the sample or reference vg-deconstruct VCF file. Required when --graph_construction_tool=pggb or cactus.")

    # optional parameters
    parser.add_argument("--spec_snarl", default=None, help="Specific snarl ID to process (e.g. '>s1>s3') for debugging/analysis. If not provided, all snarls will be extracted.")
    parser.add_argument("--force_reverse", action="store_true", help="Enable original Swave inversion detection and rescue logic for reversed contigs. Only applies when --graph_construction_tool=minigraph.")
    parser.add_argument("--remove_small", action="store_true", help="Filter out small snarls/variants below the minimum SV size threshold. Only applies when --graph_construction_tool=minigraph.")
    parser.add_argument("--min_sv_size", type=int, default=50,help="Minimum size (in base pairs) for a variant to be considered a structural variant (default: 50).")
    
    options = parser.parse_args()

    if options.graph_construction_tool == "minigraph" and options.bed is None:
        parser.error("--bed is required when --graph_construction_tool=minigraph")
    if options.graph_construction_tool in ("pggb", "cactus") and options.vcf is None:
        parser.error("--vcf is required when --graph_construction_tool=pggb or cactus")

    extract_alleles_for_sample(
        graph_construction_tool=options.graph_construction_tool,
        sample_id=options.sample_id,
        gfa_fasta_path=options.gfa_fasta,
        output_dir=options.output_dir,
        options=options,
        bed_path=options.bed,
        vcf_path=options.vcf,
    )

    sys.exit(0)


if __name__ == "__main__":
    main()
