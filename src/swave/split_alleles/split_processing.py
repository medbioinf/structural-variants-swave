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


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def split_and_interleave_fasta(fasta_path, prefix, seq_per_split):
    """   
    Reads a fasta file, sorts the alleles in descending order by length, and distributes
    them in a round-robin manner across (seq_per_fasta // seq_per_split) files.
    """
    logging.info(f"Reading fasta file: {fasta_path}")
    
    records = []
    current_header = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    seq_str = "".join(current_seq)
                    records.append((current_header, seq_str, len(seq_str)))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            seq_str = "".join(current_seq)
            records.append((current_header, seq_str, len(seq_str)))
    
    if not records:
        logging.warning("No sequences found in the fasta file")
        empty_filename = f"{prefix}.split_01.fa"
        with open(empty_filename, "w") as f:
            pass
        return
    
    num_splits = len(records) // seq_per_split
    
    if len(records) % seq_per_split != 0:
        num_splits += 1

    logging.info(f"{len(records)} sequences loaded from fasta file")
    records.sort(key=lambda x: x[2], reverse=True)  # sort by length in descending order
    
    output_files = []
    for i in range(1, num_splits + 1):
        filename = f"{prefix}.split_{i:02d}.fa"
        output_files.append(open(filename, "w"))

    for idx, (header, seq, _length) in enumerate(records):
        bucket = idx % num_splits
        output_files[bucket].write(f"{header}\n{seq}\n")

    for fh in output_files:
        fh.close()

    logging.info(f"Successfully created {num_splits} balanced fasta files")
