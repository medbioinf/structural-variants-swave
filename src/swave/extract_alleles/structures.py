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

class Node:
    def __init__(self, node_id, length):
        self.id = node_id
        self.length = length


class Snarl:
    def __init__(self, start_node_id, start_node_orient, end_node_id, end_node_orient, ref_chrom, ref_start, ref_end, reversed_mapping=False):
        self.snarl_id = f"{start_node_orient}{start_node_id}{end_node_orient}{end_node_id}"     # e.g. ">s1>s3"
        
        self.start_node_id = start_node_id
        self.start_node_orient = start_node_orient
        self.end_node_id = end_node_id
        self.end_node_orient = end_node_orient
        
        self.ref_chrom = ref_chrom
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.reversed_mapping = reversed_mapping

        self.ref_asm_path = None
        
        self.path_asm_dict = {}
