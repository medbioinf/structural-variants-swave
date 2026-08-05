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

import logging
import os
import re
import sys

from swave.utils import reverse_complement_seq
from .structures import Node


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def load_nodes_from_fasta(gfa_fasta_path):
    """
    Reads the fasta file of a pangenome graph and extracts the ID, length, and their sequences for each node.
    
    Returns: dict mapping { node_id: Node object },
             dict mapping { node_id: sequence_string }
    """
    logging.info(f"Loading GFA fasta from {gfa_fasta_path}")
    
    nodes_dict = {}
    fasta_index = {}
    
    current_node_id = None
    current_seq_parts = []
    
    with open(gfa_fasta_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                if current_node_id is not None:
                    full_sequence = "".join(current_seq_parts)
                    fasta_index[current_node_id] = full_sequence
                    nodes_dict[current_node_id] = Node(current_node_id, len(full_sequence))
                
                current_node_id = line[1:]  # node ID from header
                current_seq_parts = []
            else:
                current_seq_parts.append(line)
        
        if current_node_id is not None:
            full_sequence = "".join(current_seq_parts)
            fasta_index[current_node_id] = full_sequence
            nodes_dict[current_node_id] = Node(current_node_id, len(full_sequence))
    
    logging.info(f"Successfully loaded {len(nodes_dict)} nodes from fasta")
    return nodes_dict, fasta_index


def extract_and_write_alleles_to_fasta(snarls_dict, fasta_index, output):
    """
    Extracts the allele sequences for each snarl and writes them to a fasta file.
    """
    header_counts = {}
    
    logging.info(f"Writing alleles for {len(snarls_dict)} snarls to {output}")
    
    with open(output, 'w') as output_file:
        for snarl_id, snarl_obj in snarls_dict.items():
            for alt_path, samples in snarl_obj.path_asm_dict.items():
                
                if alt_path == "*":
                    allele_seq = ""
                else:
                    allele_parts = []
                
                    nodes_in_path = re.findall(r'[><]([a-zA-Z0-9]+)', alt_path)
                    orients_in_path = re.findall(r'([><])', alt_path)
                    
                    for i in range(len(nodes_in_path)):
                        node_id = nodes_in_path[i]
                        orient = orients_in_path[i]
                        
                        if node_id not in fasta_index:
                            logging.warning(f"Node {node_id} not found in fasta index, skipping")
                            continue
                        
                        node_seq = fasta_index[node_id]
                        
                        if orient == "<":
                            node_seq = reverse_complement_seq(node_seq)
                        
                        allele_parts.append(node_seq)
                    
                    allele_seq = "".join(allele_parts)
            
                display_seq = allele_seq if allele_seq != "" else "-"
                
                is_reversed_mapping = "true" if snarl_obj.reversed_mapping else "false"
                
                for sample in samples:
                    base_header = f">{sample}|{snarl_id}|{snarl_obj.ref_chrom}:{snarl_obj.ref_start}-{snarl_obj.ref_end}|reversed:{is_reversed_mapping}"
                    
                    # handle rare case of duplicate headers
                    if base_header not in header_counts:
                        header_counts[base_header] = 0
                        fasta_header = f">{base_header}"
                    else:
                        header_counts[base_header] += 1
                        fasta_header = f">{base_header}_{header_counts[base_header]}"   # changes e.g. "reversed:false" to "reversed:false_1"

                    output_file.write(f"{fasta_header}\n")
                    output_file.write(f"{display_seq}\n")
                
    logging.info(f"Successfully extracted snarl alleles.")


def extract_alleles_for_sample(graph_construction_tool, sample_id, gfa_fasta_path, output_dir, options, bed_path=None, vcf_path=None):
    """
    Extracts alleles from the specified graph source (minigraph, pggb, or cactus) and writes them to a fasta file.
    """
    logging.info(f"Extracting alleles for sample '{sample_id}' from {graph_construction_tool} source")
    
    nodes_dict, fasta_index = load_nodes_from_fasta(gfa_fasta_path)

    if graph_construction_tool == "minigraph":
        if bed_path is None:
            raise ValueError("bed_path is required when graph_construction_tool='minigraph'")
        from .minigraph_parsing import parse_minigraph_bed_to_snarls
        results = parse_minigraph_bed_to_snarls(bed_path, sample_id, nodes_dict, options)
        
    elif graph_construction_tool in ("pggb", "cactus"):
        if vcf_path is None:
            raise ValueError("vcf_path is required when graph_construction_tool in ('pggb', 'cactus')")
        from .vg_deconstruct_parsing import parse_vg_deconstructed_vcf_to_snarls
        results = parse_vg_deconstructed_vcf_to_snarls(vcf_path, sample_id, options)
        
    else:
        raise ValueError(f"Unknown graph_construction_tool: {graph_construction_tool}")

    os.makedirs(output_dir, exist_ok=True)
    
    output_paths = []
    for label, snarls_dict in results.items():
        output_path = os.path.join(output_dir, f"{label}_alleles.fa")
        extract_and_write_alleles_to_fasta(snarls_dict, fasta_index, output_path)
        output_paths.append(output_path)

    return output_paths
