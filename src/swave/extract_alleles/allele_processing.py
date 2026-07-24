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
import re
import logging

from src.swave.utils import reverse_complement_seq
from .structures import Node, Snarl


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def load_nodes_from_gfa_fasta(gfa_fasta_path):
    """
    Reads the fasta file of a pangenome graph and extracts the ID, length, and their sequences for each node (S).
    Returns a dictionary { node_id: Node object } and a dictionary { node_id: sequence_string }.
    """
    logging.info(f"Loading node GFA for node metadata: {gfa_fasta_path}")
    
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
                
                current_node_id = line[1:]  # Get node ID from header
                current_seq_parts = []
            else:
                current_seq_parts.append(line)
        
        # Handle the last node after the loop
        if current_node_id is not None:
            full_sequence = "".join(current_seq_parts)
            fasta_index[current_node_id] = full_sequence
            nodes_dict[current_node_id] = Node(current_node_id, len(full_sequence))
    
    logging.info(f"Successfully loaded {len(nodes_dict)} nodes from Fasta.")
    return nodes_dict, fasta_index


def retrieve_reverse_mapping_snarl(contig_mapping_dict, contig_mapping_stats):
    """
    Retrieves reverse mapping snarls based on the contig mapping info and stats.
    """
    
    def get_reverse_path_from_node_list(node_list):
        """
        Given a list of nodes with orientations (e.g. ['>s1', '<s2', '>s3']), returns the reversed path with flipped orientations (e.g. ['<s3', '>s2', '<s1']).
        """
        reversed_path = ""
        
        for i in range(len(node_list) - 1, -1, -1):
            node_with_orient = node_list[i]

            if ">" in node_with_orient:
                reversed_path += node_with_orient.replace(">", "<")
            else:
                reversed_path += node_with_orient.replace("<", ">")

        return reversed_path

    reverse_mapping_snarl_dict = {}

    for cur_contig in contig_mapping_stats:

        cur_contig_stats = contig_mapping_stats[cur_contig]
        cur_contig_list = contig_mapping_dict[cur_contig]

        cur_contig_strand = None
        if (cur_contig_stats['+'] - cur_contig_stats['-']) > 0.8 * cur_contig_stats['+']:
            cur_contig_strand = "+"

        if (cur_contig_stats['-'] - cur_contig_stats['+']) > 0.8 * cur_contig_stats['+']:
            cur_contig_strand = "-"

        if cur_contig_strand is None:
            continue

        previous_reverse_index = None
        previous_reverse_nodes = []
        previous_reverse_cords = []
        previous_reverse_chrom = None

        for sub_index in range(len(cur_contig_list)):
            sub_mapping = cur_contig_list[sub_index]
            sub_mapping_strand = sub_mapping[0]
            sub_start_node_with_orient = sub_mapping[1]
            sub_end_node_orient = sub_mapping[2]
            sub_alt_path = sub_mapping[3]
            ref_chrom = sub_mapping[4]
            ref_start = int(sub_mapping[5])
            ref_end = int(sub_mapping[6])

            if sub_mapping_strand == cur_contig_strand or sub_mapping_strand == ".":
                if previous_reverse_index is not None:
                    reverse_ref_start, reverse_ref_end = min(previous_reverse_cords), max(previous_reverse_cords)
                    reverse_start_node, reverse_end_node = previous_reverse_nodes[0], previous_reverse_nodes[-1]
                    if reverse_ref_end - reverse_ref_start > 5000:
                        reverse_mapping_snarl_dict[f"{reverse_start_node}{reverse_end_node}"] = [
                            cur_contig, get_reverse_path_from_node_list(previous_reverse_nodes),
                            previous_reverse_chrom, reverse_ref_start, reverse_ref_end
                        ]

                previous_reverse_index = None
                previous_reverse_nodes = []
                previous_reverse_cords = []
                previous_reverse_chrom = None

            else:

                if previous_reverse_index is None:
                    previous_reverse_index = sub_index

                previous_reverse_chrom = ref_chrom

                previous_reverse_cords.append(ref_start)
                previous_reverse_cords.append(ref_end)

                # append nodes from snarl and alt path
                if sub_start_node_with_orient not in previous_reverse_nodes:
                    previous_reverse_nodes.append(sub_start_node_with_orient)

                alt_path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', sub_alt_path)
                alt_path_include_nodes_orients = re.findall(r'[><]', sub_alt_path)
                for i in range(len(alt_path_include_nodes)):
                    previous_reverse_nodes.append(f"{alt_path_include_nodes_orients[i]}{alt_path_include_nodes[i]}")

                if sub_end_node_orient not in previous_reverse_nodes:
                    previous_reverse_nodes.append(sub_end_node_orient)

    return reverse_mapping_snarl_dict


def parse_minigraph_bed_to_snarls(bed_path, sample_id, snarls_dict, nodes_dict, options):
    r"""
    Parses the Minigraph --call BED file to extract snarl information and populate the snarls_dict.
    
    The BED file contains the following fields (example):
    chr1    1000    2000    >s1 >s3   <path/node>:\<offset>:\<strand>:\<contig_name>:\<contig_start>:\<contig_end>
    """
    logging.info(f"Parsing Minigraph BED file for {sample_id}: {bed_path}")
    
    with open(bed_path, 'r') as bed_file:
        
        if not os.path.exists(bed_path):
            logging.error(f"BED file not found: {bed_path}")
            sys.exit(-1)
        
        logging.info(f"Preprocessing: Loading allele info from {bed_path} for sample {sample_id}")
        
        # Initialize contig mapping structures for inversion detection
        contig_mapping_dict = {}
        contig_mapping_stats = {}
        prev_seq_source_contig = None

        for line in bed_file:
            parts = line.strip().split('\t')
            
            snarl_ref_chrom, snarl_ref_start, snarl_ref_end = parts[0], int(parts[1]), int(parts[2])
            
            snarl_start_node_with_orient = parts[3]
            snarl_start_node_id = snarl_start_node_with_orient[1:]
            snarl_start_node_orient = snarl_start_node_with_orient[0]
            
            snarl_end_node_with_orient = parts[4]
            snarl_end_node_id = snarl_end_node_with_orient[1:]
            snarl_end_node_orient = snarl_end_node_with_orient[0]
            
            snarl_id = f"{snarl_start_node_orient}{snarl_start_node_id}{snarl_end_node_orient}{snarl_end_node_id}"
            
            if options.spec_snarl and snarl_id != options.spec_snarl:
                continue
            
            alt_path_split = parts[5].split(':')
            alt_path = alt_path_split[0]

            if options.force_reverse:
                # For saving inverted mapping
                if alt_path == ".":
                    seq_source_contig, seq_source_start, seq_source_end, seq_source_strand = prev_seq_source_contig, 0, 0, "."
                else:
                    seq_source_contig, seq_source_start, seq_source_end, seq_source_strand = alt_path_split[3], int(alt_path_split[4]), int(alt_path_split[5]), alt_path_split[2]

                if seq_source_contig not in contig_mapping_dict:
                    contig_mapping_dict[seq_source_contig] = []
                    contig_mapping_stats[seq_source_contig] = {"+": 0, "-": 0, ".": 0}

                contig_mapping_dict[seq_source_contig].append([seq_source_strand, snarl_start_node_with_orient, snarl_end_node_with_orient, alt_path, snarl_ref_chrom, snarl_ref_start, snarl_ref_end])
                contig_mapping_stats[seq_source_contig][seq_source_strand] += (seq_source_end - seq_source_start)

                prev_seq_source_contig = seq_source_contig

            # If assembly has no contig covering the snarl, treat as missing allele
            if alt_path == ".":
                continue
            
            # Option to filter out small snarls/variants below the minimum SV size threshold
            if options.remove_small:
                alt_path_include_nodes = re.findall(r'([a-zA-Z0-9]+)', alt_path)
                alt_path_include_nodes_orients = re.findall(r'([><])', alt_path)

                small_node_indices = []
                for index in range(len(alt_path_include_nodes) - 1, -1, -1):
                    node_id = alt_path_include_nodes[index]
                    
                    if node_id in nodes_dict and nodes_dict[node_id].length < options.min_sv_size:
                        small_node_indices.append(index)
                
                for index in small_node_indices:
                    alt_path_include_nodes.pop(index)
                    alt_path_include_nodes_orients.pop(index)

                alt_path = "".join(f"{alt_path_include_nodes_orients[i]}{alt_path_include_nodes[i]}" for i in range(len(alt_path_include_nodes)))
            
            if alt_path == "":
                alt_path = "*"

            if snarl_id not in snarls_dict:
                snarls_dict[snarl_id] = Snarl(
                    snarl_start_node_id, snarl_start_node_orient, snarl_end_node_id, snarl_end_node_orient,
                    snarl_ref_chrom, snarl_ref_start, snarl_ref_end
                )

            # Add the alt path to the snarls
            if alt_path not in snarls_dict[snarl_id].path_asm_dict:
                snarls_dict[snarl_id].path_asm_dict[alt_path] = []
            snarls_dict[snarl_id].path_asm_dict[alt_path].append(sample_id)
        
        # Handle reversed mappings and add them to the snarls_dict
        if options.force_reverse:
            reverse_mapping_snarl_dict = retrieve_reverse_mapping_snarl(contig_mapping_dict, contig_mapping_stats)

            for snarl_id in reverse_mapping_snarl_dict:
                snarl_info = reverse_mapping_snarl_dict[snarl_id]

                alt_path, ref_chrom, ref_start, ref_end = snarl_info[1], snarl_info[2], int(snarl_info[3]), int(snarl_info[4])

                if alt_path == "":
                    alt_path = "*"

                snarl_id_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', snarl_id)
                snarl_id_include_nodes_orients = re.findall(r'[><]', snarl_id)

                snarl_start_node_id, snarl_start_node_orient = snarl_id_include_nodes[0], snarl_id_include_nodes_orients[0]
                snarl_end_node_id, snarl_end_node_orient = snarl_id_include_nodes[1], snarl_id_include_nodes_orients[1]

                if snarl_id not in snarls_dict:
                    snarls_dict[snarl_id] = Snarl(
                        snarl_start_node_id, snarl_start_node_orient, snarl_end_node_id, snarl_end_node_orient,
                        ref_chrom, ref_start, ref_end, reversed_mapping=True
                    )
                
                # Add the alt path to the snarls
                if alt_path not in snarls_dict[snarl_id].path_asm_dict:
                    snarls_dict[snarl_id].path_asm_dict[alt_path] = []
                snarls_dict[snarl_id].path_asm_dict[alt_path].append(sample_id)


def extract_and_write_alleles_to_fasta(snarls_dict, fasta_index, output):
    """
    Extracts the allele sequences for each snarl and writes them to a FASTA file.
    """
    header_counts = {}
    
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
                            logging.warning(f"Node {node_id} not found in FASTA index. Skipping.")
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
