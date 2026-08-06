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

import re

from swave.utils import reverse_complement_seq


def resolve_large_snarl_path(ref_path, alt_path, start_node_with_orient, end_node_with_orient,
                              snarl_ref_start, snarl_ref_end, node_length_fn, min_anchor_length=1000):
    """
    Given the reference and alternative interior node paths through a snarl,
    finds nodes shared between them ("anchors") and splits the snarl into
    sub-snarls at those anchors, restricted to the regions where the
    reference and alternative paths actually diverge.

    Returns: list of [sub_start_node_with_orient, sub_end_node_with_orient, sub_start_pos, sub_end_pos, sub_ref_path, sub_alt_path]
    """
    # check nodes in the ref path, determine if they exist in path
    ref_path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', ref_path)
    ref_path_include_nodes_orients = re.findall(r'[><]', ref_path)
    ref_path_split = [start_node_with_orient] + ["{}{}".format(ref_path_include_nodes_orients[i], ref_path_include_nodes[i]) for i in range(len(ref_path_include_nodes))] + [end_node_with_orient]

    alt_path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', alt_path)
    alt_path_include_nodes_orients = re.findall(r'[><]', alt_path)
    alt_path_split = [start_node_with_orient] + ["{}{}".format(alt_path_include_nodes_orients[i], alt_path_include_nodes[i]) for i in range(len(alt_path_include_nodes))] + [end_node_with_orient]

    # find the start pos for each ref node
    ref_pos_pointer = snarl_ref_start
    ref_path_include_nodes_pos = []
    for i in range(len(ref_path_include_nodes)):
        ref_path_include_nodes_pos.append(ref_pos_pointer)
        ref_pos_pointer += node_length_fn(ref_path_include_nodes[i])

    ref_path_split_pos = [snarl_ref_start - node_length_fn(start_node_with_orient[1:])] + ref_path_include_nodes_pos + [snarl_ref_end]

    # find shared nodes
    shared_nodes_with_orients = []
    for node_with_orient in ref_path_split:
        if node_with_orient in alt_path_split:
            if not (node_with_orient in [start_node_with_orient, end_node_with_orient]) and node_length_fn(node_with_orient[1:]) < min_anchor_length:
                continue
            shared_nodes_with_orients.append(node_with_orient)

    # generate sub snarls
    sub_snarls = []
    for i in range(len(shared_nodes_with_orients) - 1):

        sub_start_node_with_orient = shared_nodes_with_orients[i]
        sub_end_node_with_orient = shared_nodes_with_orients[i + 1]

        sub_start_pos = ref_path_split_pos[ref_path_split.index(sub_start_node_with_orient) + 1]   # +1: the start pos is end pos of start node, which is also start pos of next node
        sub_end_pos = ref_path_split_pos[ref_path_split.index(sub_end_node_with_orient)]

        sub_ref_path_split = ref_path_split[ref_path_split.index(sub_start_node_with_orient): ref_path_split.index(sub_end_node_with_orient) + 1]
        sub_alt_path_split = alt_path_split[alt_path_split.index(sub_start_node_with_orient): alt_path_split.index(sub_end_node_with_orient) + 1]

        if sub_ref_path_split == sub_alt_path_split:
            continue

        sub_snarls.append([
            sub_start_node_with_orient, sub_end_node_with_orient, sub_start_pos, sub_end_pos,
            "".join(sub_ref_path_split[1:-1]), "".join(sub_alt_path_split[1:-1])   # [1:-1]: remove start/end node from path
        ])

    return sub_snarls

def fetch_node_seq(gfa_fasta_file, node_with_orient):
    """
    Fetches a single graph node's sequence, applying reverse-complement if
    the node is traversed in reverse orientation.
    """
    node_id = node_with_orient[1:]
    node_orient = node_with_orient[0]
    seq = gfa_fasta_file.fetch(node_id).upper()
    return seq if node_orient == ">" else reverse_complement_seq(seq)


def fetch_path_seq(gfa_fasta_file, path_str):
    """
    Fetches and concatenates the sequence for an interior node path
    (e.g. ">s55>s86"). Returns an empty string if path_str is empty.
    """
    if path_str == "":
        return ""
    nodes_with_orients = re.findall(r'([><][a-zA-Z0-9]+)', path_str)
    return "".join(fetch_node_seq(gfa_fasta_file, n) for n in nodes_with_orients)


def build_subsnarl_sequences(gfa_fasta_file, ref_file, chrom, sub_snarl):
    """
    Builds padded reference and alternative sequences for one sub-snarl,
    using graph node sequences for both sides (including the padding/anchor).
    Falls back to the linear reference genome when the sub-snarl's reference
    interior path is empty, or when an anchor node's own sequence is shorter
    than the desired padding length.

    sub_snarl: [sub_start_node_with_orient, sub_end_node_with_orient, sub_start_pos, sub_end_pos, sub_ref_path, sub_alt_path]

    Returns: sub_ref_seq, sub_alt_seq, sub_dotplot_ref_start, sub_dotplot_ref_end
    """
    sub_start_node_with_orient, sub_end_node_with_orient, sub_start_pos, sub_end_pos, sub_ref_path, sub_alt_path = sub_snarl

    interior_ref_seq = fetch_path_seq(gfa_fasta_file, sub_ref_path)
    interior_alt_seq = fetch_path_seq(gfa_fasta_file, sub_alt_path)

    chrom_len = ref_file.get_reference_length(chrom)

    if sub_ref_path == "":
        interior_ref_seq = ref_file.fetch(chrom, sub_start_pos, sub_end_pos).upper()

    extend_len = min(10000, 2 * max([len(interior_ref_seq), len(interior_alt_seq), (sub_end_pos - sub_start_pos)]))

    start_anchor_seq = fetch_node_seq(gfa_fasta_file, sub_start_node_with_orient)[-extend_len:]
    if len(start_anchor_seq) < extend_len:
        start_anchor_seq = ref_file.fetch(chrom, max(sub_start_pos - extend_len, 0), sub_start_pos).replace("N", "").upper()

    end_anchor_seq = fetch_node_seq(gfa_fasta_file, sub_end_node_with_orient)[:extend_len]
    if len(end_anchor_seq) < extend_len:
        end_anchor_seq = ref_file.fetch(chrom, sub_end_pos, min(chrom_len, sub_end_pos + extend_len)).replace("N", "").upper()

    sub_ref_seq = start_anchor_seq + interior_ref_seq + end_anchor_seq
    sub_alt_seq = start_anchor_seq + interior_alt_seq + end_anchor_seq

    sub_dotplot_ref_start = sub_start_pos - len(start_anchor_seq)
    sub_dotplot_ref_end = sub_end_pos + len(end_anchor_seq)

    return sub_ref_seq, sub_alt_seq, sub_dotplot_ref_start, sub_dotplot_ref_end
