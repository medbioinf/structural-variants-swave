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
import re
import sys

import pysam

from .structures import Snarl


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def remove_snarl_nodes_from_path(raw_path):
    raw_path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', raw_path)
    raw_path_include_nodes_orients = re.findall(r'[><]', raw_path)

    new_path = "".join("{}{}".format(raw_path_include_nodes_orients[i], raw_path_include_nodes[i]) for i in range(1, len(raw_path_include_nodes) - 1))
    if new_path == "":
        new_path = "*"

    return new_path


def parse_vg_deconstructed_vcf_to_snarls(vcf_path, sample_id, options):
    """
    Parses a vg-deconstruct VCF (either a single-sample split or a reference
    VCF with no sample columns) and returns one Snarl dict per haplotype
    (or a single dict for the reference).

    Returns: dict mapping { output_label: snarls_dict }
             e.g. {"HG002_hap1": {...}, "HG002_hap2": {...}}
             or {"t2t_chm13v2_0": {...}} for the reference
             or {"HG002_hap0": {...}} for a haploid sample
    """
    logging.info(f"Parsing vg-deconstruct VCF for '{sample_id}': {vcf_path}")

    results = {}
    seen_snarl_ids = set()

    with pysam.VariantFile(vcf_path) as vcf_file:
        sample_order = list(vcf_file.header.samples)
        is_ref_only = len(sample_order) == 0

        if is_ref_only:
            results[sample_id] = {}

        for record in vcf_file:

            if not (len(record.ref) >= options.min_sv_size or max(len(alt) for alt in record.alts) >= options.min_sv_size):
                continue

            snarl_id = record.id
            snarl_ref_chrom = record.contig
            snarl_ref_start = record.start + 1
            snarl_ref_end = record.stop + 1

            if options.spec_snarl is not None and snarl_id != options.spec_snarl:
                continue

            if snarl_id in seen_snarl_ids:
                continue
            seen_snarl_ids.add(snarl_id)

            snarl_id_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', snarl_id)
            snarl_id_include_nodes_orients = re.findall(r'[><]', snarl_id)

            snarl_start_node_id, snarl_start_node_orient = snarl_id_include_nodes[0], snarl_id_include_nodes_orients[0]
            snarl_end_node_id, snarl_end_node_orient = snarl_id_include_nodes[1], snarl_id_include_nodes_orients[1]

            if not (snarl_start_node_orient == ">" and snarl_end_node_orient == ">"):
                continue

            # AT includes the snarl's start/end node in every path; strip them
            all_paths = [remove_snarl_nodes_from_path(path) for path in record.info["AT"]]
            ref_asm_path = all_paths[0]

            if is_ref_only:
                ref_snarls_dict = results[sample_id]
                if snarl_id not in ref_snarls_dict:
                    ref_snarls_dict[snarl_id] = Snarl(
                        snarl_start_node_id, snarl_start_node_orient, snarl_end_node_id, snarl_end_node_orient,
                        snarl_ref_chrom, snarl_ref_start, snarl_ref_end
                    )
                snarl_obj = ref_snarls_dict[snarl_id]
                if ref_asm_path not in snarl_obj.path_asm_dict:
                    snarl_obj.path_asm_dict[ref_asm_path] = []
                snarl_obj.path_asm_dict[ref_asm_path].append(sample_id)
                snarl_obj.ref_asm_path = ref_asm_path
                continue

            # sample case: one GT column per sample present in the VCF
            sample_gts = str(record).strip().split("\t")[9:]
            for sample_index in range(len(sample_gts)):
                sample_gt = sample_gts[sample_index]
                sample_name = sample_order[sample_index]

                hap_gts = sample_gt.replace("/", "|").split("|")
                ploidy = len(hap_gts)

                for hap_index in range(ploidy):
                    hap_gt = hap_gts[hap_index]

                    if hap_gt == ".":
                        continue

                    alt_path = all_paths[int(hap_gt)]
                    hap_label = f"{sample_name}_hap0" if ploidy == 1 else f"{sample_name}_hap{hap_index + 1}"

                    if hap_label not in results:
                        results[hap_label] = {}
                    hap_snarls_dict = results[hap_label]

                    if snarl_id not in hap_snarls_dict:
                        hap_snarls_dict[snarl_id] = Snarl(
                            snarl_start_node_id, snarl_start_node_orient, snarl_end_node_id, snarl_end_node_orient,
                            snarl_ref_chrom, snarl_ref_start, snarl_ref_end
                        )
                    snarl_obj = hap_snarls_dict[snarl_id]

                    if alt_path not in snarl_obj.path_asm_dict:
                        snarl_obj.path_asm_dict[alt_path] = []
                    snarl_obj.path_asm_dict[alt_path].append(hap_label)
                    
                    snarl_obj.ref_asm_path = ref_asm_path

    logging.info(f"Parsed {len(results)} haplotype/reference output(s) from {vcf_path}")
    return results
