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

import gzip
import logging
import os
import pickle
import re
import sys

import numpy as np
import pysam
from PIL import Image, ImageDraw, ImageFont

from .structures import Dotplot
from .subsnarl_building import resolve_large_snarl_path, build_subsnarl_sequences
from swave.utils import calculate_stride_size


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def process_sample_alleles_to_dotplots(alt_fasta_path, ref_fasta_path, pangenome_gfa_fasta_path, options=None):
    """
    Reads the extracted sample allele fasta (ALT), resolves each snarl into
    sub-snarls where the reference and alternative paths diverge, and
    generates dotplots for each sub-snarl. Stores the dotplots for each
    sub-snarl into a dictionary and saves it as a pickle file.
    """
    ref_file = pysam.FastaFile(ref_fasta_path)
    gfa_fasta_file = pysam.FastaFile(pangenome_gfa_fasta_path)

    with open(alt_fasta_path, 'r') as alt_file:
        total_snarl_count = sum(1 for line in alt_file if line.startswith('>'))

    snarl_count = 0
    current_header = None
    snarl_dotplot_dict = {}

    with open(alt_fasta_path, 'r') as alt_file:
        for line in alt_file:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                current_header = line[1:]
            else:
                snarl_count += 1

                header_parts = current_header.split('|')
                snarl_id = header_parts[1]
                coords = header_parts[2]
                is_reversed_mapping = header_parts[3].split(':')[1] == "true"
                alt_path = header_parts[4].split(':', 1)[1]
                ref_path = header_parts[5].split(':', 1)[1]

                chrom, pos_range = coords.split(':')
                snarl_ref_start, snarl_ref_end = map(int, pos_range.split('-'))

                process_and_plot_snarl(
                    snarl_id=snarl_id,
                    chrom=chrom,
                    snarl_ref_start=snarl_ref_start,
                    snarl_ref_end=snarl_ref_end,
                    alt_path=alt_path,
                    ref_path=ref_path,
                    gfa_fasta_file=gfa_fasta_file,
                    ref_file=ref_file,
                    snarl_dotplot_dict=snarl_dotplot_dict,
                    is_reversed_mapping=is_reversed_mapping,
                    snarl_count=snarl_count,
                    total_snarl_count=total_snarl_count,
                    options=options
                )

    ref_file.close()
    gfa_fasta_file.close()

    output_pickle_path = f"{options.pkl_out_prefix}_dotplots.pkl.gz"

    with open(output_pickle_path, 'wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb') as gz:
            pickle.dump(snarl_dotplot_dict, gz)

    logging.info(f"Saved {len(snarl_dotplot_dict)} snarl dotplot bundles")


def encode_path_for_filename(path_str):
    """
    Encodes a node path for safe use in filenames: '>' becomes 'F' (forward),
    '<' becomes 'R' (reverse), '' / '*' becomes 'none'.
    """
    if path_str in ("", "*"):
        return "none"
    return path_str.replace(">", "F").replace("<", "R")


def process_and_plot_snarl(snarl_id, chrom, snarl_ref_start, snarl_ref_end, alt_path, ref_path,
                            gfa_fasta_file, ref_file, snarl_dotplot_dict, is_reversed_mapping,
                            snarl_count, total_snarl_count, options):
    """
    Resolves one snarl into its sub-snarls (skipping regions where the
    reference and alternative paths do not diverge) and generates a dotplot
    bundle for each, storing the results in snarl_dotplot_dict.
    """

    if snarl_count % 1000 == 0:
        logging.info(f"Generating {snarl_count}th snarl. In total {total_snarl_count} snarls")

    if options.spec_path is not None and options.spec_path not in snarl_id:
        return

    if is_reversed_mapping:
        if (snarl_ref_end - snarl_ref_start) > options.max_sv_size:
            logging.info(f"Skipping reversed mapping snarl {snarl_id} due to size > {options.max_sv_size} bp")
            return

        nodes_with_orients = re.findall(r'([><][a-zA-Z0-9]+)', snarl_id)
        sub_snarls = [[nodes_with_orients[0], nodes_with_orients[-1], snarl_ref_start, snarl_ref_end,
                       ref_path, alt_path]]
    else:
        if alt_path == ref_path:
            return

        nodes_with_orients = re.findall(r'([><][a-zA-Z0-9]+)', snarl_id)
        snarl_start_node_with_orient = nodes_with_orients[0]
        snarl_end_node_with_orient = nodes_with_orients[-1]

        sub_snarls = resolve_large_snarl_path(
            ref_path, alt_path, snarl_start_node_with_orient, snarl_end_node_with_orient,
            snarl_ref_start, snarl_ref_end, gfa_fasta_file.get_reference_length
        )

        if len(sub_snarls) == 0:
            return

    for sub_snarl in sub_snarls:
        sub_start_node_with_orient, sub_end_node_with_orient, sub_ref_start, sub_ref_end, sub_ref_path, sub_alt_path = sub_snarl

        sub_snarl_key = "{}-{}".format(sub_ref_path, sub_alt_path)

        if options.spec_path is not None and options.spec_path not in sub_snarl_key:
            continue

        sub_ref_seq, sub_alt_seq, dotplot_ref_start, dotplot_ref_end = build_subsnarl_sequences(
            gfa_fasta_file, ref_file, chrom, sub_snarl
        )

        dotplot_id = "{}::{}|{}|{}|{}|{}|{}|rev_{}".format(
            snarl_id, sub_snarl_key, sub_ref_start, sub_ref_end, chrom,
            dotplot_ref_start, dotplot_ref_end, str(is_reversed_mapping).lower()
        )

        dotplot_filename = "{}_{}-{}__{}__{}_vs_{}__pad{}-{}__rev{}".format(
            chrom, sub_ref_start, sub_ref_end,
            encode_path_for_filename(snarl_id),
            encode_path_for_filename(sub_ref_path), encode_path_for_filename(sub_alt_path),
            dotplot_ref_start, dotplot_ref_end,
            "1" if is_reversed_mapping else "0"
        )
        dotplot_output_prefix = os.path.join(options.img_out_prefix, dotplot_filename)

        dotplot_stride_size = calculate_stride_size(sub_ref_seq, sub_alt_seq)

        dotplot_objects_bundle = generate_dotplots(
            ref_seq=sub_ref_seq,
            alt_seq=sub_alt_seq,
            dotplot_stride_size=dotplot_stride_size,
            dotplot_output_prefix=dotplot_output_prefix,
            options=options
        )

        snarl_dotplot_dict[dotplot_id] = dotplot_objects_bundle


def save_combined_dotplot_grid(bundle, output_path):
    """
    Renders the four dotplot matrices in a bundle (ref2ref, ref2alt,
    alt2alt, and the reverse-complement ref2alt) as one combined 2x2
    grid image.
    """
    m_ref2ref = bundle["x2x_ref2ref"].matrix
    m_ref2alt = bundle["x2y_ref2alt"].matrix
    m_alt2alt = bundle["x2x_alt2alt"].matrix
    m_ref2alt_rev = bundle["x2y_ref2alt"].matrix_rev
    
    if hasattr(m_ref2alt, "size") and m_ref2alt.size == 0:
        logging.warning(f"Skipping dotplot pngs generation for {output_path}: ref2alt matrix is empty.")
        return
    elif hasattr(m_ref2alt, "shape") and (m_ref2alt.shape[0] == 0 or m_ref2alt.shape[1] == 0):
        logging.warning(f"Skipping dotplot pngs generation for {output_path}: ref2alt matrix has zero dimension.")
        return
    
    h, w = m_ref2alt.shape
    max_side = max(h, w)
    padding = max(40, int(max_side * 0.12))
    quad_size = max_side + (2 * padding)
    
    gray_bg_val = 245
    grid_array = np.full((quad_size * 2, quad_size * 2), gray_bg_val, dtype=np.uint8)
    
    quadrants = [
        (m_ref2ref, 0, 0, "ref2ref"),
        (m_ref2alt, 0, quad_size, "ref2alt"),
        (m_alt2alt, quad_size, 0, "alt2alt"),
        (m_ref2alt_rev, quad_size, quad_size, "ref2alt_rev")
    ]
    
    for m, q_y, q_x, label in quadrants:
        if m is None:
            continue
        
        m_h, m_w = getattr(m, "shape", (0, 0))
        if m_h == 0 or m_w == 0:
            continue
        
        m_h, m_w = m.shape
        
        offset_y = q_y + padding + (max_side - m_h) // 2
        offset_x = q_x + padding + (max_side - m_w) // 2
        
        grid_array[offset_y:offset_y + m_h, offset_x:offset_x + m_w] = 255
        
        has_points = (hasattr(m, "nnz") and m.nnz > 0) or (hasattr(m, "size") and m.size > 0 and np.any(m))
        if has_points:
            if hasattr(m, "tocoo"):
                coo = m.tocoo()
                grid_array[offset_y + coo.row, offset_x + coo.col] = 0
            else:
                grid_array[offset_y:offset_y + m_h, offset_x:offset_x + m_w] = np.where(m > 0, 0, 255)
    
    img = Image.fromarray(grid_array, mode='L')
    draw = ImageDraw.Draw(img)
    
    total_size_px = quad_size * 2
    font_size = max(14, int(total_size_px / 50))

    try:
        font = ImageFont.load_default(size=font_size)
    except Exception:
        font = None

    for m, q_y, q_x, label in quadrants:
        if m is None:
            continue
        
        m_h, m_w = getattr(m, "shape", (0, 0))
        if m_h == 0 or m_w == 0:
            continue
        
        offset_y = q_y + padding + (max_side - m_h) // 2
        offset_x = q_x + padding + (max_side - m_w) // 2
        
        text = f"{label} ({m_w}x{m_h})"
        
        if hasattr(draw, "textlength") and font is not None:
            text_width = draw.textlength(text, font=font)
        else:
            text_width = len(text) * (font_size * 0.5)

        text_x = offset_x + (m_w // 2) - (text_width // 2)
        text_y = offset_y - int(padding * 0.65)
        
        draw.text((text_x, text_y), text, fill=0, font=font)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    img.save(output_path, format="PNG")


def generate_dotplots(ref_seq, alt_seq, dotplot_stride_size, dotplot_output_prefix, options):
    """
    Generates the dotplot objects for the given reference and alternative sequence and saves matrix PNG visualizations if specified.
    """
    if options.save_dotplot_images:
        os.makedirs(os.path.dirname(dotplot_output_prefix), exist_ok=True)

    x2x_ref2ref = Dotplot(ref_seq, ref_seq, options.kmer_size, out_prefix=f"{dotplot_output_prefix}.ref2ref", stride_size=dotplot_stride_size,
                          skip_forward=options.skip_forward, skip_reverse=options.skip_reverse)
    
    x2y_ref2alt = Dotplot(ref_seq, alt_seq, options.kmer_size, out_prefix=f"{dotplot_output_prefix}.ref2alt", stride_size=dotplot_stride_size,
                           given_x_kmer_index=x2x_ref2ref.get_seq_x_kmer_index(), skip_forward=options.skip_forward, skip_reverse=options.skip_reverse)
    
    x2x_alt2alt = Dotplot(alt_seq, alt_seq, options.kmer_size, out_prefix=f"{dotplot_output_prefix}.alt2alt", stride_size=dotplot_stride_size,
                          skip_forward=options.skip_forward, skip_reverse=options.skip_reverse)

    if options.save_dotplot_images:
        save_combined_dotplot_grid({
            "x2x_ref2ref": x2x_ref2ref,
            "x2y_ref2alt": x2y_ref2alt,
            "x2x_alt2alt": x2x_alt2alt
        }, f"{dotplot_output_prefix}_dotplots.png")

    return {
        "x2x_ref2ref": x2x_ref2ref,
        "x2y_ref2alt": x2y_ref2alt,
        "x2x_alt2alt": x2x_alt2alt,
        "stride_size": dotplot_stride_size
    }
