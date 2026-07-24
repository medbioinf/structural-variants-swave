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
import gzip
import pickle
import logging
import numpy as np
from matplotlib import pyplot as plt

from .diag_finder import find_and_denoise_diags, boost_diags
from .diag_extension import (
    linear_level_dotplot_extension,
    kmer10_level_dotplot_extension,
    base_level_dotplot_extension,
    rotate_diags_to_alt2ref
)


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def process_dotplot_bundles_to_projections(dotplot_bundles_pickle_path, projections_output_prefix, options=None):
    """
    Loads the precomputed dotplot bundles from a pickle file, processes each bundle to generate projections,
    and saves the resulting projections into a new pickle file.
    
    Args:
        dotplot_bundles_pickle_path (str): Path to the pickle file containing dotplot bundles
        projections_output_prefix (str): Prefix for the output pickle file containing projections
        options: Additional options for processing.
    """
    logging.info(f"Loading dotplot bundles from {dotplot_bundles_pickle_path}...")
    
    with gzip.open(dotplot_bundles_pickle_path, 'rb') as pickle_file:
        snarl_dotplot_dict = pickle.load(pickle_file)
    
    logging.info(f"Processing {len(snarl_dotplot_dict)} dotplot bundles to generate projections...")
    
    projections_dict = {}
    processed_count = 0
    
    for dotplot_id, dotplot_bundle in snarl_dotplot_dict.items():
        processed_count += 1
        
        if processed_count % 1000 == 0 or processed_count == len(snarl_dotplot_dict):
            logging.info(f"Processing dotplot bundle {processed_count}/{len(snarl_dotplot_dict)}: {dotplot_id}...")
            
        projection_matrices = generate_projections(dotplot_bundle, options)
        
        if projection_matrices["ref2alt"] == "Bad":
            logging.warning(f"Bad dotplot detected for snarl: {dotplot_id}. Skipping.")
        
        projections_dict[dotplot_id] = projection_matrices
    
    with gzip.open(projections_output_prefix + "_projections.pkl.gz", 'wb') as output_file:
        pickle.dump(projections_dict, output_file, protocol=pickle.HIGHEST_PROTOCOL)
    
    logging.info(f"Successfully processed and saved {len(projections_dict)} projections for {len(snarl_dotplot_dict)} snarls.")    
    

def generate_projections(dotplot_bundle, options):
    """
    Takes a precomputed dotplot bundle of a snarl, performs the diag-finder/extension logic,
    and returns the final projection matrices for ref2alt and alt2ref.
    """
    x2x_ref2ref = dotplot_bundle["x2x_ref2ref"]
    x2y_ref2alt = dotplot_bundle["x2y_ref2alt"]
    x2x_alt2alt = dotplot_bundle["x2x_alt2alt"]
    dotplot_stride_size = dotplot_bundle["stride_size"]

    dotplot_projection_dict = {"ref2alt": [], "alt2ref": []}
    
    unique_diags = find_and_denoise_diags(x2x_ref2ref.matrix, x2y_ref2alt.matrix, x2y_ref2alt.matrix_rev, x2x_ref2ref.stride_size)
    
    if len(unique_diags) == 0:
        return {"ref2alt": "Bad", "alt2ref": "Bad"}
    
    # --- ref2alt ---
    x2x_dotplot_project_x, augment_coeff = x2x_ref2ref.get_project_x(augment=True)
    x2x_dotplot_project_x_rev = x2x_ref2ref.get_project_x_rev(baseline=augment_coeff)
    
    x2y_dotplot_project_x, _ = x2y_ref2alt.get_project_x(augment=False)
    x2y_dotplot_project_x_rev = x2y_ref2alt.get_project_x_rev(baseline=augment_coeff)
    
    # two rounds of linear extension
    unique_diags = linear_level_dotplot_extension(unique_diags, x2y_ref2alt.matrix)
    unique_diags = linear_level_dotplot_extension(unique_diags, x2y_ref2alt.matrix)
    
    # kmer10 and base level extension
    unique_diags = kmer10_level_dotplot_extension(unique_diags, x2y_ref2alt.seq_x, x2y_ref2alt.seq_y, x2y_ref2alt.matrix, dotplot_stride_size)
    unique_diags = base_level_dotplot_extension(
        unique_diags, x2y_ref2alt.seq_x, x2y_ref2alt.seq_y, x2y_ref2alt.matrix, x2y_dotplot_project_x,
        x2y_dotplot_project_x_rev, options.kmer_size, x2x_ref2ref.stride_size
    )
    
    x2y_matrix_len_x, x2y_matrix_len_y = x2y_ref2alt.matrix.shape[1], x2y_ref2alt.matrix.shape[0]
    
    # boost diags and segment projections into matrices
    boost_diags("ref2alt", unique_diags, x2y_matrix_len_x, x2y_matrix_len_y, x2y_dotplot_project_x, x2y_dotplot_project_x_rev, augment_coeff)
    
    if options.save_projections_images:
        output_synthesized_dotplot_project(
            dotplot_bundle["x2y_ref2alt"].out_prefix, 
            x2x_dotplot_project_x, 
            x2y_dotplot_project_x, 
            x2x_dotplot_project_x_rev, 
            x2y_dotplot_project_x_rev
        )
    
    projection_matrix, bad_flag = segment_projections_into_matrix(
        x2x_dotplot_project_x, x2y_dotplot_project_x, x2x_dotplot_project_x_rev, x2y_dotplot_project_x_rev, None, augment_coeff
    )
    
    if bad_flag:
        return {"ref2alt": "Bad", "alt2ref": "Bad"}
    
    dotplot_projection_dict["ref2alt"].extend(projection_matrix)
    
    
    # --- alt2ref ---
    x2x_dotplot_project_x, augment_coeff = x2x_alt2alt.get_project_x(augment=True)
    x2x_dotplot_project_x_rev = x2x_alt2alt.get_project_x_rev(baseline=augment_coeff)

    x2y_dotplot_project_x, _ = x2y_ref2alt.get_project_y(augment=False)
    x2y_dotplot_project_x_rev = x2y_ref2alt.get_project_y_rev(baseline=augment_coeff)

    unique_diags = rotate_diags_to_alt2ref(unique_diags)
    x2y_matrix_len_x, x2y_matrix_len_y = x2y_ref2alt.matrix.shape[0], x2y_ref2alt.matrix.shape[1]

    boost_diags("alt2ref", unique_diags, x2y_matrix_len_x, x2y_matrix_len_y, x2y_dotplot_project_x, x2y_dotplot_project_x_rev, augment_coeff)
    projection_matrix, bad_flag = segment_projections_into_matrix(x2x_dotplot_project_x, x2y_dotplot_project_x, x2x_dotplot_project_x_rev, x2y_dotplot_project_x_rev, None, augment_coeff)

    if bad_flag:
        return {"ref2alt": "Bad", "alt2ref": "Bad"}

    dotplot_projection_dict["alt2ref"].extend(projection_matrix)
    
    return dotplot_projection_dict


def segment_projections_into_matrix(x2x_dotplot_project_x, x2y_dotplot_project_x, x2x_dotplot_project_x_rev, x2y_dotplot_project_x_rev, ideal_x2y_project_type, augment_coeff):
    projection_matrix = []  # format: [[segment_start, segment_end, segment_x2x, segment_x2y, segment_x2y_rev, segment_label], [], [], ....]

    dotplot_x2y_project_subtract = x2x_dotplot_project_x - x2y_dotplot_project_x
    dotplot_x2y_project_rev_subtract = x2x_dotplot_project_x_rev - x2y_dotplot_project_x_rev

    # segment using iteration
    previous_vals = [augment_coeff, dotplot_x2y_project_subtract[0], dotplot_x2y_project_rev_subtract[0]]    # [x2x_val, x2y_val, x2y_val_rev]
    previous_pointer = 0

    for i in range(1, len(x2x_dotplot_project_x) + 1):
        if i == len(x2x_dotplot_project_x):
            new_vals = [-1, -1, -1]     # meet the end
        else:
            new_vals = [augment_coeff, dotplot_x2y_project_subtract[i], dotplot_x2y_project_rev_subtract[i]]

        # find a segment
        if new_vals != previous_vals or i == len(x2x_dotplot_project_x):
            segment_start = previous_pointer
            segment_end = i

            projection_matrix.append([segment_start, segment_end, previous_vals[0], previous_vals[1], previous_vals[2], "None"])

            previous_vals = new_vals
            previous_pointer = i

    # proportion of projection to the raw dotplot (this is also the compression number)
    if np.shape(projection_matrix)[0] / np.shape(x2x_dotplot_project_x)[0] > 0.8:
        projection_matrix = [[0, len(x2x_dotplot_project_x), augment_coeff, 0, 0, 0, "None"]]
        return projection_matrix, True      # bad flag as True

    return projection_matrix, False


# TODO: adjust
def output_synthesized_dotplot_project(dotplot_output_prefix, x2x_dotplot_project_x, x2y_dotplot_project_x, x2x_dotplot_project_x_rev, x2y_dotplot_project_x_rev):
    dirname = os.path.dirname(dotplot_output_prefix)
    if dirname:
        os.makedirs(dirname, exist_ok=True)
    
    plt.figure(figsize=(10, 10))
    ax1 = plt.subplot(611)
    ax1.plot(x2x_dotplot_project_x)
    ax1.set_title("x2x Project X")
    ax1.set_xticks([])

    ax2 = plt.subplot(612)
    ax2.plot(x2y_dotplot_project_x)
    ax2.set_title("x2y Project X")
    ax2.set_xticks([])

    ax3 = plt.subplot(613)
    ax3.plot(x2x_dotplot_project_x - x2y_dotplot_project_x)
    ax3.set_title("Difference (x2x - x2y)")
    ax3.set_xticks([])

    ax4 = plt.subplot(614)
    ax4.plot(x2x_dotplot_project_x_rev)
    ax4.set_title("x2x Project X (Reverse)")
    ax4.set_xticks([])

    ax5 = plt.subplot(615)
    ax5.plot(x2y_dotplot_project_x_rev)
    ax5.set_title("x2y Project X (Reverse)")
    ax5.set_xticks([])

    ax6 = plt.subplot(616)
    ax6.plot(x2x_dotplot_project_x_rev - x2y_dotplot_project_x_rev)
    ax6.set_title("Difference (x2x - x2y) (Reverse)")

    plt.tight_layout()
    plt.savefig(dotplot_output_prefix + ".projections.png")
    plt.close()
