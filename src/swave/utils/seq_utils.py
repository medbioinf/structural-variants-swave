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

import itertools
import numpy as np


def reverse_complement_seq(seq):
    """
    Returns the reverse complement of a DNA sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def calculate_stride_size(seq_x, seq_y):
    """
    Calculates the stride size for dotplot alignment to keep the final matrix dimension manageable.
    """
    max_seq_len = max([len(seq_x), len(seq_y)])

    if max_seq_len <= 1500:
        stride_size = 1
    elif max_seq_len <= 15000:
        stride_size = 10
    elif max_seq_len <= 75000:
        stride_size = 30
    elif max_seq_len <= 100000:
        stride_size = 50
    else:
        stride_size = 100

    return stride_size


def is_kmer_similar(kmer1, kmer2, mismatch_thresh=0):
    """
    Compares two k-mers and returns True if they are similar within a specified mismatch threshold.
    """
    if len(kmer1) != len(kmer2):
        return False
    mismatch_cnt = sum(b1 != b2 for b1, b2 in zip(kmer1, kmer2))
    return mismatch_cnt <= mismatch_thresh


def find_continuous_val(val_list):
    """
    Find continuous values in a list or string.

    return format: [ [val, [indexes in raw list]], [], ..... ]

    input example:  AAABBBCCAAA
    output example : [['A', [0, 1, 2]], ['B', [3, 4, 5]], ['C', [6, 7]], ['A', [8, 9, 10]]]
    """
    continuous_res = []
    current_pointer = 0

    for val, group in itertools.groupby(val_list):
        group_len = len(list(group))
        group_start = current_pointer               # open interval
        group_end = current_pointer + group_len     # closed interval
        
        continuous_res.append([val, [i for i in range(group_start, group_end)]])

        current_pointer = group_end

    return continuous_res


def calculate_seq_similarity(seq1, seq2):
    if len(seq1) == 0 or len(seq2) == 0:
        return 0

    # use dotplot projection for sequence similarity
    from src.generate_dotplots_projections_mod.structures import Dotplot

    # use the shorter seq as x_seq, and projection to the shorter seq
    if len(seq1) < len(seq2):
        dotplot = Dotplot(seq1, seq2, 10, "tmp", stride_size=1, given_x_kmer_index=None)
    else:
        dotplot = Dotplot(seq2, seq1, 10, "tmp", stride_size=1, given_x_kmer_index=None)

    projection_x, _ = dotplot.get_project_x()
    projection_x_rev = dotplot.get_project_x_rev()

    return np.count_nonzero(projection_x + projection_x_rev) / len(projection_x)


def calculate_seq_similarity_larger_than(thresh, seq_orient, seq1, seq2, extend_orient):
    from src.generate_dotplots_projections_mod.structures import Dotplot

    if seq_orient == "forward":
        dotplot = Dotplot(seq1, seq2, 10, str(len(seq1)) + str(len(seq2)) + extend_orient, stride_size=None, given_x_kmer_index=None, skip_reverse=True)
    else:
        dotplot = Dotplot(seq1, seq2, 10, str(len(seq1)) + str(len(seq2)) + extend_orient, stride_size=None, given_x_kmer_index=None, skip_forward=True)

    detail_stride_size = dotplot.stride_size

    projection_x, _ = dotplot.get_project_x()
    projection_x_len = len(projection_x)
    projection_y, _ = dotplot.get_project_y()
    projection_x_summary = find_continuous_val(projection_x)
    
    if extend_orient == "left":
        for value, include_indexes in projection_x_summary:
            if value == 0:
                continue

            index_start = include_indexes[0]
            remaining_len = projection_x_len - index_start

            if remaining_len == 0:
                return 0

            if np.count_nonzero(projection_x[index_start: ]) / remaining_len >= thresh:
                # check the projection in y
                if seq_orient == "forward" and np.count_nonzero(projection_y[-remaining_len:]) / remaining_len >= thresh:
                    return remaining_len * detail_stride_size
                if seq_orient == "reverse" and np.count_nonzero(projection_y[: remaining_len]) / remaining_len >= thresh:
                    return remaining_len * detail_stride_size
    
    if extend_orient == "right":
        projection_x_summary.reverse()
        for value, include_indexes in projection_x_summary:
            if value == 0:
                continue

            index_end = include_indexes[-1]
            remaining_len = index_end

            if remaining_len == 0:
                return 0

            if np.count_nonzero(projection_x[: index_end + 1]) / remaining_len >= thresh:
                # check the projection in y
                if seq_orient == "forward" and np.count_nonzero(projection_y[: remaining_len]) / remaining_len >= thresh:
                    return remaining_len * detail_stride_size
                if seq_orient == "reverse" and np.count_nonzero(projection_y[-remaining_len:]) / remaining_len >= thresh:
                    return remaining_len * detail_stride_size

    return 0


def calculate_seq_repeat_ratio(seq1):
    from src.generate_dotplots_projections_mod.structures import Dotplot

    if len(seq1) == 0:
        return np.inf

    dotplot = Dotplot(seq1, seq1, 10, "tmp", stride_size=None, given_x_kmer_index=None)
    project_y, _ = dotplot.get_project_y()
    
    return np.average(project_y)
