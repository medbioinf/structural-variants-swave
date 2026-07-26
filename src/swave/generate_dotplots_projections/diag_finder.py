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

import numpy as np
from itertools import groupby

from .structures import Diag


def flip_matrix_lr(matrix):
    if sp.issparse(matrix):
        return matrix[:, ::-1]
    return np.fliplr(matrix)


def find_line_diag(matrix, min_line_len, return_format="dict", flip_lr=False, ignore_offset=None):
    rows, cols = matrix.shape

    if return_format == "dict":
        candidate_line_diags = {}
    else:
        candidate_line_diags = []

    for offset in range(-rows + 1, cols):

        if offset == ignore_offset:
            continue

        diag = matrix.diagonal(offset)

        current_pointer = 0

        for val, _ in groupby(diag):
            group_len = len(list(_))
            group_start = current_pointer
            group_end = current_pointer + group_len - 1

            if val == 1 and group_len >= min_line_len:
                if offset <= 0:
                    if flip_lr:
                        diag_obj = Diag(cols - (group_end + 1), cols - group_start, group_start + abs(offset), group_end + 1 + abs(offset), "reverse", offset)
                    else:
                        diag_obj = Diag(group_start, group_end + 1, group_start + abs(offset), group_end + 1 + abs(offset), "forward", offset)
                else:

                    if flip_lr:
                        diag_obj = Diag(cols - (group_end + 1 + offset), cols - (group_start + offset), group_start, group_end + 1, "reverse", offset)
                    else:
                        diag_obj = Diag(group_start + offset, group_end + 1 + offset, group_start, group_end + 1, "forward", offset)

                if return_format == "dict":
                    if offset not in candidate_line_diags:
                        candidate_line_diags[offset] = []

                    candidate_line_diags[offset].append(diag_obj)

                else:
                    candidate_line_diags.append(diag_obj)

            current_pointer = group_end + 1

    return candidate_line_diags


# def linear_diags_by_shift(major_offset, candidate_line_diags, neighbor_thresh, min_line_len,  n_rows, n_cols):
#     shift1_offset = major_offset + neighbor_thresh
#     shift2_offset = major_offset - neighbor_thresh

#     # do self linear
#     for offset in [major_offset, shift1_offset, shift2_offset]:
#         if offset not in candidate_line_diags:
#             continue

#         if offset != major_offset and offset == 0:
#             continue
#         offset_include_diags = candidate_line_diags[offset]

#         if len(offset_include_diags) > 1:
#             will_delete_diags = []
#             current_diag = offset_include_diags[0]
#             next_diag_pointer = 1

#             while next_diag_pointer < len(offset_include_diags):
#                 next_diag = offset_include_diags[next_diag_pointer]

#                 if next_diag.y_start - current_diag.y_end <= min_line_len:
#                     will_delete_diags.append(next_diag)

#                     current_diag.x_start = min(current_diag.x_start, next_diag.x_start)
#                     current_diag.y_start = min(current_diag.y_start, next_diag.y_start)

#                     current_diag.x_end = max(current_diag.x_end, next_diag.x_end)
#                     current_diag.y_end = max(current_diag.y_end, next_diag.y_end)
#                 else:
#                     current_diag = next_diag

#                 next_diag_pointer += 1

#             for diag in will_delete_diags:
#                 offset_include_diags.remove(diag)

#     # do shift linear
#     if major_offset in candidate_line_diags:
#         major_offset_include_diags = candidate_line_diags[major_offset]
        
#         for major_diag in major_offset_include_diags:
#             for shift_offset in [shift1_offset, shift2_offset]:
#                 if shift_offset not in candidate_line_diags:
#                     continue

#                 shift_offset_include_diags = candidate_line_diags[shift_offset]
#                 will_delete_diags = []

#                 for shift_diag in shift_offset_include_diags:
#                     # sort on reads (y)
#                     # target diag is the latter one
#                     if shift_diag.y_start >= major_diag.y_end:
#                         distance_on_y = shift_diag.y_start - major_diag.y_end
#                         distance_on_x = shift_diag.x_start - major_diag.x_end if major_diag.orient == "forward" else major_diag.x_start - shift_diag.x_end

#                     # base diag in the latter one
#                     elif major_diag.y_start >= shift_diag.y_end:
#                         distance_on_y = major_diag.y_start - shift_diag.y_end
#                         distance_on_x = major_diag.x_start - shift_diag.x_end if major_diag.orient == "forward" else shift_diag.x_start - major_diag.x_end

#                     # base and target are overlapped on y
#                     else:
#                         continue

#                     linear_thresh = min_line_len

#                     if distance_on_y <= linear_thresh and distance_on_x <= linear_thresh:
#                         major_diag.x_start = min([shift_diag.x_start, major_diag.x_start])
#                         major_diag.y_start = min([shift_diag.y_start, major_diag.y_start])
#                         major_diag.x_end = max([shift_diag.x_end, major_diag.x_end])
#                         major_diag.y_end = max([shift_diag.y_end, major_diag.y_end])

#                         major_diag_y_len = major_diag.y_end - major_diag.y_start
#                         major_diag_x_len = major_diag.x_end - major_diag.x_start

#                         if major_diag_x_len != major_diag_y_len:

#                             if major_diag_x_len > major_diag_y_len:

#                                 diff_length = major_diag_x_len - major_diag_y_len

#                                 if major_diag.y_end + diff_length >= n_rows:
#                                     major_diag.x_end -= diff_length
#                                 else:
#                                     major_diag.y_end += diff_length

#                             else:
#                                 diff_length = major_diag_y_len - major_diag_x_len

#                                 if major_diag.orient == "forward":
#                                     if major_diag.x_end + diff_length >= n_cols:
#                                         major_diag.y_end -= diff_length
#                                     else:
#                                         major_diag.x_end += diff_length
#                                 else:
#                                     if major_diag.x_start - diff_length < 0:
#                                         major_diag.y_start += diff_length
#                                     else:
#                                         major_diag.x_start -= diff_length

#                         will_delete_diags.append(shift_diag)

#                 for diag in will_delete_diags:
#                     shift_offset_include_diags.remove(diag)


def boost_diags(dotplot_type, unique_diags, x2y_matrix_len_x, x2y_matrix_len_y, raw_project_x, raw_project_x_rev, augment_coeff):
    unique_diags = sorted(unique_diags, key=lambda x: (x.orient, x.y_start))
    first_diag = sorted([diag for diag in unique_diags if diag.orient == "forward"], key=lambda x: (x.x_start + x.y_start))[0]
    last_diag = sorted([diag for diag in unique_diags if diag.orient == "forward"], key=lambda x: (x.x_end + x.y_end))[-1]
    other_diags = sorted([diag for diag in unique_diags if diag not in [first_diag, last_diag]], key=lambda x: (x.y_end - x.y_start), reverse=True)

    y_boost_flag = np.zeros(x2y_matrix_len_y)
    max_x_size = raw_project_x.shape[0]

    for diag in [first_diag, last_diag] + other_diags:  # for diag in unique_diags
        diag_x_positions = np.array(range(diag.x_start, diag.x_end))
        diag_x_positions_rev = np.flip(diag_x_positions)

        diag_y_positions = np.array(range(diag.y_start, diag.y_end))
        allowed_boost_y_positions = np.intersect1d(diag_y_positions,  np.where(y_boost_flag != 1))

        if diag.orient == "forward":
            idx = diag_x_positions[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))]
            valid_idx = idx[idx < max_x_size]
            raw_project_x[valid_idx] += augment_coeff
        else:
            idx_rev = diag_x_positions_rev[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))]
            valid_idx_rev = idx_rev[idx_rev < max_x_size]
            raw_project_x[valid_idx_rev] += augment_coeff
            
            # reverse diag (if it is also in the rev matrix, we consider it)
            if diag.true_reverse:
                raw_project_x_rev[diag_x_positions_rev[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))]] += augment_coeff

        y_boost_flag[allowed_boost_y_positions] = 1

    # last position is maintained with the previous one, since last position is commonly unmatched due to the potentially truncated sequence
    forward_diags = sorted([diag for diag in unique_diags if diag.orient == "forward"], key=lambda x: x.x_start)

    if not len(forward_diags) == 0:
        start_diag = forward_diags[0]
        start_trim_len = start_diag.x_start + 2

        raw_project_x[0: start_trim_len] = raw_project_x[start_trim_len]
        raw_project_x_rev[0: start_trim_len] = raw_project_x_rev[start_trim_len]

        forward_diags = sorted(forward_diags, key=lambda x: x.x_end)

        end_diag = forward_diags[-1]
        end_trim_len = x2y_matrix_len_x - end_diag.x_end + 2

        raw_project_x[-end_trim_len: -1] = raw_project_x[-end_trim_len]
        raw_project_x_rev[-end_trim_len: -1] = raw_project_x_rev[-end_trim_len]


def find_and_denoise_diags(x2x_matrix, x2y_matrix, x2y_matrix_rev, stride_size):
    x2y_matrix_len_x = x2y_matrix.shape[1]
    x2y_matrix_len_y = x2y_matrix.shape[0]

    max_shift_len = 50 / stride_size

    x2x_matrix_fliplr = flip_matrix_lr(x2x_matrix)
    x2y_matrix_fliplr = flip_matrix_lr(x2y_matrix)

    x2x_matrix_line_diags = find_line_diag(x2x_matrix, min_line_len=max(1, 30 / stride_size), return_format="dict", ignore_offset=0)
    x2x_matrix_line_diags_fliplr = find_line_diag(x2x_matrix_fliplr, min_line_len=max(1, 30 / stride_size), flip_lr=True, return_format="dict", ignore_offset=0)

    # allocate diags into bins to accelerate the compare process
    x2x_bin_size = 100
    x2x_matrix_line_diags_bins = {}

    for offset in x2x_matrix_line_diags:
        for diag in x2x_matrix_line_diags[offset]:
            diag_bin = int(diag.x_start / x2x_bin_size)

            if diag_bin not in x2x_matrix_line_diags_bins:
                x2x_matrix_line_diags_bins[diag_bin] = []
            x2x_matrix_line_diags_bins[diag_bin].append(diag)

    for offset in x2x_matrix_line_diags_fliplr:
        for diag in x2x_matrix_line_diags_fliplr[offset]:
            diag_bin = int(diag.x_start / x2x_bin_size)

            if diag_bin not in x2x_matrix_line_diags_bins:
                x2x_matrix_line_diags_bins[diag_bin] = []
            x2x_matrix_line_diags_bins[diag_bin].append(diag)

    unique_diags = []

    for mode in ["raw", "fliplr"]:
        if mode == "raw":
            x2y_matrix_line_diags = find_line_diag(x2y_matrix, min_line_len=max(5, 50 / stride_size), return_format="dict")
        else:
            x2y_matrix_line_diags = find_line_diag(x2y_matrix_fliplr, min_line_len=max(5, 50 / stride_size), flip_lr=True, return_format="dict")

        for offset in x2y_matrix_line_diags:
            for diag in x2y_matrix_line_diags[offset]:
                # keep first and last diags, as they are the main bone of the dotplot
                if ((diag.x_start in [0] and diag.y_start in [0])
                        or (diag.x_end in [x2y_matrix_len_x - 1, x2y_matrix_len_x] and diag.y_end in [x2y_matrix_len_y - 1, x2y_matrix_len_y])):
                    unique_diags.append(diag)
                    continue

                diag_start, diag_end = diag.x_start, diag.x_end
                match_flag = False

                diag_bin = int(diag.x_start / x2x_bin_size)
                if diag_bin in x2x_matrix_line_diags_bins:
                    for target_diag in x2x_matrix_line_diags_bins[diag_bin]:
                        if target_diag.x_start - max_shift_len <= diag_start <= diag_end <= target_diag.x_end + max_shift_len:
                            match_flag = True
                            break

                if match_flag is False:
                    unique_diags.append(diag)

    # sort
    unique_diags = sorted(unique_diags, key=lambda x: x.y_start)

    # for the reversed diags in x2y_matrix: check if they are really in the x2y_matrix_rev
    x2y_matrix_rev_line_diags = find_line_diag(flip_matrix_lr(x2y_matrix_rev), min_line_len=max(5, 50 / stride_size), flip_lr=True, return_format="dict")

    x2y_matrix_rev_line_diags_bins = {}
    for offset in x2y_matrix_rev_line_diags:
        for diag in x2y_matrix_rev_line_diags[offset]:
            diag_bin = int(diag.x_start / x2x_bin_size)
            
            if diag_bin not in x2y_matrix_rev_line_diags_bins:
                x2y_matrix_rev_line_diags_bins[diag_bin] = []
            x2y_matrix_rev_line_diags_bins[diag_bin].append(diag)

    for diag in unique_diags:
        if diag.orient == "forward":
            continue

        diag_bin = int(diag.x_start / x2x_bin_size)
        if diag_bin in x2y_matrix_rev_line_diags_bins:
            for target_diag_rev in x2y_matrix_rev_line_diags_bins[diag_bin]:
                if diag == target_diag_rev:
                    diag.true_reverse = True

    if len(unique_diags) == 0:
        return unique_diags

    # remove diags that are overlapped with the left and right-most diags (anchor diags)
    first_diag = sorted([diag for diag in unique_diags if diag.orient == "forward"], key=lambda x: (x.x_start + x.y_start))[0]
    last_diag = sorted([diag for diag in unique_diags if diag.orient == "forward"], key=lambda x: (x.x_end + x.y_end))[-1]

    overlapped_diags = []
    for diag in unique_diags:
        if diag == first_diag or diag == last_diag:
            continue
        if diag.orient == "forward" and (diag.y_start < first_diag.y_end or diag.y_end > last_diag.y_start):
            overlapped_diags.append(diag)

    for diag in overlapped_diags:
        unique_diags.remove(diag)

    return unique_diags
