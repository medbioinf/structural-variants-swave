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

import math
import numpy as np

from .structures import Diag
from src.swave.utils import reverse_complement_seq, is_kmer_similar, calculate_seq_similarity_larger_than


def base_level_dotplot_extension(unique_diags, x_seq, y_seq, x2y_matrix, raw_project_x, raw_project_x_rev, raw_kmer_size, stride_size):
    x2y_matrix_shape = x2y_matrix.shape
    x2y_matrix_len_x = x2y_matrix_shape[1]
    x2y_matrix_len_y = x2y_matrix_shape[0]

    re_check_kmer_size = 10
    re_check_length_start = 1 * raw_kmer_size
    re_check_length_end = 1 * raw_kmer_size

    # perform base level boost
    for i in range(len(unique_diags)):
        diag = unique_diags[i]
        try:
            if i == 0:
                max_augment_start = np.inf
                max_augment_end = unique_diags[i + 1].x_start - unique_diags[i].x_end

            elif i == len(unique_diags) - 1:
                max_augment_start = unique_diags[i].x_start - unique_diags[i - 1].x_end
                max_augment_end = np.inf
            else:
                max_augment_start = unique_diags[i].x_start - unique_diags[i - 1].x_end
                max_augment_end = unique_diags[i + 1].x_start - unique_diags[i].x_end

            # < 0 means there is a candidate overlap between current diag and its previous or later one
            if max_augment_start < 0:
                max_augment_start = np.inf
            if max_augment_end < 0:
                max_augment_end = np.inf

        except:
            max_augment_start = np.inf
            max_augment_end = np.inf

        # adjust the coordinate by the orientataion of the diag
        if diag.orient == "forward":
            diag_x_start = diag.x_start
            diag_x_end = diag.x_end
            diag_y_start = diag.y_start
            diag_y_end = diag.y_end
        else:
            diag_x_start = diag.x_start
            diag_x_end = diag.x_end
            diag_y_start = x2y_matrix_len_y - diag.y_end - 1
            diag_y_end = x2y_matrix_len_y - diag.y_start - 1

            y_seq = reverse_complement_seq(y_seq)

        # deal the start seq of the diag line
        if not (diag_x_start - 1 < 0 and diag_y_start - 1 < 0):

            diag_start_seq_x = x_seq[diag_x_start * stride_size - re_check_length_start: diag_x_start * stride_size]
            diag_start_seq_y = y_seq[diag_y_start * stride_size - re_check_length_start: diag_y_start * stride_size]

            re_check_length_start = min(re_check_length_start, len(diag_start_seq_x), len(diag_start_seq_y))
            diag_start_match_flag = False

            # traverse on seq_x
            max_match_length, index_on_x, index_on_y = 0, 0, 0
            for index_on_x in range(re_check_length_start, re_check_kmer_size, -1):
                # traverse on seq_y
                for index_on_y in range(re_check_length_start, re_check_kmer_size, -1):
                    if not (index_on_x == re_check_length_start or index_on_y == re_check_length_start):
                        continue

                    diag_start_seq_x_sub_seq = diag_start_seq_x[index_on_x - re_check_kmer_size: index_on_x]
                    diag_start_seq_y_sub_seq = diag_start_seq_y[index_on_y - re_check_kmer_size: index_on_y]

                    # match the sub seq
                    if is_kmer_similar(diag_start_seq_x_sub_seq, diag_start_seq_y_sub_seq, mismatch_thresh=1):
                        diag_start_match_flag = True

                        # find the max match length
                        max_match_length += re_check_kmer_size
                        max_match_flag = "meet_sub_seq_end"  # choose from meet_sub_seq_end or meet_mis_match
                        for i in range(min(index_on_x, index_on_y)):
                            if diag_start_seq_x[index_on_x - re_check_kmer_size - i] == diag_start_seq_y[index_on_y - re_check_kmer_size - i]:
                                max_match_length += 1
                            else:
                                max_match_flag = "meet_mis_match"
                                break

                        # update the diag
                        if max_match_flag == "meet_mis_match":
                            augment_length = min(math.ceil(max_match_length / stride_size), max_augment_start)
                        else:
                            # # if meet_sub_seq_end, then the whole sub seq is matched although it is not long enough, so we consider this as the whole match of the sub seq
                            augment_length = min(math.ceil(re_check_length_start / stride_size), max_augment_start)

                        for i in range(augment_length):
                            if diag.x_start - 1 >= 0:

                                # the start is an open-interval, therefore, we first -= 1 then update the project
                                diag.x_start -= 1
                                diag.y_start -= 1

                                raw_project_x[diag.x_start] += 1
                                if diag.orient == "reverse":
                                    raw_project_x_rev[diag.x_start] += 1

                        break

                if diag_start_match_flag is True:
                    break

        if not (diag_x_end + 1 >= x2y_matrix_len_x and diag_y_end + 1 >= x2y_matrix_len_y):

            diag_end_seq_x = x_seq[diag_x_end * stride_size: diag_x_end * stride_size + re_check_length_end]
            diag_end_seq_y = y_seq[diag_y_end * stride_size: diag_y_end * stride_size + re_check_length_end]

            re_check_length_end = min(re_check_length_end, len(diag_end_seq_x), len(diag_end_seq_y))
            diag_end_match_flag = False

            max_match_length, index_on_x, index_on_y = 0, 0, 0
            for index_on_x in range(re_check_length_end - re_check_kmer_size):
                # traverse on seq_y
                for index_on_y in range(re_check_length_end - re_check_kmer_size):
                    if not (index_on_x == 0 or index_on_y == 0):
                        continue

                    diag_end_seq_x_sub_seq = diag_end_seq_x[index_on_x: index_on_x + re_check_kmer_size]
                    diag_end_seq_y_sub_seq = diag_end_seq_y[index_on_y: index_on_y + re_check_kmer_size]

                    # match the sub seq
                    if is_kmer_similar(diag_end_seq_x_sub_seq, diag_end_seq_y_sub_seq, mismatch_thresh=1):
                        diag_end_match_flag = True

                        # find the max match length
                        max_match_length += re_check_kmer_size
                        max_match_flag = "meet_sub_seq_end"  # choose from meet_sub_seq_end or meet_mis_match

                        for i in range(re_check_length_end - max(index_on_x, index_on_y) - re_check_kmer_size):
                            if diag_end_seq_x[index_on_x + re_check_kmer_size + i] == diag_end_seq_y[index_on_y + re_check_kmer_size + i]:
                                max_match_length += 1
                            else:
                                max_match_flag = "meet_mis_match"
                                break

                        # update the diag
                        if max_match_flag == "meet_mis_match":
                            augment_length = min(math.ceil(max_match_length / stride_size), max_augment_end)
                        else:
                            augment_length = min(math.ceil(re_check_length_end / stride_size), max_augment_end)
                            
                        # update the diag and boost the raw project
                        for i in range(augment_length):
                            if diag.x_end < x2y_matrix_len_x:

                                # the end is a close-interval, therefore, we first update the project then += 1
                                raw_project_x[diag.x_end] += 1
                                if diag.orient == "reverse":
                                    raw_project_x_rev[diag.x_end] += 1

                                diag.x_end += 1
                                diag.y_end += 1

                        break

                if diag_end_match_flag is True:
                    break

    return unique_diags


def linear_level_dotplot_extension(unique_diags, x2y_matrix):
    x2y_matrix_shape = x2y_matrix.shape
    x2y_matrix_len_x = x2y_matrix_shape[1]
    x2y_matrix_len_y = x2y_matrix_shape[0]

    # sort
    unique_diags = sorted(unique_diags, key=lambda x: x.y_start)

    # check linear diags (why: SNPs would case break of diags, which would be predicted as del + inv)
    linear_diags = []

    for i in range(len(unique_diags) - 1, 0, -1):
        current_diag = unique_diags[i]
        current_diag_len = current_diag.y_end - current_diag.y_start

        previous_diag = unique_diags[i - 1]
        previous_diag_len = previous_diag.y_end - previous_diag.y_start

        if current_diag.orient != previous_diag.orient:
            continue

        if (current_diag.x_start == 0 and current_diag.y_start == 0) or (previous_diag.x_start == 0 and previous_diag.y_start == 0):
            continue

        linear_thresh = max([current_diag_len, previous_diag_len])

        if current_diag.orient == "reverse" and ((current_diag.y_start - previous_diag.y_end) == (previous_diag.x_start - current_diag.x_end)) and (current_diag.y_start - previous_diag.y_end) <= linear_thresh:
            previous_diag.x_start = min([previous_diag.x_start, current_diag.x_start])
            previous_diag.y_start = min([previous_diag.y_start, current_diag.y_start])
            previous_diag.x_end = max([previous_diag.x_end, current_diag.x_end])
            previous_diag.y_end = max([previous_diag.y_end, current_diag.y_end])

            linear_diags.append(current_diag)

        if current_diag.orient == "forward" and ((current_diag.y_start - previous_diag.y_end) == (current_diag.x_start - previous_diag.x_end)) and (current_diag.y_start - previous_diag.y_end) <= linear_thresh:
            previous_diag.x_start = min([previous_diag.x_start, current_diag.x_start])
            previous_diag.y_start = min([previous_diag.y_start, current_diag.y_start])
            previous_diag.x_end = max([previous_diag.x_end, current_diag.x_end])
            previous_diag.y_end = max([previous_diag.y_end, current_diag.y_end])

            linear_diags.append(current_diag)

    for diag in linear_diags:
        unique_diags.remove(diag)

    # remove diags that are fully covered by others
    unique_diags = sorted(unique_diags, key=lambda x: (x.y_end - x.y_start), reverse=True)

    full_covered_diags = []
    for i in range(len(unique_diags) - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            base_diag = unique_diags[i]

            if ((base_diag.x_start in [0] and base_diag.y_start in [0])
                    or (base_diag.x_end in [x2y_matrix_len_x - 1, x2y_matrix_len_x] and base_diag.y_end in [x2y_matrix_len_y - 1, x2y_matrix_len_y])):
                continue

            target_diag = unique_diags[j]

            if target_diag.y_start - 1 <= base_diag.y_start <= base_diag.y_end <= target_diag.y_end + 1:
                full_covered_diags.append(base_diag)
                break

    for diag in full_covered_diags:
        unique_diags.remove(diag)

    return unique_diags


def kmer10_level_dotplot_extension(unique_diags, x_seq, y_seq, x2y_matrix, stride_size, thresh=0.5):
    """
    Extends diags in repetitive regions by checking for similar matches using a detailed dotplot with small kmers and strides.
    """
    unique_diags = sorted(unique_diags, key=lambda x: x.y_start)

    x2y_matrix_shape = np.shape(x2y_matrix)
    x2y_matrix_shape_x = x2y_matrix_shape[1]
    x2y_matrix_shape_y = x2y_matrix_shape[0]

    unique_diags.insert(0, Diag(0, 0, 0, 0, "forward", 0))
    unique_diags.append(Diag(x2y_matrix_shape_x, x2y_matrix_shape_x, x2y_matrix_shape_y, x2y_matrix_shape_y, "forward", 0))

    for diag_index in range(1, len(unique_diags) - 1):
        current_diag = unique_diags[diag_index]

        if current_diag.x_start == 0 and current_diag.y_start == 0:
            continue

        # deal with the previous gap
        previous_diag = unique_diags[diag_index - 1]
        gap_length = current_diag.y_start - previous_diag.y_end

        if gap_length > 0:
            gap_start_on_y = previous_diag.y_end
            gap_end_on_y = current_diag.y_start

            if current_diag.orient == "forward":
                gap_start_on_x = max(current_diag.x_start - gap_length, 0)
                gap_end_on_x = current_diag.x_start
            else:
                gap_start_on_x = current_diag.x_end
                gap_end_on_x = min(current_diag.x_end + gap_length, x2y_matrix_shape_x)

            gap_seq_on_y = y_seq[gap_start_on_y * stride_size: gap_end_on_y * stride_size]
            gap_seq_on_x = x_seq[gap_start_on_x * stride_size: gap_end_on_x * stride_size]

            if gap_start_on_y == gap_end_on_y or gap_start_on_x == gap_end_on_x:
                continue

            if current_diag.orient == "forward":
                similar_len = calculate_seq_similarity_larger_than(thresh, current_diag.orient, gap_seq_on_x, gap_seq_on_y, "left")
            else:
                similar_len = calculate_seq_similarity_larger_than(thresh, current_diag.orient, gap_seq_on_x, gap_seq_on_y, "right")

            if current_diag.orient == "forward":
                current_diag.x_start -= int(similar_len / stride_size)
                current_diag.y_start -= int(similar_len / stride_size)
            else:
                current_diag.x_end += int(similar_len / stride_size)
                current_diag.y_start -= int(similar_len / stride_size)

        # deal with the latter gap
        latter_diag = unique_diags[diag_index + 1]
        gap_length = latter_diag.y_start - current_diag.y_end

        if gap_length > 0:
            gap_start_on_y = current_diag.y_end
            gap_end_on_y = latter_diag.y_start

            if current_diag.orient == "forward":
                gap_start_on_x = current_diag.x_end
                gap_end_on_x = min(current_diag.x_end + gap_length, x2y_matrix_shape_x)
            else:
                gap_start_on_x = max(current_diag.x_start - gap_length, 0)
                gap_end_on_x = current_diag.x_start

            gap_seq_on_y = y_seq[gap_start_on_y * stride_size: gap_end_on_y * stride_size]
            gap_seq_on_x = x_seq[gap_start_on_x * stride_size: gap_end_on_x * stride_size]

            if gap_start_on_y == gap_end_on_y or gap_start_on_x == gap_end_on_x:
                continue

            if current_diag.orient == "forward":
                similar_len = calculate_seq_similarity_larger_than(thresh, current_diag.orient, gap_seq_on_x, gap_seq_on_y, "right")
            else:
                similar_len = calculate_seq_similarity_larger_than(thresh, current_diag.orient, gap_seq_on_x, gap_seq_on_y, "left")

            if current_diag.orient == "forward":
                current_diag.x_end += int(similar_len / stride_size)
                current_diag.y_end += int(similar_len / stride_size)
            else:
                current_diag.x_start -= int(similar_len / stride_size)
                current_diag.y_end += int(similar_len / stride_size)

    return unique_diags[1: -1]


def rotate_diags_to_alt2ref(unique_diags):
    rotated_diags = []

    for diag in unique_diags:
        new_diag = Diag(diag.y_start, diag.y_end, diag.x_start, diag.x_end, diag.orient, diag.offset)
        new_diag.true_reverse = diag.true_reverse

        rotated_diags.append(new_diag)

    return rotated_diags
