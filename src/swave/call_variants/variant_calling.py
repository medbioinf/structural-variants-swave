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

from .structures import SSV, CSV
from src.swave.utils import calculate_seq_similarity, calculate_seq_repeat_ratio, find_continuous_val

label2index = {"None": 0, "REF": 1, "DEL": 2, "INV": 3, "DUP": 4, "invDUP": 5}
index2label = {0: "None", 1: "REF", 2: "DEL", 3: "INV", 4: "DUP", 5: "invDUP"}


def fetch_ssv_comp_seq(ssv_comp, dotplot_ref_start, dotplot_bundle, max_seq_len=None): 
    x2y_ref2alt = dotplot_bundle["x2y_ref2alt"]

    if ssv_comp.type == "INS":
        seq_start = ssv_comp.raw_start - dotplot_ref_start
        seq_end = ssv_comp.raw_start + ssv_comp.length - dotplot_ref_start
        
        if max_seq_len is not None:
            seq_end = min(seq_end, seq_start + max_seq_len)
            
        return x2y_ref2alt.seq_y[seq_start:seq_end]
    else:
        seq_start = ssv_comp.start - dotplot_ref_start
        seq_end = ssv_comp.end - dotplot_ref_start + 1

        if max_seq_len is not None:
            seq_end = min(seq_end, seq_start + max_seq_len)
            
        return x2y_ref2alt.seq_x[seq_start:seq_end]


def prediction_to_ssv_comps(matrix, prediction, dotplot_ref_chr, dotplot_ref_start, dotplot_stride_size, options):
    ssv_components = []
    for pred_label_index, include_indexes in find_continuous_val(prediction):
        # 0: None, 1: REF, they are not SVs, so we skip
        if pred_label_index == 0 or pred_label_index == 1:
            continue

        ssv_type = index2label[pred_label_index]
        ssv_start = int(matrix[include_indexes[0]][0])  # first seg's start pos
        ssv_end = int(matrix[include_indexes[-1]][1]) - 1  # last seg's end pos. why -1: the matrix is 'Open left and close right'

        ssv_obj = SSV(ssv_type, dotplot_ref_chr, dotplot_ref_start + (ssv_start * dotplot_stride_size), dotplot_ref_start + (ssv_end * dotplot_stride_size), )

        if options.max_sv_size >= ssv_obj.length >= max(options.min_sv_size, 2 * dotplot_stride_size):
            ssv_obj.set_raw_projection([matrix[i] for i in include_indexes])
            ssv_components.append(ssv_obj)

    return ssv_components


def eval_csv_quality(csv, dotplot_bundle):    
    e = 0.0001
    x2y_ref2alt = dotplot_bundle["x2y_ref2alt"]
    
    ref_seq = x2y_ref2alt.seq_x
    alt_seq = x2y_ref2alt.seq_y

    seq_effective_length = len(alt_seq) - len(ref_seq)

    csv_effective_length = 0
    for ssv in csv.ssv_components:
        if ssv.type in ["INS", "DUP", "invDUP"]:
            csv_effective_length += ssv.length
        elif ssv.type in ["DEL"]:
            csv_effective_length -= ssv.length

    if -50 <= abs(csv_effective_length) - abs(seq_effective_length) <= 50:
        return 1
    else:
        return round((csv_effective_length + e) / (seq_effective_length + e), 3)


def handle_continuous_ssv_types(ssv_comps):
    ssv_comps = sorted(ssv_comps, key=lambda x: x.start)
    will_removed_ssv_comp = set()
    
    for ssv_comp_index in range(len(ssv_comps) - 1, 0, -1):
        ssv_comp = ssv_comps[ssv_comp_index]
        previous_ssv_comp = ssv_comps[ssv_comp_index - 1]
        
        if ssv_comp.type != previous_ssv_comp.type:
            continue

        is_dup_or_inv = ssv_comp.type in ["INV", "DUP", "invDUP"] and abs(ssv_comp.start - previous_ssv_comp.end) < 1000
        
        is_del = ssv_comp.type in ["DEL"] and abs(ssv_comp.start - previous_ssv_comp.end) < 100
        
        if is_dup_or_inv or is_del:
            previous_ssv_comp.end = ssv_comp.end
            previous_ssv_comp.length = previous_ssv_comp.end - previous_ssv_comp.start + 1
            will_removed_ssv_comp.add(ssv_comp)

    return [comp for comp in ssv_comps if comp not in will_removed_ssv_comp]


def merge_ref2alt_alt2ref_csv(ref2alt_comps, alt2ref_comps, dotplot_id, dotplot_bundle, options):
    id_parts = dotplot_id.split("|")
    dotplot_ref_start = int(id_parts[4])
    dotplot_ref_end = int(id_parts[5])
    
    dotplot_stride_size = dotplot_bundle["stride_size"]

    ref2alt_comps = [comp for comp in handle_continuous_ssv_types(ref2alt_comps) if comp.length >= options.min_sv_size]

    ref2alt_comps = deal_with_dup(ref2alt_comps, dotplot_ref_start, dotplot_bundle, options)

    ref2alt_comps = deal_with_ins_from_alt2ref(ref2alt_comps, alt2ref_comps, dotplot_ref_start, dotplot_ref_end, dotplot_bundle, options)

    ref2alt_comps = find_dup_insert_pos(ref2alt_comps, dotplot_stride_size, options)

    ref2alt_comps = deal_with_multi_ins_del(ref2alt_comps, dotplot_ref_start, dotplot_bundle, options)

    csv = CSV(handle_continuous_ssv_types(ref2alt_comps))
    
    csv.set_score(eval_csv_quality(csv, dotplot_bundle))

    if len(csv.ssv_components) > options.max_sv_comps and is_scarred_inv(csv.type) is False:
        ssv_type_list = [ssv.type for ssv in csv.ssv_components]

        if "INV" in ssv_type_list:
            csv.type = "hyperCPX_INV"
        elif "invDUP" in ssv_type_list:
            csv.type = "hyperCPX_invDUP"
        elif "DUP" in ssv_type_list:
            csv.type = "hyperCPX_DUP"
        elif "DEL" in ssv_type_list and "INS" in ssv_type_list:
            csv.type = "hyperCPX_INS+DEL"
        elif "DEL" in ssv_type_list:
            csv.type = "hyperCPX_DEL"
        elif "INS" in ssv_type_list:
            csv.type = "hyperCPX_INS"
        else:
            csv.type = "hyperCPX"

        csv.ssv_components = [SSV(csv.type, csv.chr, csv.start, csv.end)]

    return csv


def is_scarred_inv(sv_type):
    sv_type_split = sv_type.split("+")
    if not len(sv_type_split) >= 3:
        return False

    if sv_type_split[0] == "INV":
        for i in range(len(sv_type_split)):
            if i % 2 == 0 and sv_type_split[i] != "INV":
                return False
            if i % 2 != 0 and sv_type_split[i] not in ["INS", "DEL"]:
                return False
    else:
        for i in range(len(sv_type_split)):
            if i % 2 != 0 and sv_type_split[i] != "INV":
                return False
            if i % 2 == 0 and sv_type_split[i] not in ["INS", "DEL"]:
                return False

    return True


def deal_with_ins_from_alt2ref(ref2alt_comps, alt2ref_comps, dotplot_ref_start, dotplot_ref_end, dotplot_bundle, options):
    dup_exist = False
    for ref2alt_ssv in ref2alt_comps:
        if ref2alt_ssv.type in ["DUP"]:
            dup_exist = True
            break

    ins_within_inv = []

    for i in range(len(alt2ref_comps)):
        # ins is del in alt2ref
        if alt2ref_comps[i].type == "DEL":
            if i == 0 or i == len(alt2ref_comps) - 1:
                continue
            try:
                if alt2ref_comps[i - 1].type in ["INV", "invDUP"] and alt2ref_comps[i + 1].type in ["INV", "invDUP"]:
                    ins_within_inv.append(alt2ref_comps[i])
            except IndexError:
                pass

    alt2ref_comps = [comp for comp in handle_continuous_ssv_types(alt2ref_comps) if comp.type == "DEL" and comp.length >= options.min_sv_size]
    
    # traverse remained ssv in alt2ref, for capturing the INSs that are not reflected in ref2alt
    for alt2ref_ssv in alt2ref_comps:
        new_start = alt2ref_ssv.start
        raw_start = alt2ref_ssv.start
        
        for tmp_ssv in alt2ref_comps:
            if tmp_ssv.type == "DEL" and raw_start >= tmp_ssv.end:
                new_start -= tmp_ssv.length

        raw_start = new_start
        for tmp_ssv in ref2alt_comps:
            if tmp_ssv.type == "DEL" and raw_start >= tmp_ssv.end:
                new_start += tmp_ssv.length

        if options.min_sv_size <= alt2ref_ssv.length <= options.max_sv_size:
            new_ssv = SSV("INS", alt2ref_ssv.chr, new_start, new_start, alt_seq="N" * alt2ref_ssv.length)

            # check if this insertion is covered by a inversion
            if alt2ref_ssv in ins_within_inv:
                new_ssv = SSV("INS", alt2ref_ssv.chr, dotplot_ref_start + (dotplot_ref_end - new_start), dotplot_ref_start + (dotplot_ref_end - new_start), alt_seq="N" * alt2ref_ssv.length)

                # split the whole inv that is cutted by the ins
                for target_ssv in ref2alt_comps:
                    if target_ssv.type == "INV" and target_ssv.start + 1 <= new_ssv.start <= new_ssv.end <= target_ssv.end - 1:
                        split_inv_left = SSV("INV", target_ssv.chr, target_ssv.start, new_ssv.start - 1)
                        split_inv_right = SSV("INV", target_ssv.chr, new_ssv.end + 1, target_ssv.end)
                        
                        if options.min_sv_size <= split_inv_left.length <= options.max_sv_size:
                            ref2alt_comps.append(split_inv_left)
                        if options.min_sv_size <= split_inv_right.length <= options.max_sv_size:
                            ref2alt_comps.append(split_inv_right)
                        ref2alt_comps.remove(target_ssv)
                        break

                new_ssv.set_raw_cords(alt2ref_ssv.chr, alt2ref_ssv.start, alt2ref_ssv.end)
                new_ssv.within_inv = True

                if dup_exist:
                    new_ssv_seq = fetch_ssv_comp_seq(new_ssv, dotplot_ref_start, dotplot_bundle, max_seq_len=5000)
                    repeat_ratio = calculate_seq_repeat_ratio(new_ssv_seq)
                    if repeat_ratio > 1.2:
                        new_ssv.within_repeat = True
                    # print(repeat_ratio, len(new_ssv_seq), new_ssv.to_string())

                ref2alt_comps.append(new_ssv)

            else:
                covered_flag = False
                for target_ssv in ref2alt_comps:
                    if target_ssv.type == "INV" and target_ssv.start + 1 <= new_ssv.start <= new_ssv.end <= target_ssv.end - 1:
                        covered_flag = True
                        break

                if not covered_flag:
                    new_ssv.set_raw_cords(alt2ref_ssv.chr, alt2ref_ssv.start, alt2ref_ssv.end)

                    if dup_exist:
                        new_ssv_seq = fetch_ssv_comp_seq(new_ssv, dotplot_ref_start, dotplot_bundle, max_seq_len=5000)
                        repeat_ratio = calculate_seq_repeat_ratio(new_ssv_seq)
                        if repeat_ratio > 1.2:
                            new_ssv.within_repeat = True

                    ref2alt_comps.append(new_ssv)

    return ref2alt_comps


def find_dup_insert_pos(ref2alt_comps, dotplot_stride_size, options):
    dup_comp_no_insert = []
    consumed_ins = []
    
    for dup_comp in ref2alt_comps:
        if "DUP" not in dup_comp.type:
            continue

        # determine the inserted pos of this dup by finding the most likely insertion in alt2ref
        best_match_ins_index = -1
        best_match_ins_score = -1

        for index in range(len(ref2alt_comps)):
            if ref2alt_comps[index].type != "INS":
                continue

            if ref2alt_comps[index].length < max(options.min_sv_size, 2 * dotplot_stride_size):
                continue

            if ref2alt_comps[index].within_repeat is True:
                continue

            ins_comp = ref2alt_comps[index]
            length_thresh = max(options.min_sv_size, 2 * dotplot_stride_size)
            
            if dup_comp.type == "DUP":
                # exact length similarity
                if abs(dup_comp.length - ins_comp.length) < length_thresh:
                    cur_match_score = 1 - abs(1 - dup_comp.length / ins_comp.length)
                    if cur_match_score > best_match_ins_score:
                        best_match_ins_index = index
                        best_match_ins_score = cur_match_score
                else:
                    # is a tandem dup
                    if abs(dup_comp.start - ins_comp.start) < length_thresh or abs(dup_comp.end - ins_comp.start) < length_thresh:
                        cur_match_score = 1 - abs(1 - dup_comp.length / ins_comp.length)
                        if cur_match_score > best_match_ins_score:
                            best_match_ins_index = index
                            best_match_ins_score = cur_match_score
            else:
                cur_match_score = 1 - abs(1 - dup_comp.length / ins_comp.length)
                if cur_match_score > best_match_ins_score:
                    best_match_ins_index = index
                    best_match_ins_score = cur_match_score

        if best_match_ins_index != -1:
            dup_comp.refine_dup_events(ref2alt_comps[best_match_ins_index].chr, ref2alt_comps[best_match_ins_index].start)
            
            ref2alt_comps[best_match_ins_index].length -= dup_comp.length

            if ref2alt_comps[best_match_ins_index].length < options.min_sv_size:
                consumed_ins.append(ref2alt_comps[best_match_ins_index])
        else:
            dup_comp_no_insert.append(dup_comp)

    for dup_comp in dup_comp_no_insert:
        ref2alt_comps.remove(dup_comp)

    for ins_comp in consumed_ins:
        ref2alt_comps.remove(ins_comp)

    return ref2alt_comps


def deal_with_dup(ref2alt_comps, dotplot_ref_start, dotplot_bundle, options):
    # check invdup between inversions
    for comp_index in range(len(ref2alt_comps)):
        if ref2alt_comps[comp_index].type == "invDUP":
            if comp_index == 0 or comp_index == len(ref2alt_comps) - 1:
                continue
            try:
                if ((ref2alt_comps[comp_index - 1].type == "INV" and ref2alt_comps[comp_index + 1].type == "INV") or
                        (ref2alt_comps[comp_index - 1].type == "invDUP" and ref2alt_comps[comp_index + 1].type == "INV") or
                        (ref2alt_comps[comp_index - 1].type == "INV" and ref2alt_comps[comp_index + 1].type == "invDUP")):
                    ref2alt_comps[comp_index].type = "INV"
            except IndexError:
                pass

    # collect and merge continuous variants
    if options.dup_to_ins:
        ref2alt_comps = [comp for comp in handle_continuous_ssv_types(ref2alt_comps) if comp.type not in ["DUP", "invDUP"] and comp.length >= options.min_sv_size]
    else:
        ref2alt_comps = [comp for comp in handle_continuous_ssv_types(ref2alt_comps) if comp.length >= options.min_sv_size]
        
        dup_cnt = len([ssv.type for ssv in ref2alt_comps if 'DUP' in ssv.type])

        # filter massive dups
        if dup_cnt >= 3:
            dup_comp_at_mass = []
            for dup_comp in ref2alt_comps:

                if dup_comp.type in ["DUP", "invDUP"]:
                    dup_comp_at_mass.append(dup_comp)

            for dup_comp in dup_comp_at_mass:
                ref2alt_comps.remove(dup_comp)
        else:

            dup_comp_at_repeats = []
            for dup_comp in ref2alt_comps:

                if dup_comp.type != "DUP":
                    continue

                # determine whether the dup is located at repeat regions
                dup_comp_seq = fetch_ssv_comp_seq(dup_comp, dotplot_ref_start, dotplot_bundle, max_seq_len=None)
                repeat_ratio = calculate_seq_repeat_ratio(dup_comp_seq)

                if repeat_ratio > 1.2:  # the dup locates at a repeat region
                    dup_comp_at_repeats.append(dup_comp)
                    continue

            for dup_comp in dup_comp_at_repeats:
                ref2alt_comps.remove(dup_comp)

    return ref2alt_comps


def deal_with_multi_ins_del(ref2alt_comps, dotplot_ref_start, dotplot_bundle, options):
    ref2alt_comps = sorted(ref2alt_comps, key=lambda x: (x.start, x.end))

    ref2alt_comps_types = [ssv.type for ssv in ref2alt_comps]

    if ("DEL" in ref2alt_comps_types and "INS" in ref2alt_comps_types) or ref2alt_comps_types.count("INS") > 1 or ref2alt_comps_types.count("DEL") > 1:
        allowed_seq_len = 5000

        # deal with DEL INS from unlinear mapping due to snps
        identical_comps = []
        ins_comps = sorted([ssv for ssv in ref2alt_comps if ssv.type in ["INS"]], key=lambda x: x.length, reverse=True)
        del_comps = sorted([ssv for ssv in ref2alt_comps if ssv.type in ["DEL"]], key=lambda x: x.length, reverse=True)

        for ins_comp in ins_comps:
            if ins_comp in identical_comps:
                continue
            
            ins_seq = fetch_ssv_comp_seq(ins_comp, dotplot_ref_start, dotplot_bundle, max_seq_len=allowed_seq_len)

            for del_comp in del_comps:
                if del_comp in identical_comps:
                    continue

                if abs(ins_comp.length - del_comp.length) <= options.min_sv_size:
                    del_seq = fetch_ssv_comp_seq(del_comp, dotplot_ref_start, dotplot_bundle, max_seq_len=allowed_seq_len)

                    seq_similarity = calculate_seq_similarity(ins_seq, del_seq)

                    if seq_similarity >= 0.5:
                        identical_comps.append(del_comp)
                        identical_comps.append(ins_comp)
                        break

        for identical_comp in identical_comps:
            ref2alt_comps.remove(identical_comp)

        # choose the largest comp as a anchor, compare the others with it
        del_ins_comps = sorted([ssv for ssv in ref2alt_comps if ssv.type in ["DEL", "INS"]], key=lambda x: x.length, reverse=True)

        if len(del_ins_comps) != 0:
            largest_del_ins_comp = del_ins_comps[0]
            identical_comps = [largest_del_ins_comp]

            # compare with the largest del or ins comp, and calculate the sequence identity
            # fetch the sequence for largest del_ins_comp
            largest_seq = fetch_ssv_comp_seq(largest_del_ins_comp, dotplot_ref_start, dotplot_bundle, max_seq_len=allowed_seq_len)

            for index in range(1, len(del_ins_comps)):
                cur_del_ins_comp = del_ins_comps[index]
                cur_seq = fetch_ssv_comp_seq(cur_del_ins_comp, dotplot_ref_start, dotplot_bundle, max_seq_len=allowed_seq_len)

                seq_similarity = calculate_seq_similarity(largest_seq, cur_seq)

                if seq_similarity >= 0.5:
                    identical_comps.append(cur_del_ins_comp)

            # find the final type and length
            final_del_ins_chr = identical_comps[0].chr
            final_del_ins_pos = identical_comps[0].start  # use the leftmost cord
            final_del_ins_length = 0

            for identical_comp in identical_comps:
                if identical_comp.type == "DEL":
                    final_del_ins_length -= identical_comp.length
                else:
                    final_del_ins_length += identical_comp.length

                ref2alt_comps.remove(identical_comp)

            if final_del_ins_length >= options.min_sv_size:  # insertion
                ref2alt_comps.append(SSV("INS", final_del_ins_chr, final_del_ins_pos, final_del_ins_pos, alt_seq="N" * final_del_ins_length))
            elif final_del_ins_length <= -options.min_sv_size:  # deletion
                final_del_ins_length = abs(final_del_ins_length)
                ref2alt_comps.append(SSV("DEL", final_del_ins_chr, final_del_ins_pos, final_del_ins_pos + final_del_ins_length))
            else:
                pass

    return ref2alt_comps


def manual_interpret_inv(segment, raw_prediction):
    # we calculate a manual type
    x2x_value = float(segment[2])
    x2y_value = float(segment[3])
    x2y_rev_value = float(segment[4])
    allowed_thresh = x2x_value * 0.1

    if x2x_value - allowed_thresh <= abs(x2y_rev_value) <= x2x_value + allowed_thresh:
        manual_prediction = label2index["INV"]
        return manual_prediction

    else:
        return raw_prediction