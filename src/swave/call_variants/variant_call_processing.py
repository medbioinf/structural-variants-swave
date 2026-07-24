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
import gzip
import pickle
import logging

from .structures import SSV
from .variant_calling import prediction_to_ssv_comps, merge_ref2alt_alt2ref_csv, manual_interpret_inv

logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def process_predictions_to_tsv(predictions_pkl_path, projections_pkl_path, dotplots_pkl_path, output_tsv_path, options):
    """
    Processes the predictions and dotplot bundles to generate a TSV file of called variants.
    """
    logging.info(f"Loading predictions from {predictions_pkl_path}")
    with gzip.open(predictions_pkl_path, 'rb') as pkl_file:
        raw_prediction_dict = pickle.load(pkl_file)
    
    logging.info(f"Loading projection matrices from {projections_pkl_path}")
    with gzip.open(projections_pkl_path, 'rb') as pkl_file:
        projections_dict = pickle.load(pkl_file)
    
    logging.info(f"Loading dotplot bundles from {dotplots_pkl_path}")
    with gzip.open(dotplots_pkl_path, 'rb') as pkl_file:
        dotplot_dict = pickle.load(pkl_file)    
    
    grouped_predictions = {}
    for full_key, prediction in raw_prediction_dict.items():
        if full_key.endswith("|ref2alt") or full_key.endswith("|alt2ref"):
            base_dotplot_id, dotplot_type = full_key.rsplit("|", 1)
        else:
            base_dotplot_id, dotplot_type = full_key, "ref2alt"
            
        if base_dotplot_id not in grouped_predictions:
            grouped_predictions[base_dotplot_id] = {"ref2alt": None, "alt2ref": None}
        grouped_predictions[base_dotplot_id][dotplot_type] = prediction

    with open(output_tsv_path, 'w') as tsv_out:
        tsv_out.write("SAMPLE_ID\tSNARL_ID\tSNARL_REF_START\tSNARL_REF_END\tCHROM\tDOTPLOT_REF_START"
                      "\tDOTPLOT_REF_END\tREVERSED_MAPPING\tSV_TYPE\tSV_LENGTH\tSCORE\tBKPS\n")
        
        called_count = 0
        total_processed_snarls = len(grouped_predictions)
        
        for dotplot_id, directions in grouped_predictions.items():
            if dotplot_id not in dotplot_dict or dotplot_id not in projections_dict:
                logging.warning(f"Dotplot bundle {dotplot_id} not found. Skipping.")
                continue
                
            dotplot_bundle = dotplot_dict[dotplot_id]
            proj_bundle = projections_dict[dotplot_id]
            dotplot_stride_size = dotplot_bundle["stride_size"]
            
            id_parts = dotplot_id.split("|")
            snarl_id_clean = id_parts[0]
            snarl_ref_start = int(id_parts[1])
            snarl_ref_end = int(id_parts[2])
            chrom = id_parts[3]
            dotplot_ref_start = id_parts[4]
            dotplot_ref_end = id_parts[5]
            reversed_mapping = id_parts[6].replace("rev_", "")
            
            comps = {"ref2alt": [], "alt2ref": []}
            
            for d_type in ["ref2alt", "alt2ref"]:
                pred = directions[d_type]
                if pred is None:
                    continue
                
                matrix = proj_bundle[d_type]
                
                for index in range(len(pred)):
                    if index == 0 or index == len(pred) - 1:
                        continue
                    if pred[index] == 1:
                        cur_segment = matrix[index]
                        allowed_thresh = float(cur_segment[2]) * 0.1
                        if pred[index - 1] != 1:
                            prev_seg = matrix[index - 1]
                            if abs(float(cur_segment[3]) - float(prev_seg[3])) <= allowed_thresh and abs(float(cur_segment[4]) - float(prev_seg[4])) <= allowed_thresh:
                                pred[index] = pred[index - 1]
                        if pred[index + 1] != 1:
                            latt_seg = matrix[index + 1]
                            if abs(float(cur_segment[3]) - float(latt_seg[3])) <= allowed_thresh and abs(float(cur_segment[4]) - float(latt_seg[4])) <= allowed_thresh:
                                pred[index] = pred[index + 1]
                
                comps[d_type] = prediction_to_ssv_comps(
                    matrix=matrix,
                    prediction=pred,
                    dotplot_ref_chr=chrom,
                    dotplot_ref_start=snarl_ref_start,
                    dotplot_stride_size=dotplot_stride_size,
                    options=options
                )

            if len(comps["ref2alt"]) == 0 and directions["ref2alt"] is not None:
                matrix = proj_bundle["ref2alt"]
                for i in range(len(directions["ref2alt"])):
                    directions["ref2alt"][i] = manual_interpret_inv(matrix[i], directions["ref2alt"][i])

                comps["ref2alt"] = prediction_to_ssv_comps(
                    matrix=matrix,
                    prediction=directions["ref2alt"],
                    dotplot_ref_chr=chrom,
                    dotplot_ref_start=snarl_ref_start,
                    dotplot_stride_size=dotplot_stride_size,
                    options=options
                )

            # fallback if ref2alt failed the prediction, then we retrieve ssvs from alt2ref
            if len(comps["ref2alt"]) == 0 and len(comps["alt2ref"]) != 0:
                for comp in comps["alt2ref"]:
                    if comp.type == "INV":
                        retrieved_ssv = SSV(
                            "INV", 
                            comp.chr, 
                            snarl_ref_start + (snarl_ref_end - comp.end),
                            snarl_ref_start + (snarl_ref_end - comp.start)
                        )
                        comps["ref2alt"] = [retrieved_ssv]
                        break

            csv_obj = merge_ref2alt_alt2ref_csv(
                ref2alt_comps=comps["ref2alt"],
                alt2ref_comps=comps["alt2ref"],
                dotplot_id=dotplot_id,
                dotplot_bundle=dotplot_bundle,
                options=options
            )
            
            if csv_obj and len(csv_obj.ssv_components) > 0:
                bkp_list = [
                    f"{ssv.type}_{ssv.length}_{ssv.source_chr.replace(':', '-').replace('_', '-')}_{ssv.source_start}_{ssv.source_end}_{ssv.insert_chr.replace(':', '-').replace('_', '-')}_{ssv.insert_start}_{ssv.insert_end}"
                    if ssv.type in ["DUP", "invDUP"]
                    else f"{ssv.type}_{ssv.length}_{ssv.chr.replace(':', '-').replace('_', '-')}_{ssv.start}_{ssv.end}"
                    for ssv in csv_obj.ssv_components
                ]
                bkp_string = ",".join(bkp_list)
                
                tsv_out.write(
                    f"{options.sample_id}\t{snarl_id_clean}\t{snarl_ref_start}\t{snarl_ref_end}\t"
                    f"{chrom}\t{dotplot_ref_start}\t{dotplot_ref_end}\t{reversed_mapping}\t"
                    f"{csv_obj.type}\t{csv_obj.length}\t{csv_obj.score}\t{bkp_string}\n"
                )
                called_count += 1
    
    if called_count == 0:
        logging.info(f"Processed {total_processed_snarls} snarls: No variants found (all predictions were REF/None).")
    else:
        logging.info(f"Processed {total_processed_snarls} snarls: Successfully wrote {called_count} called variants to TSV.")
