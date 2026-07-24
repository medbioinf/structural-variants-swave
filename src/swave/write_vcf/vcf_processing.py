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
import pysam
import logging

from src import __version__


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def hap_gt_to_sample_gt(sample_gt, mode="raw"):
    if mode == "split":
        return f"{'|'.join(sample_gt)}"
    else:
        gt_flags = []
        sv_types = []
        sv_lengths = []
        sv_bkps = []

        for hap_gt in sample_gt:
            hap_gt_split = hap_gt.split(":")
            gt_flags.append(hap_gt_split[0])
            sv_types.append(hap_gt_split[1])
            sv_lengths.append(hap_gt_split[2])
            sv_bkps.append(hap_gt_split[3])

        return f"{'|'.join(gt_flags)}:{'|'.join(sv_types)}:{'|'.join(sv_lengths)}:{'|'.join(sv_bkps)}"


def calculate_ac_af_an_ns(gt_list, asm_names, sample_names):
    AC = 0
    non_called_cnt = 0
    for gt_flag in gt_list:
        flag = gt_flag.split(":")[0]
        if flag == ".":
            non_called_cnt += 1
        else:
            if flag != "0":
                AC += 1

    AN = len(asm_names) - non_called_cnt
    AF = round(AC / AN, 5) if AN != 0 else 0.0
    NS = len(sample_names)

    return AC, AF, AN, NS


def output_vcf_header(vcf_fout, ref_asm_path, sample_names):
    """
    Generates a standardized VCF header.
    """
    print("##fileformat=VCFv4.3", file=vcf_fout)
    print(f"##source=xxx v{__version__}", file=vcf_fout)

    # add chromosome info
    ref_file = pysam.FastaFile(ref_asm_path)
    chroms = ref_file.references
    for chrom in chroms:
        chr_length = ref_file.get_reference_length(chrom)
        print(f"##contig=<ID={chrom},length={chr_length}>", file=vcf_fout)
    ref_file.close()

    print("##ID=<ID=Path,Description=\"Snarl path\">", file=vcf_fout)
    print("##REF=<ID=Path,Description=\"Reference path\">", file=vcf_fout)
    print("##ALT=<ID=Path,Description=\"Alternative paths\">", file=vcf_fout)

    print("##FILTER=<ID=PASS,Description=\"Passed\">", file=vcf_fout)
    print("##FILTER=<ID=LowQual,Description=\"LowQual\">", file=vcf_fout)
    print("##FILTER=<ID=MediumQual,Description=\"MediumQual\">", file=vcf_fout)

    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV/CSV\">", file=vcf_fout)
    print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV/CSV\">", file=vcf_fout)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV/CSV\">", file=vcf_fout)
    print("##INFO=<ID=BKPS,Number=.,Type=String,Description=\"Breakpoints of the SV/CSV (length-start-end-insert)\">", file=vcf_fout)

    print("##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">", file=vcf_fout)
    print("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">", file=vcf_fout)
    print("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">", file=vcf_fout)
    print("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">", file=vcf_fout)
    print("##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">", file=vcf_fout)

    # add gt info
    print("##FORMAT=<ID=GT,Number=.,Type=String,Description=\"SV genotype for each haplotype\">", file=vcf_fout)
    print("##FORMAT=<ID=TYPE,Number=.,Type=String,Description=\"SV type for each haplotype\">", file=vcf_fout)
    print("##FORMAT=<ID=LENGTH,Number=.,Type=String,Description=\"SV length for each haplotype\">", file=vcf_fout)
    print("##FORMAT=<ID=QUAL,Number=.,Type=String,Description=\"SV quality for each haplotype\">", file=vcf_fout)
    print("##FORMAT=<ID=BKPS,Number=.,Type=String,Description=\"SV breakpoint for each haplotype\">", file=vcf_fout)

    alt_asms = "\t".join(sample_names)
    print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{alt_asms}", file=vcf_fout)


def output_vcf_records(vcf_fout, records):
    # sort records by chrom and start positions
    records = sorted(records, key=lambda x: (x[0], int(x[1])))
    for record in records:
        print("\t".join(record), file=vcf_fout)


def output_at_snarl_level_tsv(vcf_records, vcf_records_split, snarl_id, snarl_coords, snarl_sample_calls, asm_names):
    if snarl_coords["reversed_mapping"]:
        snarl_gt_list = ["0:.:.:." for _ in range(len(asm_names))]
    else:
        snarl_gt_list = [".:.:.:." for _ in range(len(asm_names))]
    
    alt_variants = []
    for call in snarl_sample_calls.values():
        sv_type = call["type"]
        sv_length = call["length"]
        sv_bkps = call["bkps"]
        
        if sv_type == "REF" or not sv_type or sv_type == ".":
                continue
        
        variant_tuple = (sv_type, sv_length, sv_bkps)
        if variant_tuple not in alt_variants:
            alt_variants.append(variant_tuple)
    
    snarl_include_alt_variants = []
    for i, asm_name in enumerate(asm_names):
        if asm_name in snarl_sample_calls:
            call = snarl_sample_calls[asm_name]
            sv_type = call["type"]
            sv_length = call["length"]
            sv_bkps = call["bkps"]
            
            if sv_type == "REF" or not sv_type or sv_type == ".":
                snarl_gt_list[i] = "0:.:.:."
            else:
                variant_idx = alt_variants.index((sv_type, sv_length, sv_bkps)) + 1
                snarl_gt_list[i] = f"{variant_idx}:{sv_type}:{sv_length}:{sv_bkps}"
                
                alt_tag = f"<{sv_type}>"
                if alt_tag not in snarl_include_alt_variants:
                    snarl_include_alt_variants.append(alt_tag)
        else:
            snarl_gt_list[i] = "0:.:.:."
                
        
    sample_mapping = {} 
    for hap in asm_names:
        sample_name = hap.rsplit("_", 1)[0] if "_" in hap else hap
        if sample_name not in sample_mapping:
            sample_mapping[sample_name] = []
        sample_mapping[sample_name].append(hap)
        
    sample_names = sorted(list(sample_mapping.keys()))
    
    snarl_AC, snarl_AF, snarl_AN, snarl_NS = calculate_ac_af_an_ns(
        snarl_gt_list, 
        asm_names=asm_names, 
        sample_names=sample_names
    )
    snarl_NS = len(sample_names)
    
    diploid_sample_gts = []
    for sample in sample_names:
        sample_haps = sample_mapping[sample]
        sample_hap_gts = [snarl_gt_list[asm_names.index(hap)] for hap in sample_haps]
        
        # generates diploid genotypes (e.g., "1|0:INV|.:1500|.:...")
        diploid_sample_gts.append(hap_gt_to_sample_gt(sample_hap_gts, mode="raw"))
    
    snarl_qual = "MediumQual" if snarl_coords["reversed_mapping"] else "PASS"
    info_string = f"END={snarl_coords['end']};AC={snarl_AC};AF={snarl_AF};AN={snarl_AN};NS={snarl_NS};LV=0"
    
    vcf_records.append([
        str(snarl_coords["chrom"]),
        str(snarl_coords["start"]),
        str(snarl_id),  # ID
        str(snarl_id),  # REF
        ",".join(snarl_include_alt_variants) if snarl_include_alt_variants else ".",    # ALT
        ".",    # QUAL
        snarl_qual, # FILTER
        info_string,    # INFO
        "GT:TYPE:LENGTH:BKPS",  # FORMAT
        "\t".join(diploid_sample_gts)  # multi-sample GTs
    ])
    
    # handle the split VCF records for each alternative variant
    if len(alt_variants) <= 1:
        split_diploid_gts = []
        for sample in sample_names:
            sample_haps = sample_mapping[sample]
            sample_hap_gts = [snarl_gt_list[asm_names.index(hap)].split(":")[0] for hap in sample_haps]
            split_diploid_gts.append(hap_gt_to_sample_gt(sample_hap_gts, mode="split"))
        
        single_info = info_string
        if len(alt_variants) == 1:
            single_info = f"END={snarl_coords['end']};SVLEN={alt_variants[0][1]};SVTYPE={alt_variants[0][0]};AC={snarl_AC};AF={snarl_AF};AN={snarl_AN};NS={snarl_NS};LV=0"
        
        vcf_records_split.append([
            str(snarl_coords["chrom"]),
            str(snarl_coords["start"]),
            str(snarl_id),
            str(snarl_id),
            ",".join(snarl_include_alt_variants) if snarl_include_alt_variants else ".",
            ".",
            snarl_qual,
            single_info,
            "GT",
            "\t".join(split_diploid_gts)
        ])
    else:
        for alt_idx, (sv_type, sv_length, sv_bkps) in enumerate(alt_variants):
            split_gt_list = []
            for gt_info in snarl_gt_list:
                gt_flag = gt_info.split(":")[0]
                if gt_flag == str(alt_idx + 1):
                    split_gt_list.append("1")
                elif gt_flag == ".":
                    split_gt_list.append(".")
                else:
                    split_gt_list.append("0")
            
            split_AC, split_AF, split_AN, split_NS = calculate_ac_af_an_ns(split_gt_list, asm_names, sample_names)
            
            split_diploid_gts = []
            for sample in sample_names:
                sample_haps = sample_mapping[sample]
                sample_hap_gts = [split_gt_list[asm_names.index(hap)] for hap in sample_haps]
                split_diploid_gts.append(hap_gt_to_sample_gt(sample_hap_gts, mode="split"))
                
            split_info = f"END={snarl_coords['end']};SVLEN={sv_length};SVTYPE={sv_type};BKPS={sv_bkps};AC={split_AC};AF={split_AF};AN={split_AN};NS={split_NS};LV=0"
            
            vcf_records_split.append([
                str(snarl_coords["chrom"]),
                str(snarl_coords["start"]),
                str(snarl_id),
                str(snarl_id),
                f"<{sv_type}>",
                ".",
                snarl_qual,
                split_info,
                "GT",
                "\t".join(split_diploid_gts)
            ])


def generate_vcf_records_from_tsv(interpret_res_dict, snarl_coords_dict, asm_names):
    vcf_records = []
    vcf_records_split = []

    for snarl_id in interpret_res_dict:
        snarl_coords = snarl_coords_dict[snarl_id]
        snarl_sample_calls = interpret_res_dict[snarl_id]
        
        output_at_snarl_level_tsv(
            vcf_records=vcf_records,
            vcf_records_split=vcf_records_split,
            snarl_id=snarl_id,
            snarl_coords=snarl_coords,
            snarl_sample_calls=snarl_sample_calls,
            asm_names=asm_names
        )

    return vcf_records, vcf_records_split
    

def process_tsv_files_to_vcf(tsv_files, output_vcf_path, ref_asm_path):
    interpret_res_dict = {}
    snarl_coords_dict = {}
    found_samples = set()
    
    for tsv_path in tsv_files:
        with open(tsv_path, 'r') as tsv_file:
            header = tsv_file.readline().strip().split('\t')
            header_idx = {col: i for i, col in enumerate(header)}
            
            for line in tsv_file:
                parts = line.strip().split('\t')
                if not parts or len(parts) < len(header):
                    continue
                
                sample_id = parts[header_idx['SAMPLE_ID']]
                snarl_id = parts[header_idx['SNARL_ID']]
                sv_type = parts[header_idx['SV_TYPE']]
                
                found_samples.add(sample_id)
                
                if snarl_id not in interpret_res_dict:
                    interpret_res_dict[snarl_id] = {}
                    snarl_coords_dict[snarl_id] = {
                        "chrom": parts[header_idx["CHROM"]],
                        "start": int(parts[header_idx["SNARL_REF_START"]]),
                        "end": int(parts[header_idx["SNARL_REF_END"]]),
                        "reversed_mapping": parts[header_idx["REVERSED_MAPPING"]].lower() == "true"
                    }
                
                interpret_res_dict[snarl_id][sample_id] = {
                    "type": sv_type,
                    "length": parts[header_idx["SV_LENGTH"]],
                    "score": parts[header_idx["SCORE"]],
                    "bkps": parts[header_idx["BKPS"]]
                }
                
    asm_names = sorted(list(found_samples))
    
    sample_names = sorted(set([asm.rsplit("_", 1)[0] if "_" in asm else asm for asm in asm_names]))
    
    vcf_records, vcf_records_split = generate_vcf_records_from_tsv(interpret_res_dict, snarl_coords_dict, asm_names=asm_names)
    
    with open(output_vcf_path, "w") as fout:
        output_vcf_header(fout, ref_asm_path, sample_names)
        output_vcf_records(fout, vcf_records)
        
    output_split_path = output_vcf_path.replace(".vcf", ".split.vcf")
    with open(output_split_path, "w") as fout_split:
        output_vcf_header(fout_split, ref_asm_path, sample_names)
        output_vcf_records(fout_split, vcf_records_split)
