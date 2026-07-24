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


class SSV:
    def __init__(self, type, chr, start, end, alt_seq=None):
        self.id = f"{chr}-{start}-{end}-{type}"

        self.type = type

        self.chr = chr
        self.start = start
        self.end = end

        self.alt_seq = alt_seq

        if self.type == "INS" or self.type == "insertion":
            self.length = len(self.alt_seq)
        else:
            self.length = self.end - self.start + 1

        self.source_chr = None
        self.source_start = None
        self.source_end = None

        self.insert_chr = None
        self.insert_start = None
        self.insert_end = None

        self.raw_dotplot = None
        self.raw_chr = None
        self.raw_start = None
        self.raw_end = None

        self.projection_matrix = None

        self.id = None  # TODO: find out why they do this, probably old code, remove?

        self.within_inv = False
        self.within_repeat = False

    def set_id(self, id):
        self.id = id

    def refine_dup_events(self, inserted_chr, inserted_pos):
        """
        Set the source and insert coordinates for duplication events.
        """
        self.source_chr = self.chr
        self.source_start = self.start
        self.source_end = self.end

        self.insert_chr = inserted_chr
        self.insert_start = inserted_pos
        self.insert_end = inserted_pos + 1

    def set_raw_dotplot(self, raw_dotplot):
        self.raw_dotplot = raw_dotplot

    def set_raw_projection(self, projection_matrix):
        self.projection_matrix = projection_matrix

    def set_raw_cords(self, raw_chr, raw_start, raw_end):
        self.raw_chr = raw_chr
        self.raw_start = raw_start
        self.raw_end = raw_end

    def string_format(self):
        return f"{self.type},{self.chr}-{self.start}-{self.end},{self.length}"

    def to_string(self):
        if "DUP" in self.type or "duplication" in self.type:
            return f"{self.type}: {self.source_start}-{self.source_end}-{self.insert_start}-{self.insert_end},{self.length}"

        else:
            return f"{self.type}: {self.start}-{self.end},{self.length}"


class CSV:
    def __init__(self, ssv_components):
        self.ssv_components = ssv_components
        self.update_csv_info()
        self.score = "NA"

    def set_score(self, score):
        if score < 0 or score > 1.5:
            self.score = "LowQual"
        else:
            self.score = "PASS"

    def update_csv_info(self):
        if len(self.ssv_components) == 0:
            self.chr = -1
            self.start = -1
            self.end = -1
            self.length = -1
            self.type = "+".join([ssv.type for ssv in self.ssv_components])
        else:
            self.ssv_components = sorted(self.ssv_components, key=lambda x: (x.start, x.end))

            # generate CSV info
            self.chr = self.ssv_components[0].chr
            self.start = self.ssv_components[0].start
            self.end = self.ssv_components[-1].end
            self.length = int(np.sum([ssv.length for ssv in self.ssv_components]))

            self.type = "+".join([ssv.type for ssv in self.ssv_components])

    def add_ssv(self, new_ssv):
        self.ssv_components.append(new_ssv)
        self.update_csv_info()

    def string_format(self):
        return "+".join([f"{ssv.type},{ssv.chr}-{ssv.start}-{ssv.end},{ssv.length}" for ssv in self.ssv_components])

    def to_string(self):
        return self.chr, self.start, self.end, self.length, self.score, [ssv.to_string() for ssv in self.ssv_components]
