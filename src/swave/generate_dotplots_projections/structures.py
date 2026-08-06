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

import os
import sys
import logging
import numpy as np
import scipy.sparse as sp
from PIL import Image

from swave.utils.seq_utils import reverse_complement_seq, calculate_stride_size


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


class Dotplot:
    def __init__(self, seq_x, seq_y, kmer_size, out_prefix, against="auto", stride_size=None, given_x_kmer_index=None, given_y_kmer_index=None, skip_forward=False, skip_reverse=False):

        self.seq_x = seq_x.upper()
        self.seq_y = seq_y.upper()

        self.seq_x_len = len(seq_x)
        self.seq_y_len = len(seq_y)

        self.kmer_size = kmer_size
        
        self.skip_forward = skip_forward
        self.skip_reverse = skip_reverse

        if stride_size is None:
            self.stride_size = calculate_stride_size(self.seq_x, self.seq_y)
        else:
            self.stride_size = stride_size

        self.out_prefix = out_prefix

        n_rows = int(self.seq_y_len / self.stride_size) + 1
        n_cols = int(self.seq_x_len / self.stride_size) + 1

        self._rows, self._cols = [], []
        self._rows_rev, self._cols_rev = [], []

        if against == "auto":
            if self.seq_x_len >= self.seq_y_len:
                self.create_matrix_against_x(given_x_kmer_index)
            else:
                self.create_matrix_against_y(given_y_kmer_index)

        elif against == "x":
            self.create_matrix_against_x(given_x_kmer_index)

        elif against == "y":
            self.create_matrix_against_y(given_y_kmer_index)

        else:
            logging.error("No such against axis: {}. Choose from [auto, x, y]".format(against))

        data = np.ones(len(self._rows), dtype=np.uint8)
        self.matrix = sp.csr_matrix((data, (self._rows, self._cols)), shape=(n_rows, n_cols))[:-1, :-1]
        self.matrix.data = np.minimum(self.matrix.data, 1)

        data_rev = np.ones(len(self._rows_rev), dtype=np.uint8)
        self.matrix_rev = sp.csr_matrix((data_rev, (self._rows_rev, self._cols_rev)), shape=(n_rows, n_cols))[:-1, :-1]
        self.matrix_rev.data = np.minimum(self.matrix_rev.data, 1)
        
        del self._rows, self._cols, self._rows_rev, self._cols_rev

    def create_matrix_against_x(self, given_x_kmer_index=None):
        """
        Creates the dotplot matrix by iterating over the y sequence and finding
        all matching k-mers in the x sequence.
        """
        if given_x_kmer_index is None:
            self.seq_x_kmer_index = KmerIndex(self.seq_x, kmer_size=self.kmer_size, stride_size=self.stride_size)
        else:
            self.seq_x_kmer_index = given_x_kmer_index

        pos_on_y = 0
        while pos_on_y < self.seq_y_len:
            kmer_str = self.seq_y[pos_on_y: pos_on_y + self.kmer_size]
            index_on_y = int(pos_on_y / self.stride_size)

            # for original kmer
            if not self.skip_forward:
                indexes_on_x = self.seq_x_kmer_index.find_all(kmer_str)
                if indexes_on_x is not None:
                    for idx_x in indexes_on_x:
                        self._rows.append(index_on_y)
                        self._cols.append(idx_x)

            # for reversed kmer
            if not self.skip_reverse:
                kmer_str_reversed = reverse_complement_seq(kmer_str)
                indexes_on_x = self.seq_x_kmer_index.find_all(kmer_str_reversed)
                if indexes_on_x is not None:
                    for idx_x in indexes_on_x:
                        self._rows.append(index_on_y)
                        self._cols.append(idx_x)
                        self._rows_rev.append(index_on_y)
                        self._cols_rev.append(idx_x)

            pos_on_y += self.stride_size

    def create_matrix_against_y(self, given_y_kmer_index=None):
        """
        Creates the dotplot matrix by iterating over the x sequence and finding
        all matching k-mers in the y sequence.
        """
        if given_y_kmer_index is None:
            self.seq_y_kmer_index = KmerIndex(self.seq_y, kmer_size=self.kmer_size, stride_size=self.stride_size)
        else:
            self.seq_y_kmer_index = given_y_kmer_index

        pos_on_x = 0
        while pos_on_x < self.seq_x_len:
            kmer_str = self.seq_x[pos_on_x: pos_on_x + self.kmer_size]
            index_on_x = int(pos_on_x / self.stride_size)
            
            # for original kmer
            if not self.skip_forward:
                indexes_on_y = self.seq_y_kmer_index.find_all(kmer_str)
                if indexes_on_y is not None:
                    for idx_y in indexes_on_y:
                        self._rows.append(idx_y)
                        self._cols.append(index_on_x)

            # for reversed kmer
            if not self.skip_reverse:
                kmer_str_reversed = reverse_complement_seq(kmer_str)
                indexes_on_y = self.seq_y_kmer_index.find_all(kmer_str_reversed)
                if indexes_on_y is not None:
                    for idx_y in indexes_on_y:
                        self._rows.append(idx_y)
                        self._cols.append(index_on_x)
                        self._rows_rev.append(idx_y)
                        self._cols_rev.append(index_on_x)

            pos_on_x += self.stride_size

    def get_seq_x_kmer_index(self):
        return self.seq_x_kmer_index

    def get_seq_y_kmer_index(self):
        return self.seq_y_kmer_index
    
    def rotate_to_alt2ref(self):
        """
        Rotates the dotplot matrices and meta data to switch from ref2alt to alt2ref or vice versa.
        """
        self.matrix = self.matrix.T[::-1, :]
        self.matrix_rev = self.matrix_rev.T[::-1, :]

        tmp_seq_x = self.seq_x
        self.seq_x = self.seq_y
        self.seq_y = tmp_seq_x

        tmp_seq_x_len = self.seq_x_len
        self.seq_x_len = self.seq_y_len
        self.seq_y_len = tmp_seq_x_len

    def get_project_x(self, augment=False):
        """
        Projects the matrix ont the x-axis (sum of columns).
        """
        project_x = self.matrix.sum(axis=0).A1.astype(np.float64)
        augment_coeff = int(100 * np.average(project_x)) if len(project_x) > 0 else 0
        
        if augment:
            diag = self.matrix.diagonal()
            project_x[diag == 1] += augment_coeff

        return project_x, augment_coeff

    def get_project_x_rev(self, baseline=0):
        """
        Projects the reverse matrix onto the x-axis with a baseline.
        """
        project_x_rev = baseline + self.matrix_rev.sum(axis=0).A1.astype(np.float64)
        return project_x_rev

    def get_project_y(self, augment=False):
        """
        Projects the matrix onto the y-axis (sum of rows).
        """
        project_y = self.matrix.sum(axis=1).A1.astype(np.float64)
        augment_coeff = int(100 * np.average(project_y)) if len(project_y) > 0 else 0

        if augment:
            diag = self.matrix.diagonal()
            project_y[diag == 1] += augment_coeff

        return project_y, augment_coeff

    def get_project_y_rev(self, baseline=0):
        """
        Projects the reverse matrix onto the y-axis with a baseline.
        """
        project_y_rev = baseline + self.matrix_rev.sum(axis=1).A1.astype(np.float64)
        return project_y_rev

    def to_png(self, reverse=False, out_img=False):
        """Saves the dotplot matrix as a PNG image."""
        if not out_img:
            return
        
        self.dotplot_file = self.out_prefix + ".dotplot.png"
        target_matrix = self.matrix_rev if reverse else self.matrix
        
        m_h, m_w = getattr(target_matrix, "shape", (0, 0))
        if m_h == 0 or m_w == 0:
            logging.warning(f"Skipping PNG generation for {self.dotplot_file}: Matrix is empty.")
            return
        
        if hasattr(target_matrix, "tocoo"):
            img_array = np.full((m_h, m_w), 255, dtype=np.uint8)
            coo = target_matrix.tocoo()
            img_array[coo.row, coo.col] = 0
        else:
            img_array = np.where(target_matrix > 0, 0, 255).astype(np.uint8)

        os.makedirs(os.path.dirname(self.dotplot_file), exist_ok=True)
        img = Image.fromarray(img_array, mode='L')
        img.save(self.dotplot_file)


class KmerIndex:
    def __init__(self, seq, kmer_size, stride_size):
        self.seq = seq
        self.index_table = {}

        pos_on_seq = 0
        while pos_on_seq < len(seq):
            index_on_seq = int(pos_on_seq / stride_size)

            for i in range(stride_size):
                kmer_str = seq[pos_on_seq + i: pos_on_seq + i + kmer_size]

                if kmer_str not in self.index_table:
                    self.index_table[kmer_str] = []

                if not self.index_table[kmer_str] or self.index_table[kmer_str][-1] != index_on_seq:
                    self.index_table[kmer_str].append(index_on_seq)

            pos_on_seq += stride_size

    def find_all(self, seq_str):
        """
        Returns all matrix indices where this k-mer exists.
        """
        return self.index_table.get(seq_str)


class Diag:
    def __init__(self, x_start, x_end, y_start, y_end, orient, offset):
        self.x_start = x_start
        self.x_end = x_end
        self.y_start = y_start
        self.y_end = y_end

        self.orient = orient
        self.offset = offset

        self.true_reverse = False

    def __eq__(self, other):
        if not isinstance(other, Diag):
            return False
        
        return self.x_start == other.x_start and self.x_end == other.x_end and self.y_start == other.y_start and self.y_end == other.y_end and self.orient == other.orient

    def to_string(self):
        return f"{self.x_start}-{self.x_end}, {self.y_start}-{self.y_end}, {self.orient}, {self.offset}"
