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
import logging
import numpy as np
from matplotlib import pyplot as plt

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

        self.matrix = np.zeros((int(self.seq_y_len / self.stride_size) + 1, int(self.seq_x_len / self.stride_size) + 1))
        self.matrix_rev = np.zeros((int(self.seq_y_len / self.stride_size) + 1, int(self.seq_x_len / self.stride_size) + 1))

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

        self.matrix = self.matrix[:-1, :-1]
        self.matrix_rev = self.matrix_rev[:-1, :-1]

    def create_matrix_against_x(self, given_x_kmer_index=None):
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
                    self.matrix[index_on_y, indexes_on_x] = 1

            # for reversed kmer
            if not self.skip_reverse:
                kmer_str_reversed = reverse_complement_seq(kmer_str)
                indexes_on_x = self.seq_x_kmer_index.find_all(kmer_str_reversed)
                if indexes_on_x is not None:
                    self.matrix[index_on_y, indexes_on_x] = 1
                    self.matrix_rev[index_on_y, indexes_on_x] = 1

            pos_on_y += self.stride_size

    def create_matrix_against_y(self, given_y_kmer_index=None):
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
                    self.matrix[indexes_on_y, index_on_x] = 1

            # for reversed kmer
            if not self.skip_reverse:
                kmer_str_reversed = reverse_complement_seq(kmer_str)
                indexes_on_y = self.seq_y_kmer_index.find_all(kmer_str_reversed)
                if indexes_on_y is not None:
                    self.matrix[indexes_on_y, index_on_x] = 1
                    self.matrix_rev[indexes_on_y, index_on_x] = 1

            pos_on_x += self.stride_size

    def get_seq_x_kmer_index(self):
        return self.seq_x_kmer_index

    def get_seq_y_kmer_index(self):
        return self.seq_y_kmer_index
    
    def rotate_to_alt2ref(self):
        """
        Rotates the dotplot matrices and meta data to switch from ref2alt to alt2ref or vice versa.
        """
        self.matrix = np.fliplr(np.rot90(self.matrix, k=-1))
        self.matrix_rev = np.fliplr(np.rot90(self.matrix_rev, k=-1))

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
        project_x = np.sum(self.matrix, axis=0)
        augment_coeff = int(100 * np.average(project_x))
        
        if augment:
            project_x[np.diag(self.matrix) == 1] += augment_coeff

        return project_x, augment_coeff

    def get_project_x_rev(self, baseline=0):
        """
        Projects the reverse matrix onto the x-axis with a baseline.
        """
        project_x_rev = baseline + np.sum(self.matrix_rev, axis=0)
        return project_x_rev

    def get_project_y(self, augment=False):
        """
        Projects the matrix onto the y-axis (sum of rows).
        """
        project_y = np.sum(self.matrix, axis=1)
        augment_coeff = int(100 * np.average(project_y))

        if augment:
            project_y[np.diag(self.matrix) == 1] += augment_coeff

        return project_y, augment_coeff

    def get_project_y_rev(self, baseline=0):
        """
        Projects the reverse matrix onto the y-axis with a baseline.
        """
        project_y_rev = baseline + np.sum(self.matrix_rev, axis=1)
        return project_y_rev

    def to_png(self, reverse=False, out_img=False):
        self.dotplot_file = self.out_prefix + ".dotplot.png"
        target_matrix = self.matrix_rev if reverse else self.matrix
        
        if target_matrix.size == 0:
            logging.warning(f"Skipping PNG generation for {self.dotplot_file}: Matrix is empty.")
            return
        
        if np.max(target_matrix) == 0:
            matrix_resize_norm = 255 * np.ones(np.shape(target_matrix))
        else:
            matrix_resize_norm = 255 * abs(target_matrix - np.max(target_matrix)) / (np.max(target_matrix) - np.min(target_matrix))

        if out_img:
            plt.imsave(self.dotplot_file, matrix_resize_norm, cmap='gray')
            plt.close()


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
