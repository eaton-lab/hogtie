#! usr/bin/env python

"""
Runs Pagel on matrix derived from kmerkit
"""

#I'll start by making this separately, but I think these can be functions
#of my Pagel class object

import numpy as np
import pandas as pd #assuming matrix will be a pandas df
from hogtie import Pagel

class MatrixParser:
    """
    Runs Pagel on matrix columns, returns a likelihood score for each column, flags
    likelihood scores that do not meet a certain threshold
    """
    def __init__(self, tree, matrix):
        self.tree = tree
        self.matrix = matrix

    def column_to_list(self):
        """
        takes column in matrix and returns a list
        """
        pass

    def pagel_run(self):
        """
        runs Pagel on the list
        """
        pass

    def total_matrix_run(self):
        """
        uses for loop to run the two previous functions on the whole matrix. Stores
        a dict of likelihoods, key corresponds to column idx from matrix and value
        is the likelihood
        """
        pass

    def threshold(self):
        """
        Pulls out likelihoods that do not meet a minimum likelihood threshold
        """
        pass
