#! usr/bin/env python

"""
Runs Pagel on matrix derived from kmerkit
"""

#I'll start by making this separately, but I think these can be functions
#of my Pagel class object

import numpy as np
import pandas as pd #assuming matrix will be a pandas df
from hogtie import binary_state_model

DATA = my_input_file

class MatrixParser:
    """
    Runs BinaryStateModel on matrix columns, returns a likelihood score for each column, flags
    likelihood scores that do not meet a certain threshold
    """
    def __init__(self, tree, matrix):
        self.tree = tree
        self.matrix = matrix

    def csv_to_df:
        """
        """
        data = pd.read_csv(DATA)


    def pagel_run(self):
        """
        runs Pagel on the array
        """
        pass

    def total_matrix_run(self):
        """
        iterates over df to get likelihoods, (object.log_lik, where object is an instance
        of BinaryStateModel
        """
        pass

    def threshold(self):
        """
        Pulls out likelihoods that do not meet a minimum likelihood threshold
        """
        pass
