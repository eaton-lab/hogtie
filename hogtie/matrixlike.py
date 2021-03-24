#! usr/bin/env python

"""
Runs Pagel on matrix derived from kmerkit
"""

#I'll start by making this separately, but I think these can be functions
#of my Pagel class object

import numpy as np
import pandas as pd #assuming matrix will be a pandas df
from loguru import logger
from hogtie.binary_state_model import BinaryStateModel


#DATA = pd.read_csv(my_input_file)

class MatrixParser:
    """
    Runs BinaryStateModel on matrix columns, returns a likelihood score for each column, flags
    likelihood scores that do not meet a certain threshold
    """
    def __init__(self, 
        tree,               #must be Toytree class object
        file):       #must be .csv filetype - as filepath

        self.tree = tree
        self.matrix = pd.read_csv(file)

        #check that data only contains 1's and 0's
        # ...

    def matrix_likelihoods(self):
        """
        Gets likelihoods for each column of the matrix
        """
        likelihoods = []
        for column in matrix:
            data = matrix[column]
            out = BinaryStateModel(tree, data)
            out.optimize()
            
            likelihoods.append(out.log_lik)
            
        return likelihoods
        logger.info(f"Likelihoods: {likelihoods}")

    def threshold(self):
        """
        Pulls out likelihoods that do not meet a minimum likelihood threshold
        """
        pass

if __name__ == "__main__":
    import toytree
    import os
    HOGTIEDIR = os.path.dirname(os.getcwd())
    tree = toytree.rtree.unittree(ntips=10)
    file = os.path.join(HOGTIEDIR, "sampledata", "testmatrix.csv")
    testmatrix = MatrixParser(tree=tree, file=file)
    testmatrix.column_to_list()
    testmatrix.total_matrix_run()

