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

        #for i in self.matrix:
        #  if i != 1 or 0:
        #        raise ValueError('Only valid trait values are 0 and 1')

    def matrix_likelihoods(self):
        """
        Gets likelihoods for each column of the matrix
        """
        likelihoods = np.empty((0,len(self.matrix.columns)),float)
        for column in self.matrix:
            data = self.matrix[column]
            out = BinaryStateModel(self.tree, data)
            out.optimize()
        
            lik = out.log_lik
            likelihoods = np.append(likelihoods, lik)
         
            self.likelihoods = likelihoods

    def threshold(self):
        """
        Pulls out likelihoods that do not meet a minimum likelihood threshold
        """
        pass

if __name__ == "__main__":
    import toytree
    import os
    HOGTIEDIR = os.path.dirname(os.getcwd())
    tree1 = toytree.rtree.unittree(ntips=10)
    file1 = os.path.join(HOGTIEDIR, "sampledata", "testmatrix.csv")
    testmatrix = MatrixParser(tree=tree1, file=file1)
    testmatrix.matrix_likelihoods()
    print(testmatrix.likelihoods)
