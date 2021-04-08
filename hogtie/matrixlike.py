#! usr/bin/env python

"""
Runs BinaryStateModel on matrix of binary character state data for gwas data for the input tree.
The matrix should correspond to presence/absence data corresponding to sequence variants (this
could be kmers, snps, transcripts, etc.).

Parameters
    ----------
    tree: newick string or toytree object
        species tree to be used for the GWAS
    matrix: pandas.dataframe object, csv
    model: str
        Either equal rates ('ER') or all rates different ('ARD')
    prior: float
        Prior probability that the root state is 1 (default=0.5). Flat, uniform prior is assumed.
"""


import numpy as np
import toytree
import pandas as pd #assuming matrix will be a pandas df
from loguru import logger
from hogtie.binary_state_model import BinaryStateModel


class MatrixParser:
    """
    Runs BinaryStateModel on matrix columns, returns a likelihood score for each column, flags
    likelihood scores that do not meet a certain threshold

    """
    def __init__(self, 
        tree,               #must be Toytree class object
        matrix = None,      #must be pandas DataFrame class object
        model = None,
        prior = 0.5    
        ):       #must be .csv filetype - as filepath


        if isinstance(tree, toytree.tree):
            self.tree = tree
        elif isinstance(tree, str):
            self.tree = toytree.tree(tree, tree_format=0)
        else: 
            raise Exception('tree must be either a newick string or toytree object')


        if isinstance(matrix, pd.DataFrame):
            self.matrix = matrix  
        else:
            self.matrix = pd.read_csv(matrix)

        self.model = model
        self.prior = prior

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
            out = BinaryStateModel(self.tree, data, self.model, self.prior)
            out.optimize()
        
            lik = out.log_lik
            likelihoods = np.append(likelihoods, lik)
         
            self.likelihoods = pd.DataFrame(likelihoods)

        logger.debug(f'Likelihoods for each column: {self.tree.get_node_values("likelihood",True,True)}')

    
if __name__ == "__main__":
    import toytree
    import os
    HOGTIEDIR = os.path.dirname(os.getcwd())
    tree1 = toytree.rtree.unittree(ntips=10)
    file1 = os.path.join(HOGTIEDIR, "sampledata", "testmatrix.csv")
    testmatrix = MatrixParser(tree=tree1, matrix=file1, model='ARD')
    testmatrix.matrix_likelihoods()
    print(testmatrix.likelihoods)
