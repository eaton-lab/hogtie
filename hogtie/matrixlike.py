#! usr/bin/env python

"""
Runs Pagel on matrix derived from kmerkit
"""

#I'll start by making this separately, but I think these can be functions
#of my Pagel class object

import numpy as np
import pandas as pd #assuming matrix will be a pandas df
from loguru import logger
from hogtie import binary_state_model


#DATA = pd.read_csv(my_input_file)

class MatrixParser:
    """
    Runs BinaryStateModel on matrix columns, returns a likelihood score for each column, flags
    likelihood scores that do not meet a certain threshold
    """
    def __init__(self, 
        tree,               #must be Toytree class object
        matrix = None,      #must be pandas DataFrame class object
        file = None):       #must be .csv filetype - as filepath
        

        self.tree = tree    

        if isinstance(matrix, pd.DataFrame):
            self.matrix = matrix  

        elif file != None:
            self.matrix = pd.read_csv(file)

        else:
            raise ValueError("please supply matrix")

    def column_to_list(self):
        """
        takes column in matrix and returns a list
        """
        datalists = []
        for i in range(0, len(self.matrix.columns)):
            data = []
            for index, row in self.matrix.iterrows():
                data.append(row[i])
            datalists.append(data)
        self.datalists=datalists

    def pagel_run(self):
        """
        runs Pagel on the list
        """
        ## - the way I wrote the code, this function is unnecessary ~ elissa
        pass

    def total_matrix_run(self):
        """
        uses for loop to run the two previous functions on the whole matrix. Stores
        a dict of likelihoods, key corresponds to column idx from matrix and value
        is the likelihood
        """

        likelihoods = []
        for item in self.datalists:
            data = item
            out = binary_state_model.BinaryStateModel(tree, data)
            out.optimize()
            
            likelihoods.append(out.log_lik)
        
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

