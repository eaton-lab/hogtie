#! usr/bin/env python

"""
generating data using ipcoal
"""

import ipcoal
import toytree
from hogtie import BinaryStateModel, MatrixParser

def create_null(tree):
    """
    Simulates SNPs across the input tree to create the null expectation for likelihood scores. 
    Deviation from the null will be flagged. 
    """
    
    #high ILS
    mod = ipcoal.Model(tree=tree, Ne=(tree.treeheight + 1e3))
    mod.sim_loci(nloci=1, nsites=1000)
    vcf = mod.write_vcf()
    null_genos = vcf.iloc[:, 9:].T

    #run the Binary State model on the matrix and get likelihoods
    null = MatrixParser(tree=tree, matrix=null_genos, model='ARD')
    null.matrix_likelihoods()

    #null likelihood expectation
    null.likelihoods

    