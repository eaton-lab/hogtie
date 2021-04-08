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
    
    #simulate snps across the input tree
    #note: snp sims take a lot of computational power

    mod = ipcoal.Model(tree=tree, Ne=1e6)
    mod.sim_snps(nsnps=10000)
    null_genos = mod.write_vcf()

    #run the Binary State model on the matrix and get likelihoods
    null = MatrixParser(tree=tree, matrix=null_genos, model='ARD')
    null.matrix_likelihoods()

    #null likelihood expectation
    null.likelihoods

    