#! usr/bin/env python

"""
generating null data and testing observations against it
"""

import ipcoal
import toytree
import toyplot
import numpy as np
from hogtie import MatrixParser

class SimulateNull():
    """
    Compare data to expectations

    TO DO: Not sure standard deviation-based comparison is the best for
    log-likelihoods from a statistical point of view. Maybe something more like
    an AIC-type comparison or a likelihood-ratio test?
    """
    def __init__(self, tree, matrix, model=None, prior=0.5):
        self.model = model
        self.prior = prior
        self.matrix = matrix

        if isinstance(tree, toytree.tree):
            self.tree = tree
        elif isinstance(tree, str):
            self.tree = toytree.tree(tree, tree_format=0)
        else: 
            raise Exception('tree must be either a newick string or toytree object')

        self.treeheight = float(self.tree.treenode.height)

    #@property
    def null(self):
        """
        Simulates SNPs across the input tree to create the null expectation for likelihood
        scores and compares
        """
        
        #high ILS
        mod = ipcoal.Model(tree=self.tree, Ne=(self.treeheight + 1e3))
        mod.sim_snps(nsnps=1)
        null_genos = mod.write_vcf().iloc[:, 9:].T

        #run Binary State model on the matrix and get likelihoods
        null = MatrixParser(tree=self.tree, matrix=null_genos, model=self.model)
        null.matrix_likelihoods()

        #get z -scores null likelihood expectations
        std = null.likelihoods.std()
        mean = null.likelihoods.mean()

        #return mean
        #eturn std
        
    # def significance_test(self):
        #"""
        #Identifies kmer phylogenetic patterns that deviate from the null distribution
        #"""
        #get the likelihood value that corresponds 2 standard deviations above the null mean
        high_lik = mean + 2*std

        lik_calc = MatrixParser(tree=self.tree,
                               model=self.model,
                               prior=self.prior,
                               matrix=self.matrix
                               )

        lik_calc.matrix_likelihoods()

        devs = [] #would prefer to append to an empty np.array
        for like in list(lik_calc.likelihoods[0]):
            if like >= high_lik:
                devs.append(1)
            else:
                devs.append(0)

        lik_calc.likelihoods['deviation_score'] = np.array(devs)
        self.likelihoods = lik_calc.likelihoods

    def genome_graph(self):
        """
        Graphs rolling average of likelihoods along the linear genome, identifies 
        regions that deviate significantly from null expectations

        TO DO: change color of outliers, integrate into Hogtie class instead of MatrixParser
        class
        """

        self.likelihoods['rollingav']= self.likelihoods[1].rolling(10, win_type='triang').mean()
        
        plot = toyplot.plot(
            self.likelihoods['rollingav'],
            width = 500,
            height=500,
            color = 'blue'
        )

        return plot

if __name__ == "__main__":
    testtree = toytree.rtree.unittree(ntips=10, treeheight=1e5)
    mod1 = ipcoal.Model(tree=testtree, Ne=1e3, admixture_edges=[(3, 8, 0.5, 0.5)], nsamples=1)
    mod1.sim_snps(nloci=1, nsites=100000)
    genos1 = mod1.write_vcf()
    data=genos1.iloc[:, 9:].T
    print(data)
    test = SimulateNull(tree=testtree, model='ARD', matrix=data)
    test.significance_test()
    print(test.likelihoods)