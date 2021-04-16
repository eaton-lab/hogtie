#! usr/bin/env python

"""
generating data using ipcoal
"""

import ipcoal
import toytree
import pandas as pd
from scipy.stats import zscore
from hogtie import MatrixParser

class Hogtie:
    """
    Compare data to expectations
    """
    def __init__(self, tree, matrix, model=None, prior=0.5):
        self.tree = tree
        self.model = model
        self.prior = prior
        self.matrix = matrix

    def create_null(self):
        """
        Simulates SNPs across the input tree to create the null expectation for likelihood scores. 
        Deviation from the null will be flagged. 
        """

        #high ILS
        mod = ipcoal.Model(tree=self.tree)
        mod.sim_loci(nloci=1, nsites=100000)
        genos = mod.write_vcf()
        
        null_genos = genos.iloc[:, 9:].T

        #run the Binary State model on the matrix and get likelihoods
        null = MatrixParser(tree=self.tree, matrix=null_genos, model=self.model)
        null.matrix_likelihoods()

        #get z-scores null likelihood expectations
        sigma = null.likelihoods.std()
        mu = null.likelihoods.mean()
        self.sigma = sigma
        self.mu = mu


    def cutoff(self):
        """
        identifies k-mer likelihoods that are +/-3 outside of zero from simulated likelihood z-score
        """
        #get the likelihood value that corresponds to the z-score of 3 in simulations
        high_lik = self.mu + 2*self.sigma

        lik_calc = MatrixParser(tree=self.tree,
                               model=self.model,
                               prior=self.prior,
                               matrix=self.matrix
                               )

        lik_calc.matrix_likelihoods()
        
        deviations = []
        for like in lik_calc.likelihoods:
            if (like < high_lik).bool:
                deviations.append(0)
            else:
                deviations.append(1)

        self.deviations = deviations

    def genome_graph(self):
        """
        Graphs rolling average of likelihoods along the linear genome, identifies regions that deviate
        significantly from null expectations

        TO DO: change color of outliers 
        """

        self.likelihoods['rollingav']= self.likelihoods.rolling(10, win_type='triang').mean()
        
        plot = toyplot.plot(
            self.likelihoods['rollingav'],
            width = 500,
            height=500,
            color = 'blue'
        )

        return plot

if __name__ == "__main__":
    testtree = toytree.rtree.unittree(ntips=10)
    mod1 = ipcoal.Model(tree=testtree, Ne=1e3, admixture_edges=[(3, 8, 0.5, 0.5)], nsamples=1)
    mod1.sim_loci(nloci=1, nsites=100000)
    genos = mod1.write_vcf()
    data=genos.iloc[:, 9:].T
    test = Hogtie(tree=testtree, model='ARD', matrix=data)
    test.create_null()
    test.cutoff()
    print(test.deviations)