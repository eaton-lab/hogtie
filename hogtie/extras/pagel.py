#! usr/bin/env python

"""
Implementation of Pagel's method for ancestral state reconstruction
"""

#need a total likelihood score, assume value of 0 at root?

import numpy as np
import toytree
from scipy.optimize import minimize
from scipy.linalg import expm
from loguru import logger

class Pagel:
    """
    Pagel class object, implements Pagel's (1994) method for ancestral character
    state reconstruction.

    Parameters:
    -----------
    tree: ToyTree
        A toytree object of the phylogeny of interest

    data: list
        List of data representing binary character states (0's and 1's) along the
        tips of the tree. The length of the list must be equal to the number of 
        tips on the tree

    Returns:
    -----------
    Right now: A ToyTree object with estimated likelihood scores at each node
    NOTE: I want it to return likelihood of the whole tree given the data
    """


    def __init__(self, tree, data):
        self.tree = tree
        self.data = data

    def assign_tip_like_values(self):
        """
        Assigns likelihood values to tree tips based on data. 
        """
        values = [{0:-(i-1),1:i} for i in self.data]
        keys = list(range(0, len(self.data), 1))
        valuesdict = dict(zip(keys,values))
        self.tree = self.tree.set_node_values(feature = "likelihood", values = valuesdict)
        return self.tree


    def cond_like(self, likeleft0, likeleft1, likeright0, likeright1, tL, tR, alpha, beta):
        """
        Calculates conditional likelihood of character states at each node
        """
        Q = np.array([[-alpha, alpha], [beta, -beta]])
        probleft = expm(Q*tL)
        probright = expm(Q*tR)

        #ancestor is 0
        left0 = probleft[0, 0] * likeleft0 + probleft[0, 1] * likeleft1
        right0 = probright[0, 0] * likeright0 + probright[0, 1] * likeright1
        like_zero = left0*right0

        #ancestor is 1
        left1 = probleft[1, 0] * likeleft0 + probleft[1, 1] * likeleft1
        right1 = probright[1, 0] * likeright0 + probright[1, 1] * likeright1
        like_one = left1*right1

        return {0: like_zero, 1: like_one}

    def pruning_alg_before(self, alpha=0.5, beta=0.5):
        """
        Runs Felsenstein's pruning algorithm before max likelihood estimation
        of the instantaneous rates of transition. Allows for these rates to be 
        written to features of the trees. Likelihoods are subsequently over-
        written by pruning_alg_after
        """
        for node in self.tree.treenode.traverse("postorder"):
            if len(node.children) == 2:
                child1 = node.children[0]
                child2 = node.children[1]
                likedict = self.cond_like(likeright0 = child1.likelihood[0],
                                     likeright1 = child1.likelihood[1],
                                     likeleft0 = child2.likelihood[0],
                                     likeleft1 = child2.likelihood[1],
                                     tR = child1.dist,
                                     tL = child2.dist,
                                     alpha = alpha,
                                     beta = beta)
                #print(likedict)
                node.likelihood = likedict


    def node_like(self, x0, likeleft0, likeleft1, likeright0, likeright1, tL, tR, anca):
        """
        Calculates total likelihood (likelihood of the node having the character state
        0 or the character state 1) at each node.
        """

        condlik = self.cond_like(likeleft0, likeleft1, likeright0, likeright1, tL, tR, x0[0], x0[1])

        # get full likelihood
        lik = (1 - anca) * condlik[0] + (anca) * condlik[1]
        if anca in [0., 1.]:
            lik /= 2

        return -lik #np.log(lik)

    def model_fit(self, likeleft0, likeleft1, likeright0, likeright1, tL, tR, anca):
        """
        Find the maximum likelihood estimate of the two
        rate model parameters given the data.
        """
        args = (likeleft0, likeleft1, likeright0, likeright1, tL, tR, anca)

        # ML estimate
        estimate = minimize(
            fun=self.node_like, 
            x0=np.array([1., 1.]),
            args=args,
            method='L-BFGS-B',
            bounds=((0, 10), (0, 10))
            )

        score = -1 * self.node_like(estimate.x, *args)
        result = {
            "alpha": round(estimate.x[0], 3),
            "beta": round(estimate.x[1], 3), 
            "lik": round(score, 3),
            "convergence": estimate.success,
            }
        return result

    def fit_model_at_nodes(self, anca = 0.5):
        """
        Fits the Markov state model at each individual node, traversing the tree in postorder.
        Stores alpha and beta values (instantaneous transition rates) at each node.
        """
        self.tree = self.tree.set_node_values('alpha')
        self.tree = self.tree.set_node_values('beta')
        for node in self.tree.treenode.traverse("postorder"):
            if len(node.children) == 2:
                child1 = node.children[0]
                child2 = node.children[1]
                model = self.model_fit(likeright0 = child1.likelihood[0],
                                  likeright1 = child1.likelihood[1],
                                  likeleft0 = child2.likelihood[0],
                                  likeleft1 = child2.likelihood[1],
                                  tR = child1.dist,
                                  tL = child2.dist,
                                  anca = anca)
                node.alpha = model['alpha']
                node.beta = model['beta']
        return self.tree

    def pruning_alg_after(self):
        """
        Runs Felsenstein's pruning algorithm on an input tree, given instantaneous transition
        rates alpha and beta. Assigns likelihood scores for characters states at each node.
        Specifically for binary state model. 
        """
        self.tree = self.fit_model_at_nodes()
        for node in self.tree.treenode.traverse("postorder"):
            if len(node.children) == 2:
                child1 = node.children[0]
                child2 = node.children[1]
                likedict = self.cond_like(likeright0 = child1.likelihood[0],
                                     likeright1 = child1.likelihood[1],
                                     likeleft0 = child2.likelihood[0],
                                     likeleft1 = child2.likelihood[1],
                                     tR = child1.dist,
                                     tL = child2.dist,
                                     alpha = float(node.alpha),
                                     beta = float(node.beta))
                node.likelihood = likedict
        return self.tree
        logger.debug(self.tree)

    def run(self):
        """
        Assigns likelihood values to tips and runs Felstein's pruning algorithm for 
        a list of data along the provided tree
        """
        self.assign_tip_like_values()
        #I don't understand why it needs like values at nodes to run the pruning alg
        #after model_fit_at_nodes
        self.pruning_alg_before()
        self.tree = self.pruning_alg_after()


if __name__ == "__main__":
    tree1 = toytree.rtree.unittree(ntips=10)
    data1 = [0,1,1,0,1,1,0,0,0,1]
    testobject = Pagel(tree=tree1, data=data1)
    testobject.run()
    print(testobject.tree.get_node_values('likelihood',True,True)) #works!
