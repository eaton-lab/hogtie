#! usr/bin/env python

"""
test data module
"""

import ipcoal
import toytree


# if this is meant to be executed then you should put it in main
# otherwise it would be run if this file was imported. 
if __name__ == "__main__":

    # sample a 10 tip imbalanced species tree
    tree = toytree.rtree.imbtree(ntips=10, treeheight=1e5)

    # setup a demographic model with high ILS
    mod = ipcoal.Model(tree=tree, Ne=1e5, nsamples=1)

    # simulate 1 long chromosome
    mod.sim_loci(nloci=1, nsites=100000)

    # write haploid genotypes to VCF file
    genos = mod.write_vcf(outdir="../examples", name="highILS")
    # genos.iloc[:, 9:].T
