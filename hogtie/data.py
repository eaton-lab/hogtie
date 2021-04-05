#! usr/bin/env python

"""
test data module
"""

import ipcoal
import toytree


# simulating a tree with admixture edges between divergent tips close to the present 
if __name__ == "__main__":

    # sample a 10 tip imbalanced species tree
    tree = toytree.rtree.imbtree(ntips=10, treeheight=1e5)

    # setup a demographic model with high ILS
    mod = ipcoal.Model(tree=tree, Ne=1e5, nsamples=1, admixture_edges=(3,5,0.5,0.01))

    # simulate 1 long chromosome
    mod.sim_snps(1000)

    # write haploid genotypes to VCF file
    genos = mod.write_vcf(outdir="../examples", name="highILS")
