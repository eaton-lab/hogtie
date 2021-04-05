#! usr/bin/env python

"""
generating data using ipcoal
"""

import ipcoal
import toytree

if __name__ == "__main__":

    # sample a 10 tip imbalanced species tree
    tree = toytree.rtree.imbtree(ntips=10, treeheight=1e5)

    # set up a Model with an admixture edge between divergent tips close to the present
    mod = ipcoal.Model(tree=tree, Ne=1e5, nsamples=1, admixture_edges=(3,5,0.5,0.01))

    # simulate 1 long chromosome
    mod.sim_snps(1000)

    # write haploid genotypes to VCF file
    genos = mod.write_vcf(outdir="../examples", name="recent_admixture")
