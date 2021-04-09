#! usr/bin/env python

"""
Command line interface for hogtie
"""

import argparse
import sys
import os
import subprocess
import pandas as pd
#import pandas as pd #eventually change this to convert matrix to df when I can parse matrix
import toytree
from hogtie import BinaryStateModel, MatrixParser, genomegraph



def parse_command_line():
    """
    Parses the args for the Pagel class
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('--matrix',
        nargs='?',
        type=argparse.FileType('r'),
        help='Input is binary character trait data in a csv file',
        default=sys.stdin
        )

    parser.add_argument('--tree',
        nargs='?',
        type=argparse.FileType('r'),
        help='tree in newick format with edge lengths and support values'
        )

    args = parser.parse_args()
    return args

def main():
    """
    Runs Pagel on parsed args
    """
    print('Reading in data and tree...')
    args = parse_command_line()
   
    #mydata = args.matrix.read()
    mytree = args.tree.read()
   
    print('Calculating likelihoods...')
    liketree = MatrixParser(tree=mytree, matrix=args.matrix)
    liketree.matrix_likelihoods()

    HOGTIEDIR = os.path.dirname(os.getcwd())
    subprocess.run(["mkdir", f"{HOGTIEDIR}/hogtie_output"], check=True)
    subprocess.run(["touch", f"{HOGTIEDIR}/hogtie_output/results.csv"], check=True)

    liketree.likelihoods.to_csv(path_or_buf=f"{HOGTIEDIR}/hogtie_output/results.csv", mode='b')

#if __name__ == "__main__":
