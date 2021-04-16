#! usr/bin/env python

"""
Command line interface for hogtie
"""

import argparse
import sys
import os
import pandas as pd
from hogtie import BinaryStateModel, MatrixParser, genomegraph



def parse_command_line():
    """
    Parses the args CLI inputs
    """

    parser = argparse.ArgumentParser('Welcome to hogtie!')

    parser.add_argument('-d', '--data',
        nargs='?',
        type=argparse.FileType('r'),
        help='Input is binary character trait data in a csv file',
        default=sys.stdin
        )

    parser.add_argument('-t', '--tree',
        nargs='?',
        type=argparse.FileType('r'),
        help='tree in newick format with edge lengths and support values'
        )

    parser.add_argument('-m', '--model',
        nargs='?',
        type=str,
        help='User must specify either ER (equal rates) or ARD (all rates different)'
        )

    parser.add_argument('-p', '--prior',
        nargs='?',
        type=float,
        help='Prior probability that the root state is 1 (default=0.5). Flat, uniform prior is assumed.'
        )

    args = parser.parse_args()
    return args

def main():
    """
    Runs Pagel on parsed args
    """
    args = parse_command_line()
   
    print('Reading in data and tree...')
    #mydata = args.matrix.read()
    mytree = args.tree.read()
   
    print('Calculating likelihoods...')
    liketree = MatrixParser(tree=mytree, matrix=args.data, model=args.model)
    liketree.matrix_likelihoods()


    HOGTIEDIR = os.path.dirname(os.getcwd())
    path = f"{HOGTIEDIR}/hogtie_output"
    os.makedirs(path, exist_ok=True)

    result = open(f"{HOGTIEDIR}/hogtie_output/result.csv", "w")

    liketree.likelihoods.to_csv(result)

    print(f'Wrote log-likelihoods to {HOGTIEDIR}/hogtie_output/result.csv.')

if __name__ == "__main__":
    HOGTIEDIR1 = os.path.dirname(os.getcwd())
    print(HOGTIEDIR1)
    path1 = f"{HOGTIEDIR1}/hogtie_output"
    print(path1)
