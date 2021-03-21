#! usr/bin/env python

"""
Command line interface for hogtie
"""

import argparse
import os
import toytree
import pandas as pd
from hogtie import Pagel



def parse_command_line():
    """
    Parses the args for the Pagel class
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('--data',
        nargs='?',
        type=argparse.FileType('r'),
        help='Input is binary character trait data in a csv file'
        )

    parser.add_argument('--tree',
        nargs='?',
        type=argparse.FileType('r'),
        help='tree in newick format with edge lengths and support vaues'
        )

    args = parser.parse_args()
    return args

def main():
    """
    Runs Pagel on parsed args
    """
    args = parse_command_line()
    print('Reading in data...')
    mydata = args.data.read()
    mydata = pd.DataFrame(mydata)
    print('Data read in as a pandas dataframe')

    #mydata = list(mydata.split(','))

    #mytree = args.tree.read()
    #mytree = toytree.tree(mytree, tree_format=0)

    #liketree = Pagel(tree=mytree, data=mydata)
    #liketree.run()


#if __name__ == "__main__":
