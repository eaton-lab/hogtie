#! usr/bin/env python

"""
Command line interface for hogtie
"""

import argparse
import sys
#import pandas as pd #eventually change this to convert matrix to df when I can parse matrix
import toytree
from hogtie import Pagel



def parse_command_line():
    """
    Parses the args for the Pagel class
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('--data',
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
    mydata = args.data.read()
    mydata = mydata.split(',')
    mydata = [int(i) for i in mydata] 
 
    mytree = args.tree.read()
    mytree = toytree.tree(mytree, tree_format=0)
 
    print('Getting conditional likelihoods...')
    liketree = Pagel(tree=mytree, data=mydata)
    liketree.run()
    print('The conditional likelihoods are:')
    print(liketree.tree.get_node_values('likelihood',True,True))


#if __name__ == "__main__":
