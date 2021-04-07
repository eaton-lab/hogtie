#! usr/bin/env python

import toyplot
import numpy as np
import pandas as pd
from scipy import stats
from hogtie import BinaryStateModel, MatrixParser


#there should be a horizontal line at the expected value
#outliers different colors?

class GenomeVisualization(MatrixParser):
    """"
    Takes data input of likelihood scores in a pandas dataframe and computes
    rolling window average of likelihoods then and plots them along
    linear genome. If the likelihood deviates more than expected in the positive direction,
    is flagged through coloring it red.
    """
    def __init__(self):
        pass

    def get_graph(self):
        """
        Graphs rolling average of likelihoods along the linear genome, identifies regions that deviate
        significantly from null expectations
        """
        
        self.likelihoods['rollingav']=self.likelihoods.rolling(50,win_type='triang').mean()
        #data['z_score']=stats.zscore(data['rollingav'],nan_policy='omit')

        plot = toyplot.plot(
            self.likelihoods['rollingav'],
            width = 500,
            height=500,
            color = 'blue'
        )

        return plot


if __name__ == "__main__":

    # Create test dataset
    data1 = np.random.uniform(low=0.0, high=1.0, size=1000)
    df = pd.DataFrame(data1)
    df['rollingmean']=df.rolling(50,win_type='triang').mean()

    # Compute rolling mean as 'z_score'
    df['z_score']=stats.zscore(df['rollingmean'],nan_policy='omit')
    scores = list(df['z_score'])
