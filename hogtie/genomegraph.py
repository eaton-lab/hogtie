#! usr/bin/env python

import numpy as np
import pandas as pd
from scipy import stats
import toyplot


#there should be a horizontal line at the expected value
#outliers different colors?

def genome_vis(data):
    """"
    Takes data input of likelihood scores in a pandas dataframe and computes
    rolling window average of likelihoods then and plots them along
    linear genome. If the likelihood deviates more than expected in the positive direction,
    is flagged through coloring it red.
    """
    
    data['rollingav']=data.rolling(50,win_type='triang').mean()
    data['z_score']=stats.zscore(data['rollingav'],nan_policy='omit')
    
    plot = toyplot.plot(
        data['rollingav'],
        width = 500,
        height=500
    );
    return plot

