import pandas as pd
import numpy as np
#import math
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
#import sys
#import os

## objective functions
def NSE(simulated,observed):
    num = np.power(observed - simulated,2).sum()
    den = np.power(observed - observed.mean(),2).sum()
    return 1 - num/den

def KGE(simulated,observed):
    pearson = observed.astype('float64').corr(simulated.astype('float64'),method='pearson')
    t1 = (pearson-1.0)**2.0
    t2 = (simulated.std()/observed.std()-1)**2
    t3 = (simulated.mean()/observed.mean()-1)**2
    return 1 - (t1+t2+t3)**0.5

def PBIAS(simulated,observed):
    return 100*(simulated-observed).sum()/simulated.sum()

def pkerr(simulated, observed):
    sim_max = np.max(simulated)
    obs_max = np.max(observed)
    
    return 100*(sim_max - obs_max) / obs_max

# call function
def call_func(x,y,func,fn_mapper):
    try:
        return fn_mapper[func](x,y)
    except:
        print(f"Invalid objective function requested: {func}")
        return np.inf


