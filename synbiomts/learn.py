"""
Machine learning tools for model improvement
*In development*

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.

"""

import numpy as np
import sklearn
import scipy
import scipy.stats

# Feature selection
# http://scikit-learn.org/stable/modules/feature_selection.html#feature-selection
from sklearn.feature_selection import SelectKBest # select k best features
from sklearn.feature_selection import chi2 # use chi-squared to find n best
from sklearn.feature_selection import RFE # recursive feature eliminiation

# Remove features with low variance (i.e. low functional diversity)
# http://scikit-learn.org/stable/modules/feature_selection.html#removing-features-with-low-variance
from sklearn.feature_selection import VarianceThreshold

from sklearn.feature_selection import f_regression
def ANOVA(X,y):
    '''Univariate linear regression tests
    Quick linear model for sequentially testing the effect of many regressors
    Using scikit learn's Feature selection toolbox
    Returns:
        F (array) = F-values for regressors
        pvalues (array) = p-values for F-scores'''

    # (m,n) = np.shape(X)
    # M = np.transpose( [list(coeffs)*m] )
    # mX = np.multipy(M,X)

    (F,pvalues) = f_regression(X,y)
    return (F,pvalues)


# Pipeline
# http://scikit-learn.org/stable/modules/generated/sklearn.pipeline.Pipeline.html#sklearn.pipeline.Pipeline
from sklearn.pipeline import Pipeline

# Principle Component Analysis
# http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
from sklearn.decomposition import PCA
