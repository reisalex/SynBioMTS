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

# Pipeline
# http://scikit-learn.org/stable/modules/generated/sklearn.pipeline.Pipeline.html#sklearn.pipeline.Pipeline
from sklearn.pipeline import Pipeline

# Principle Component Analysis
# http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
from sklearn.decomposition import PCA
