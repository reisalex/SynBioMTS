from scipy import stats,polyfit
import numpy as np

import sklearn

class statistics(object):

	def __init__(self):
		pass

	def corr(self,x,y,name="Pearson"):
		# Linear or rank correlation
		# Returns:
		#   r,rho, or tau (float) = the test statistic
		#   pvalue (float) = the p-value for the test

		assert len(x)==len(y), "Arrays x & Y must be of equal length to calculate correlation statistics."

		if name == "Pearson":
			(r,pvalue) = stats.pearsonr(x,y)
			return r,pvalue

		elif name == "Spearman":
			(rho,pvalue) = stats.spearmanr(x,y)
			return rho,pvalue

		elif name == "Kendall":
			(tau,pvalue) = stats.kendalltau(x,y)
			return tau,pvalue

		else:
			error = "The {} correlation is not available. \
			Please use 'Pearson', 'Spearman', or 'Kendall.'".format(name)
			ValueError(error)

	def vartest2(self,x,y,alpha=0.05,test="F"):
		# Two-sample F-test/Barlett-test/Levene-test for equal variances
		# Returns:
		# 	h (bool)  = test decision, True if the test rejects the null hypothesis
		# 	S (float) = the test statistic
		# 	pvalue (float) = the p-value for significance of test decision

		if test == "F":
			df1 = len(x)-1
			df2 = len(y)-1
			statistic = np.var(x)/np.var(y) # stat = F
			pvalue = stats.f.sf(stat, df1, df2)

		elif test == "Barlett":
			(statistic,pvalue) = stats.bartlett(x,y)

		elif test == "Levene":
			(statistic,pvalue) = stats.levene(x,y)
		else:
			error = "The {}-test for equal variances is not available. \
			Please use 'F', 'Barlett', or 'Levene'.".format(test)
			ValueError(error)

		h = pvalue < alpha
		# add confidence interval calculation later
		return (h,statistic,pvalue)

	def ttest2(self,x,y,alpha=0.05):
		# Two sample t-test
		# Returns:
		# 	h (bool)  = test decision, True if the test rejects the null hypothesis
		# 	t (float) = the test statistic
		#	pvalue (float) = the p-value for significance of test decision

		(t,pvalue) = stats.ttest_ind(x,y,equal_var=False,nan_policy='omit')
		h = pvalue < alpha
		return (h,t,pvalue)

	def fit_linear_model(self,x,y,slope=None):
		# Linear least squares (LSQ) linear regression (polynomial fit of degree=1)
		# Returns:
		# 	m (float) = slope of linear regression line of form (y = m*x + b)
		# 	b (float) = intercept of linear regression line

		assert len(x)==len(y), "Arrays x & Y must be of equal length to fit linear regression model."

		if slope == None:
			(m,b) = polyfit(x,y,deg=1)

		else:
			LSQ = lambda b: np.sum( (y-(slope*x+b))**2.0 )
			res = minimize(LSQ,x0=1,bounds=None)
			(m,b) = (slope,res.x[0])

		# calculate significance of linear regression
		# self.ANOVA()
		
		return (m,b)

	def ANOVA(self,X,y,coeffs):
		# Univariate linear regression tests
		# Quick linear model for sequentially testing the effect of many regressors
		# Using scikit learn's Feature selection toolbox
		# Returns:
		# 	F (array) = F-values for regressors
		# 	pvalues (array) = p-values for F-scores

		# Compute regressors as product of predictors and the coeffs of a model
		(m,n) = np.shape(X)
		M = np.transpose( [list(coeffs)*m] )
		mX = np.multipy(M,X)

		(F,pvalues) = sklearn.feature_selection.f_regression(mX,y)
		return (F,pvalues)

	def mad(self,x):
		# Median absolute deviation
		# Returns:
		# 	MAD (float) = the median absolute deviation of the values in X
		# 	diff (list) = a list of deviations of the values in X from the median of X
		median = np.median(x)
		diff = np.absolute(x - median)
		MAD = np.median(diff)
		return (MAD,diff)

	'''Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and Handle Outliers",
	The ASQC Basic References in Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.'''	
	def find_outliers(self,x,threshold=2):
		# Outlier detection method using Iglewicz & Hoaglin's modified Z-score
		# Returns a list of bools, where True is an outlier
		(MAD,diff) = self.mad(x)
		M = 0.6745 * diff / MAD
		outliers = M > threshold
		return outliers

	def histogram(self):

		# pdf
		# number of elements in each bin

		pass

	def empirical_cdf(self):
		# both raw count and probability
		pass

	def roc(self):
		# calculate roc curve
		pass

	def auroc(self):
		# calculate area under the roc curve
		pass





class information_theory(obj):

	def __init__(self):
		pass

	def kldiv(self):
		pass

	def shannon_entropy(self):
		pass

	def seq_entropy(self):
		pass

	def total_seq_entropy(self):
		pass

	def model_capacity(self,H_seq):
		pass

	def 












