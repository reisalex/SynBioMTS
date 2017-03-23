"""
Main classes for test-sfm module

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.

"""

# import sys
# import os
# import re
# import time
import cPickle as pickle
import multiprocessing as mp
import scipy
import pandas
import copy_reg

# Using copy_reg to allow the pickle module to instance methods
# Replicates multiprocessing module's ForkingPickler
def _reduce_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)
copy_reg.pickle(types.MethodType, _reduce_method)


class ModelTest(object):
    '''
    ModelTest is the core of testsfm
    '''
    
    def __init__(self,models,datasets,database,nprocesses=mp.cpu_count(),recalc=False):
        '''Inputs:
        models (interface.Models obj)   = see interface.Models
        datasets (list)                 = list of datasets to study
        database (testsfm database obj) = see dbms.DatabaseManager
        nprocesses (int)                = number of processes to use with
                                          multiprocessing if 1, ModelTest
                                          does not use multiprocessing
        recalc (bool)                   = boolean to tell the testsystem
                                          to recalcualte model predictions
                                          on existing datasets'''

        if models.__name__ != "Models":
            raise Exception("Not an interface.Models object: {}.".format(models))

        available = db["PAPER"].cat.categories
        unavailable = [x for x in datasets if x not in available]
        if unavailable:
            error = "These datasets are unavailable: " + ", ".join(unavailable)
            raise ValueError(error)
        
        assert nprocesses > 0,          "nprocesses should be an int > 0"
        assert isinstance(recalc,bool), "recalc should be boolean"

        self.models = models
        self.datasets = datasets
        self.database = database
        self.nprocesses = nprocesses
        self.recalc = True

    def run(self):




        # pool.imap
        # chunksize = 15

        pass

    def add_datasets(self,datasets):
        self.datasets = self.datasets + datasets

    def _export(self):
        pass

    def calcstats(self):
        pass

    def compare2models(self):
        pass

    def propertyerror(self):
        pass

    def export2mat(self):
        '''Export data to mat file store in table format
        (we use MATLAB for plotting figures)'''
        
        # scipy.io.savemat()
        pass


if __name__ == "__main__":
    pass