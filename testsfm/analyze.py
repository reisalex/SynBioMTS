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

class ModelTest(object):
    
    nprocesses = mp.cpu_count()
    recalc = False
    pickle_output = True

    def __init__(self,models,datasets,database,nprocesses=-1,recalc=False):
        '''Inputs:
        models (class)    = object model.interface that defines IO for each
        datasets (list)   = list of datasets to run calculations on
        database (string) = string that describes the pickled database file
        nprocesses (int)  = number of processes to use with multiprocessing
                            if 1, testsystem does not use multiprocessing
        recalc (bool)     = boolean to tell the testsystem to recalcualte
                            model predictions on existing datasets'''

        # set models interface
        if models.__name__ == "models.interface": self.models = models

        # load database
        if isinstance(database,str): self.database = database
        handle = open(self.database,'r')
        db = pickle.load(handle)
        handle.close()

        # get data for specified datasets only
        bad_datasets = [x for x in datasets if x not in db["PAPER"].cat.categories]
        if bad_datasets:
            message = "These datasets are not available: " + ", ".join(bad_datasets)
            raise ValueError(message)
        self.db = db[db['PAPER'].isin(datasets)]

        # set other options
        if nprocesses > 1: self.nprocesses = nprocesses
        if recalc: self.recalc = True

    def run(self):
        pass

    def _import(self):
        pass

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