"""
Main classes for test-sfm module

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.

"""

import copy_reg
import types
import cPickle as pickle
from itertools import product
import multiprocessing as mp
import scipy
import pandas

import dbms

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
    
    def __init__(self,models,datasets,dbfileName,nprocesses=mp.cpu_count(),recalc=False):
        '''Inputs:
        models (interface.Models obj) = see interface.Models
        datasets (list)               = list of datasets to study
        dbfileName (string)           = filename of the geneticsystems database
        nprocesses (int)              = number of processes to use with
                                        multiprocessing if 1, ModelTest
                                        does not use multiprocessing
        recalc (bool)                 = boolean to tell the testsystem
                                        to recalcualte model predictions
                                        on existing datasets'''

        if models.__name__ != "Models":
            raise Exception("Not an interface.Models object: {}.".format(models))

        try:
            handle = open(dbfileName,'r')
            database = pickle.load(handle)
            handle.close()
            assert str(type(database)) == "<class 'pandas.core.frame.DataFrame'>"
        except:
            raise Exception("Database filename: {} is not valid.".format(dbfileName))

        listed = database["PAPER"].cat.categories
        unlisted = [x for x in datasets if x not in listed]
        if unlisted:
            error = "These datasets are unlisted: " + ", ".join(unlisted)
            raise ValueError(error)
        
        assert nprocesses > 0,          "nprocesses should be an int > 0"
        assert isinstance(recalc,bool), "recalc should be boolean"

        self.models = models
        self.datasets = datasets
        self.database = database
        self.partialdb = dbms.select_datasets(database,datasets)
        self.dbfileName = dbfileName
        self.nprocesses = nprocesses
        self.recalc = True

    def run(self):

        # generate iterable

        # for model in self.models._listed():
        # for entry in self.partialdb
        # print len(self.partialdb)
        # for i in xrange(len(self.partialdb)):

        # pool.imap
        # chunksize = 15

        entries = [{k: self.partialdb[k][i] for k in list(self.partialdb)} for i in xrange(len(self.partialdb))]
        pool = mp.Pool(processes=self.nprocesses)
        output = pool.imap(self._wrap,product(self.models._listed(),entries))

        # convert model output to pandas dataframe!

        pass

    def _wrap(self,bundle):
        (name,entry) = bundle
        # the model wrapper will interpret the inputs of the interface.Model
        # and pull those values from the database

        # print self.models[name].__code__.co_varnames
        # self.models[name]()
        return 1

    def add_datasets(self,datasets):
        self.datasets = list(set(self.datasets+datasets))
        self.partialdb = dbms.select_datasets(self.database,self.datasets)

    def remove_datasets(self,datasets):
        for ds in datasets:
            self.datasets.remove(ds)
        self.partialdb = dbms.select_datasets(self.database,self.datasets)

    def add_model(self,alias,model):
        pass

    def remove_model(self,alias,model):
        pass

    def calcstats(self):
        pass

    def compare2models(self):
        pass

    def propertyerror(self):
        pass

    def _export(self):
        pass

    def export2mat(self):
        '''Export data to mat file store in table format
        (we use MATLAB for plotting figures)'''
        
        # scipy.io.savemat()
        pass


if __name__ == "__main__":
    pass