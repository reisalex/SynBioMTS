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
import shelve
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

    def __init__(self,models,datasets,dbfilename,nprocesses=mp.cpu_count(),recalc=False,verbose=False):
        '''Inputs:
        models (interface.Models obj) = see interface.Models
        datasets (list)               = list of datasets to study
        dbfilename (string)           = filename of the geneticsystems database
        nprocesses (int)              = number of processes to use with
                                        multiprocessing if 1, ModelTest
                                        does not use multiprocessing
        recalc (bool)                 = boolean to tell the testsystem
                                        to recalcualte model predictions
                                        on existing datasets'''

        if models.__name__ != "Models":
            raise Exception("Not an interface.Models object: {}.".format(models))

        try:
            handle = open(dbfilename,'r')
            database = pickle.load(handle)
            handle.close()
            assert str(type(database)) == "<class 'pandas.core.frame.DataFrame'>"
        except:
            raise Exception("Database filename: {} is not valid.".format(dbfilename))

        # Add read_Excel pandas code so users can have a database in a spreadsheet

        listed = database["DATASET"].cat.categories
        unlisted = [x for x in datasets if x not in listed]
        if unlisted:
            error = "These datasets are unlisted: " + ", ".join(unlisted)
            raise ValueError(error)
        
        assert nprocesses > 0,          "nprocesses should be an int > 0"
        assert isinstance(recalc,bool), "recalc should be boolean"

        self.models     = models
        self.datasets   = datasets
        self.dbfilename = dbfilename
        self.nprocesses = nprocesses
        self.recalc     = recalc
        self.verbose    = verbose
        self.database   = database
        self.partialdb  = dbms.select_datasets(database,datasets)
        self.results    = {}

    def run(self,filename=None):

        entries = ( {k: self.partialdb[k][i] for k in list(self.partialdb)} \
                for i in xrange(len(self.partialdb)) )
        bundles = product(self.models.available,entries)
        if self.nprocesses > 1:
            pool = mp.Pool(processes=self.nprocesses)
            output = pool.map(self._wrap,product(self.models.available,bundles))
            # Consider using (chunksize = n)
            pool.close()
            pool.join()
        else:
            output = [self._wrap(bundle) for bundle in bundles]

        # Convert model output to pandas dataframe!
        # And pickle
        chunksize = len(output)/len(self.models)
        output_by_model = [output[x:x+chunksize] for x in xrange(0,len(output),chunksize)]

        for i in xrange(len(self.models.available)):
            name = self.models.available[i]
            self.results[name] = pandas.DataFrame(output_by_model[i])
            self.results[name]['ID'] = self.partialdb['ID'] # ADD IDs

        # write to shelve
        if filename is not None:
            d = shelve.open(filename)
            d.update(self.results)
            d.close()

    def _wrap(self,bundle):
        # the model wrapper will interpret the inputs of the interface.Model
        # and pull those values from the database

        (name,entry) = bundle
        # entry = {'ORGANISM':"Escherichia coli",'SEQUENCE':"ACTCGATCTT",...}

        # dev notes
        #('ACTGTAC',) # args
        #{'organism': 'E. coli', 'temp': 37.0} # keywords
        #['sequence', 'organism', 'temp', 'startpos'] # variables

        # Remove args and keywords from variables list
        vrs = self.models[name].variables[len(self.models[name].args):]
        vrs = [k for k in vrs[:] if k not in self.models[name].keywords.keys()]

        # Exception handling when data is not available
        if any(k.upper() not in entry.keys() for k in vrs):
            err = "One of {}'s arguments is not in the database.".format(name)
            print "Model requested arguments: " + str(vrs)
            print "Database available values: " + str(entry.keys())
            raise KeyError(err)
        else:
            kargs = {k: entry[k.upper()] for k in vrs}

        # Print if verbose
        if self.verbose:
            fmt = "model={m:s}, ID={ID:s}, seq={seq:s}"
            print fmt.format(m=name,ID=entry['ID'],seq=entry['SEQUENCE'][:50])

        # Run model
        # Using keyword argument unpacking to pass only requested args
        output = self.models[name](**kargs)

        return output

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