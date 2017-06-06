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

    # identifiers is used to uniquely identify sequence entries
    identifiers = ["SEQUENCE","SUBGROUP"]

    def __init__(self,models,dbfilename,filters={},recalc=False,add_data=False,nprocesses=mp.cpu_count(),verbose=False):
        '''Inputs:
        models (interface.Models obj) = see interface.Models
        dbfilename (string)           = filename of the geneticsystems database
        filters (dictionary)          = dictionary to filter dataset using dbms
        nprocesses (int)              = number of processes to use with
                                        multiprocessing if 1, ModelTest
                                        does not use multiprocessing
        recalc (bool)                 = boolean to tell the testsystem
                                        to recalcualte model predictions
                                        on existing datasets
        verbose (bool)                = talk to me'''

        if models.__name__ != "Models":
            raise Exception("Not an interface.Models object: {}.".format(models))

        try:
            handle = open(dbfilename,'r')
            database = pickle.load(handle)
            handle.close()
            assert str(type(database)) == "<class 'pandas.core.frame.DataFrame'>"
        except:
            raise Exception("Database filename: {} is not valid.".format(dbfilename))

        for i in identifiers:
            assert i in database.keys(), "{} isn't a database label.".format(i)

        # Add read_Excel pandas code so users can have a database in a spreadsheet

        if "DATASET" in filters:
            listed = database["DATASET"].cat.categories
            unlisted = [x for x in filters["DATASET"] if x not in listed]
            if unlisted:
                error = "These datasets are unlisted: " + ", ".join(unlisted)
                raise ValueError(error)
        
        assert nprocesses > 0,          "nprocesses should be an int > 0"
        assert isinstance(recalc,bool), "recalc should be boolean"
        assert isinstance(filters,dict), "filters should be a dictionary"

        self.models     = models
        self.dbfilename = dbfilename
        self.nprocesses = nprocesses
        self.recalc     = recalc
        self.verbose    = verbose
        self.database   = database
        self.filters    = filters
        self.results    = {}

    def run(self,filename=None):

        # Use self.filters to keep only genetic systems of interest
        if self.filters: db = dbms.filter(self.database,self.filters,False)
        else: db = self.database

        # Remove sequences that have already been calculated if self.recalc is False,
        # convert pandas dataframe into entries, a list of records (dictionaries)
        # then bundle these entries with the models that are available
        num_entries = []
        if not self.recalc and not filename is None:
            d = shelve.open(filename)
            for model in self.models.available:
                kargs = {i: d[mode][i] for i in self.identifiers}
                db = dbms.remove(db,kargs,ordered=True)
                entries = db.T.to_dict().values()
                num_entries.append((model,len(entries)))
                bundles += [(model,entry) for entry in entries]
            d.close()
        else:
            entries = db.T.to_dict().values()
            bundles = product(self.models.available,entries)

        print len(entries)
        print len(bundles)
        wait = input(" ")

        if self.nprocesses > 1:
            pool = mp.Pool(processes=self.nprocesses) # consider chunking
            output = pool.map(self._wrap,bundles)
            pool.close()
            pool.join()
        else:
            output = [self._wrap(bundle) for bundle in bundles]

        total = 0
        for model,n in num_entries:
            modelcalcs = pandas.DataFrame(output[total:total+n])
            if add_data:
                dfsave = pd.concat([db, modelcalcs], axis=1)
            else:
                dfsave = pd.concat([db[self.identifiers], modelcalcs], axis=1)
            self.results[model] = dfsave
            print dfsave
            total += n

        wait = input(" ")

        # write to shelve
        if not filename is None:
            d = shelve.open(filename)
            if not d:
                d.update(self.results)
            else:
                for model in self.models.available:
                    d[model] = dbms.udpate(d[model],self.results[model],self.identifiers)
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

    def remove_datasets(self,datasets):
        self.datasets = [ds for ds in self.datasets[:] if ds not in datasets]

    def add_filters(self,filters):
        self.filters.update(filters)

    def add_model(self,alias,model,args,kargs):
        self.models.add(alias, model, args, kargs)

    def remove_model(self,alias):
        self.models.remove(alias)

    def calcstats(self):
        pass

    def compare2models(self):
        pass

    def propertyerror(self):
        pass

    def export(self):
        pass

    def export2mat(self):
        '''Export data to mat file store in table format
        (we use MATLAB for plotting figures)'''
        
        # scipy.io.savemat()
        pass


if __name__ == "__main__":
    pass