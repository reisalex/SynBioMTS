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
import numpy as np
import pandas
import shelve
import dbms
import stats

# Using copy_reg to allow the pickle module to instance methods
# Replicates multiprocessing module's ForkingPickler
def _reduce_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)
copy_reg.pickle(types.MethodType, _reduce_method)


class ModelTest(object):

    # identifiers is used to uniquely identify sequence entries
    identifiers = ["SEQUENCE","SUBGROUP"]

    def __init__(self,models,dbfilename,filters={},recalc=False,add_data=True,nprocesses=mp.cpu_count(),verbose=False):
        '''Inputs:
        models (interface.Container obj) = see interface.Container
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
            raise Exception("Not an interface.Container object: {}.".format(models))
        
        assert nprocesses > 0,            "nprocesses should be an int > 0"
        assert isinstance(recalc,bool),   "recalc should be boolean"
        assert isinstance(add_data,bool), "add_data should be boolean"
        assert isinstance(filters,dict),  "filters should be a dictionary"

        self.models      = models
        self.dbfilename  = dbfilename
        self.recalc      = recalc
        self.add_data    = add_data
        self.nprocesses  = nprocesses
        self.verbose     = verbose
        self.filters     = filters
        self.predictions = {}
        self.statistics  = {}

        # import sequences from genetic systems database
        # based on specified dbfilename and filters
        self._update_database()        

    def run(self,calcsFilename=None,statsFilename=None):
        self.predict(calcsFilename)
        self.calc_stats(statsFilename)

    def predict(self,filename=None):

        db = self.database

        # Remove sequences that have already been calculated if self.recalc is False,
        # convert pandas dataframe into entries, a list of records (dictionaries)
        # then bundle these entries with the models that are available
        n_entries = []
        if (not self.recalc) and (not filename is None):
            bundles = []
            d = shelve.open(filename)
            for model in self.models.available:
                if model in d.keys():
                    kargs = {i: d[model][i] for i in self.identifiers}
                    db = dbms.remove(db,kargs,ordered=True)
                    dict_list = dbms.to_dict_list(db)
                    n_entries.append((model,len(dict_list)))
                    bundles += [(model,entry) for entry in dict_list]
                else:
                    dict_list = dbms.to_dict_list(db)
                    bundles += product([model],dict_list)
                    n_entries.append((model,len(dict_list)))
            d.close()
        else:
            dict_list = dbms.to_dict_list(db)
            bundles = product(self.models.available,dict_list)
            n = len(dict_list)
            n_entries = [(m,n) for m in self.models.available]

        if self.nprocesses > 1:
            pool = mp.Pool(processes=self.nprocesses) # consider chunking
            output = pool.map(self._wrap,bundles)
            pool.close()
            pool.join()
        else:
            output = [self._wrap(bundle) for bundle in bundles]
        
        db.reset_index(drop=True, inplace=True)
        total = 0
        for model,n in n_entries:
            modelcalcs = pandas.DataFrame(output[total:total+n])
            if self.add_data:
                dfsave = pandas.concat([db, modelcalcs], axis=1)
            else:
                dfsave = pandas.concat([db[self.identifiers], modelcalcs], axis=1)
            self.predictions[model] = dfsave
            total += n

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

        if self.verbose:
            fmt = "model={m:s}, SUBGROUP={sbgrp:s}, seq={seq:s}..."
            print fmt.format(m=name,sbgrp=entry['SUBGROUP'],seq=entry['SEQUENCE'][:50])

        # Run model
        return self.models[name](**kargs)

    def _update_database(self):

        try:
            handle = open(self.dbfilename,'r')
            database = pickle.load(handle)
            handle.close()
            assert isinstance(database,pandas.DataFrame)
        except:
            raise Exception("Database filename: {} is not valid.".format(self.dbfilename))

        for i in self.identifiers:
            assert i in database.keys(), "{} isn't a database label.".format(i)

        # dbms.get_indexes checks that filters are labels in database
        # This code block code be removed if dbms.get_indexes gets no values
        if "DATASET" in self.filters:
            listed = database["DATASET"].cat.categories
            unlisted = [x for x in self.filters["DATASET"] if x not in listed]
            if unlisted:
                error = "These datasets are unlisted: " + ", ".join(unlisted)
                raise ValueError(error)

        # Use filters specified by user on database
        if self.filters:
            database = dbms.filter(database,self.filters,False)
        self.database = database

    def calc_stats(self,filename=None):
        
        # assert that model predictions are available
        # run statistics calculations for only datasets that have been calculated?

        for m in self.models.available:
            df = self.predictions[m]
            allError = np.inf*np.ones(len(df))
            entries = []
            allPredicted = np.array([])
            yScale = self.models[m].yScale

            for subgroup in df["SUBGROUP"].unique():

                # Extract subgroup data & predictions
                indx = df["SUBGROUP"] == subgroup
                x = np.array(df[self.models[m].x][indx])
                y = np.array(df[self.models[m].y][indx])

                # Run statistics and information theory calcs
                data,yError = stats.linear_complete(x,y,yScale,self.models[m].a1)
                data["Sequence entropy"],_ = stats.sequence_entropy(df["SEQUENCE"][indx],\
                                                                          positions=df["STARTPOS"][indx])
                data["SUBGROUP"] = subgroup
                entries.append(data)

                # Add yError to model predictions
                allError[indx] = yError

                # Calculate yPredicted based on linear model fit
                yPredicted = data['slope']*x + data['intercept']
                if yScale == 'ln':
                    yPredicted = np.exp(yPredicted)
                elif yScale == 'log10':
                    yPredicted = np.power(10,yPredicted)
                elif yScale == 'linear':
                    pass
                else:
                    raise Exception('Bad yScale set for {}'.format(m))
                allPredicted = np.concatenate((allPredicted,yPredicted))

            # Calculate statistics for model {m} on all data
            y = np.array(df[self.models[m].y])
            data,_ = stats.linear_simple(allPredicted,y,yScale)
            entries.append(data)

            self.statistics[m] = pandas.DataFrame(entries)
            self.predictions[m]['yError'] = pandas.Series(allError, index=df.index)

        # write statistics to shelve if filename given
        if not filename is None:
            pass
    
    def to_shelve(self,filename):
        
        if not filename is None:
            d = shelve.open(filename)
            if not d:
                d.update(self.predictions)
            else:
                for model in self.models.available:
                    d[model] = dbms.udpate(d[model],self.predictions[model],self.identifiers)
            d.close()

    def to_excel(self,filename,predictColumns=[],statsColumns=[],models=[]):
        assert isinstance(filename,str), "Filename provided, {}, for export() needs to be a string.".format(filename)
        if not models:
            models = self.predictions.keys()
        else:
            for m in models:
                assert m in self.predictions.keys(), "Model {} was not tested.".format(m)
        if filename[-4:] == ".xlsx":
            fn = filename
        else:
            fn = filename + ".xlsx"
        writer = pandas.ExcelWriter(fn)
        if predictColumns:
            for model in models:
                self.predictions[model].to_excel(writer,sheet_name=model,columns=predictColumns)
        if statsColumns:        
            for model in models:
                self.statistics[model].to_excel(writer,sheet_name="{}-stats".format(model),columns=statsColumns)
        writer.save()

    def to_csv(self):
        pass

    def add_datasets(self,datasets):
        self.datasets = list(set(self.datasets+datasets))
        self._update_database()

    def remove_datasets(self,datasets):
        self.datasets = [ds for ds in self.datasets[:] if ds not in datasets]
        self._update_database()

    def add_filters(self,filters={}):
        self.filters.update(filters)

    def remove_filters(sef,filters=[]):
        for key in filters:
            self.filters.pop(key,None)
        self._update_database()

    def add_model(self,alias,model,args,kargs):
        self.models.add(alias, model, args, kargs)

    def remove_model(self,alias):
        self.models.remove(alias)


if __name__ == "__main__":
    pass