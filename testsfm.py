"""
Main classes for test-sfm module

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Please cite:

  Alexander C. Reis, and Howard M. Salis
  An automated model test system for systematic development and improvement of
  gene expression models, Nature Methods (2017)

"""

# import sys
# import os
# import re
# import time
import cPickle as pickle
import multiprocessing as mp
import scipy
import pandas

class testsystem(object):

    database = 'geneticsystems.db'
    nprocesses = mp.cpu_count()
    recalc = False
    pickle_output = True

    def __init__(self,models,datasets,database=None,nprocesses=-1,recalc=False):
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

    '''
    class class_B(object):
        def __init__(self,name="Alex"):
            self.name = name
        def say_my_name(self):
            print self.name

    class class_A(class_B):
        def __init__(self,name,code=100):
            class_B.__init__(self,name)
            self.code = code
        def say_goodbye(self):
            print "goodbye"

    A = class_A(name="Sean")
    A.say_my_name()
    A.say_goodbye()
    print A.code
    '''

    import models

    transl_rate_models = models.interface()

    datasets = ['EspahBorujeni_NAR_2013',
                'EspahBorujeni_NAR_2015',
                'EspahBorujeni_JACS_2016',
                'EspahBorujeni_Footprint',
                'Salis_Nat_Biotech_2009',
                'Farasat_MSB_2014',
                'Tian_NAR_2015',
                'Mimee_Cell_Sys_2015',
                'Bonde_NatMethods_IC_2016'
                ]

    testsystem(transl_rate_models,datasets)