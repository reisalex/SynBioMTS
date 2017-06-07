
import sys
sys.path.append('../')
sys.path.append('../models')
sys.path.append('../datasets')
sys.path.append('/usr/local/lib/python2.7/site-packages/')

import testsfm
import cPickle as pickle
import TranslationRateModels as tl

# We're going to use shelve to store model predictions
# as a dictionary-like persistance object of pandas dataframes
import shelve


if __name__ == "__main__":

    # add models to interface.Models
    transl_rate_models = testsfm.interface.Models()
    transl_rate_models.add("RBSCalc_v1",tl.RBSCalc_v1)
    transl_rate_models.add("RBSCalc_v1_1",tl.RBSCalc_v1_1)
    transl_rate_models.add("RBSCalc_v2",tl.RBSCalc_v2)
    transl_rate_models.add("RBSCalc_v2_1",tl.RBSCalc_v2)
    transl_rate_models.add("UTRDesigner",tl.wrap_UTR_Designer)
    transl_rate_models.add("RBSDesigner",tl.wrap_RBS_Designer)
    transl_rate_models.add("EMOPEC",tl.EMOPEC)

    # define database filters
    filters = { "DATASET" :  ['Farasat_MSB_2014']
    }
    # filters = { "DATASET": ['EspahBorujeni_NAR_2013',
#                             'EspahBorujeni_NAR_2015',
#                             'EspahBorujeni_JACS_2016',
#                             'EspahBorujeni_Footprint',
#                             'Salis_Nat_Biotech_2009',
#                             'Farasat_MSB_2014',
#                             'Tian_NAR_2015',
#                             'Mimee_Cell_Sys_2015',
#                             'Bonde_NatMethods_IC_2016',
#                             'Egbert_Spacers_PNAS_2012']
#                 }

    # Provide the pickled database file name
    dbfilename = '../geneticsystems.db'

    # customtest = testsfm.analyze.ModelTest(transl_rate_models,dbfilename,filters,nprocesses=1,verbose=True)
    customtest = testsfm.analyze.ModelTest(transl_rate_models,dbfilename,filters,nprocesses=10, verbose=True)
    customtest.run(filename='model_calcs.db')