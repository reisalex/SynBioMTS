import pandas as pd
import numpy as np

# For more general usage, use filter_by_label with label="DATASET"
def select_datasets(database,datasets):
    return database[database["DATASET"].isin(datasets)]

def remove_datasets(database,datasets):
    return database[~database["DATASET"].isin(datasets)]

def filter_by_label(database,kargs):
    '''Use to filter a pandas dataframe to get subset with kargs
    Input:  database (pandas dataframe)
            kargs (dictionary)  :: keys=database label, values=a list of filter values
    Output: Filtered database by kargs
    NOTE: Replaces use of function select_datasets'''
    kargs = {k.upper(): v for k,v in kargs.iteritems()}
    for key in kargs.keys():
        assert key in database.keys(), "{} is not a label in the database".format(key)
    getindx = database[kargs.keys()].isin(kargs).all(1)
    return database[getindx]