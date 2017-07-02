import pandas as pd
import numpy as np

# Methods to update model calcs dataframe when new calculations are run
def update(db1,db2,identifiers):

    assert db1.keys() in db2.keys(), "Dataframes are not compatible on using dmbs.update(). \
                                      db1 keys: {} and db2 keys: {}".format(db1.keys(),db2.keys())
    for label in identifiers:
        assert label in db2.keys(), "{} is not a label in db2 in dbms.update.".format(label)

    indexes = get_indexes(db1,kargs,ordered=True,allqueries=True)
    db1length = len(db1)
    counter = 0
    for i in range(len(indexes)):
        if indexes[i] is None:
            indexes[i] = db1length+counter
            counter += 1
    
    assert len(db2)==len(indexes), "Length of indexes to add does not match length of db2 in dbms.udpate()"
    db2['indexes'] = indexes
    db2 = db2.set_index('indexes')
    db1.update(db2)

    return db1

def filter(database,kargs,ordered):
    '''Use to filter a pandas dataframe to get subset with kargs
    Input:  database (pandas dataframe)
            kargs (dictionary)  :: keys=database label, values=a list of filter values
            ordered (bool)      :: specifies if kargs contains lists with 
                                   of record attributes that are ordered
    Output: Database filtered to only include kargs'''
    return database[get_indexes(database,kargs,ordered)]

def remove(database,kargs,ordered):
    '''Use to filter a pandas dataframe to exclude entries that match kargs
    Input:  database (pandas dataframe)
            kargs (dictionary)  :: keys=database label, values=a list of filter values
            ordered (bool)      :: specifies if kargs contains lists with 
                                   of record attributes that are ordered            
    Output: Database filtered to remove kargs'''
    return database[~get_indexes(database,kargs,ordered)]

def get_indexes(database,kargs,ordered,allqueries=False):
    '''Gets indexes for values in database that match kargs'''
    assert isinstance(ordered,bool), "Argument, ordered, should be a boolean. Type given = {}".format(type(ordered))
    kargs = {k.upper(): v for k,v in kargs.iteritems()}
    keylist = []
    for key in kargs.keys():
        assert key in database.keys(), "{} is not a label in the database".format(key)
        keylist.append(key)
    if ordered:
        # Can probably be written using pandas functionalities like query()
        indexes = []
        for valuetup in zip(kargs.itervalues()):
            boolarray = np.all([database[keylist[i]]==valuetup[i] for i in range(len(valuetup))],axis=1)
            
            for i,x in enumerate(boolarray):
                if x:
                    print i
            indx = [i for i,x in enumerate(boolarray) if x]
            # print indx
            if not indx: indexes.append(None)
            else:        indexes += indx
    else:
        indexes = database[kargs.keys()].isin(kargs).all(1)
    if not allqueries:
        indexes = [indx for indx in indexes[:] if not indx is None]
    # print sum(indexes)
    return indexes

def to_dict_list(database):
    '''Convert pandas dataframe to list of dictionaries, where each entry (row)
    is a dictionary with the dataframe labels as keys (sorted by index).'''
    db_as_dict = database.T.to_dict()
    indxs = sorted(db_as_dict.keys())
    return [db_as_dict[i] for i in indxs]

#=== DEPRECATED ===#

def select_datasets(database,datasets):
    # For more general usage, use filter() with label="DATASET"
    return database[database["DATASET"].isin(datasets)]

def remove_datasets(database,datasets):
    # For more general usage, use remove() with label="DATASET"
    return database[~database["DATASET"].isin(datasets)]