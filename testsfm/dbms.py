import pandas as pd


# class Database(object):
#     def __init__(self,dbfileName):
#         pass

def select_datasets(database,datasets):
    return database[database["DATASET"].isin(datasets)]

def remove_datasets(database,datasets):
    return database[~database["DATASET"].isin(datasets)]