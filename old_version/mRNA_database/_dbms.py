


import sqlite3


createDb = sqlite3.connect('mRNA.db')

queryCurs = createDb.cursor()



def createTable():
    queryCurs.execute('''CREATE TABLE datasets
    (id TEXT PRIMARY KEY, 





def main():
    pass







if __name__ == '__main__':
    main()
