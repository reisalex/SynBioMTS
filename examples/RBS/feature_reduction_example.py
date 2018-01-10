
import pandas as pd
import numpy as np
import csv
from synbiomts import learn,dbms

# Consider patial database
filters = { "DATASET": ['EspahBorujeni_NAR_2013',
                        'EspahBorujeni_NAR_2015',
                        'EspahBorujeni_JACS_2016',
                        'EspahBorujeni_Footprint',
                        'EspahBorujeni_Bsubtilis_2016',
                        'Salis_Nat_Biotech_2009',
                        'Farasat_MSB_2014',
                        'Tian_NAR_2015',
                        'Mimee_Cell_Sys_2015',
                        'Bonde_NatMethods_IC_2016',
                        'Egbert_PNAS_2012']
}
                        # 'Hecht_NAR_2017',
                        # 'Beck_PLoS_2016'

# RBS Calculator v2.0
db = pd.read_excel(open('spreadsheets/all_models_1014IC.xlsx','rb'), sheetname="RBSCalc_v2")
db = dbms.filter(db,filters)

X = db.filter(['dG_mRNA','dG_mRNA_rRNA','dG_start','dG_standby','dG_spacing'])
X['dG_mRNA'] = -1.0*X['dG_mRNA'] # Make sure dG_mRNA is negative (initial state)
y = np.log(db['PROT.MEAN'].values)

reduced_models = learn.feature_reduction(X,y)

combs = sorted(reduced_models.keys(),reverse=True)
table = [[c,reduced_models[c][0],reduced_models[c][1]] for c in combs]
with open('quick_feature_reduct.csv','w') as csvfile:
    writer = csv.writer(csvfile)
    [writer.writerow(r) for r in table]

