import testsfm

transl_rate_models = testsfm.models.Interface()

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

testsfm.analyze.ModelTest(transl_rate_models,datasets,database='../geneticsystems.db')