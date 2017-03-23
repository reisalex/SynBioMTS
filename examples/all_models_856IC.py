
import testsfm
import cPickle as pickle

transl_rate_models = testsfm.interface.Models()

def RBS_Calculator_v2_0_wrapper(SEQUENCE,ORGANISM,TEMP,STARTPOS):
    print SEQUENCE,ORGANISM,TEMP,STARTPOS

transl_rate_models.register("RBSCalcv1",RBS_Calculator_v2_0_wrapper,"ORGANISM","STARTPOS",SEQUENCE="ACTGTAC",TEMP=37.0)
transl_rate_models.RBSCalcv1("E. coli",10)

transl_rate_models.register("RBSCalcv2",RBS_Calculator_v2_0_wrapper)
transl_rate_models.RBSCalcv2("ACTGTAC","E. coli",37.0,10)

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

# load database
fileName = '../geneticsystems.db'

customtest = testsfm.analyze.ModelTest(transl_rate_models,datasets,fileName)
customtest.run()