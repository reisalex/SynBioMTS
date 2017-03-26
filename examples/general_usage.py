
import testsfm

transl_rate_models = testsfm.interface.Models()

def RBS_Calculator_v2_0_wrapper(sequence,organism,temp,startpos):
    print sequence,organism,temp,startpos

# simple usage
transl_rate_models.register("RBSCalcv2",RBS_Calculator_v2_0_wrapper)
# transl_rate_models.RBSCalcv2("ACTGTAC","E. coli",37.0,10)

# more comlex usage where you define args/kargs
transl_rate_models.register("RBSCalcv1",RBS_Calculator_v2_0_wrapper,"ACTGTAC",organism="E. coli",temp=37.0)
# transl_rate_models.RBSCalcv1(startpos=10)

datasets = ['EspahBorujeni_NAR_2013']

# load database
fileName = '../geneticsystems.db'

customtest = testsfm.analyze.ModelTest(transl_rate_models,datasets,fileName,nprocesses=1)
customtest.run()

# Wrapper working well!