''' Let's identify the most predictive anti-SD sequence for future use with
RBS Calculator v2.0. Currently the model uses the 9 3'-most nucleotides of
the 16S rRNA, however there are 13 total unbase-paired nucleotides, more or
fewer of which may be interacting with the mRNA ribosome binding site.'''

import testsfm
import cPickle as pickle

# We're going to use shelve to store model predictions
# as a dictionary-like persistance object of pandas dataframes
import shelve

# Import RBS Calculator
import sys
sys.path.append('/home/alex/Private-Code')

import RBS_Calculator_v2_0
def RBSCalc_v2(sequence,rRNA,temp,startpos):
    start_range = [0,startpos+1]
    model = RBS_Calculator_v2_0.RBS_Calculator(sequence,start_range,rRNA)
    model.temp = temp
    model.run()
    output = model.output()
    (RBS,best_start_pos) = find_best_start(sequence,startpos,output)
    
    results = {
    'TIR' : RBS.tir,
    'dG_total' : RBS.dG_total,
    'dG_mRNA_rRNA' : RBS.dG_mRNA_rRNA,
    'dG_mRNA' : RBS.dG_mRNA,
    'dG_start' : RBS.dG_start,
    'dG_standby' : RBS.dG_standby,
    'dG_spacing' : RBS.dG_spacing,
    'used_mRNA_sequence' : RBS.sequence,
    'dG_UTR_fold' : RBS.dG_UTR_folding,
    'dG_SD_16S_hybrid' : RBS.dG_SD_hybridization,
    'dG_after_footprint' : RBS.dG_after_footprint_folding,
    'spacing_length' : RBS.spacing_length,
    'final_post_cutoff' : RBS.adaptive_post_cutoff,
    'binding_fragment' : RBS.optimum_hairpin_fragment,
    'surface_area' : RBS.optimum_available_RNA_surface_area,
    'dG_distortion' : RBS.optimum_distortion_energy,
    'dG_sliding' : RBS.optimum_sliding_energy,
    'dG_unfolding' : RBS.optimum_unfolding_energy,  
    'most_5p_paired_SD_mRNA' : RBS.most_5p_paired_SD_mRNA,
    'most_3p_paired_SD_mRNA' : RBS.most_3p_paired_SD_mRNA,
    'aligned_most_5p_paired_SD_mRNA' : RBS.aligned_most_5p_SD_border_mRNA,
    'aligned_most_3p_paired_SD_mRNA' : RBS.aligned_most_3p_SD_border_mRNA,
    '5pUTR' : sequence[:best_start_pos],
    'CDS' : sequence[best_start_pos:],
    }

    return results

# last 13 nt of Escherichia coli 16S rRNA sequence
rRNA_16S_end = "GAUCACCUCCUUA"

# Generate all possible anti-SD sequences of length 6-13 nt
candidates = [(rRNA_16S_end[i:i+j],i,j) for i in xrange(10) for j in xrange(4,14) if i+j<(len(rRNA_16S_end)+1)]
candidates = sorted(candidates, key=lambda x: (x[2], x[1]))

# Use subset of sequences with E. coli only



if __name__ == "__main__":

    # add models to interface.Models
    transl_rate_models = testsfm.interface.Models()
    transl_rate_models.add("RBSCalc_v2",RBSCalc_v2)

    # define database filters
    filters = {
                "DATASET": ['EspahBorujeni_NAR_2013',
                            'EspahBorujeni_NAR_2015',
                            'EspahBorujeni_JACS_2016',
                            'EspahBorujeni_Footprint',
                            'Salis_Nat_Biotech_2009',
                            'Farasat_MSB_2014',
                            'Tian_NAR_2015',
                            'Bonde_NatMethods_IC_2016'],

                "ORGANISM": ['Escherichia coli str. K-12 substr. MG1655',
                             'Escherichia coli str. K-12 substr. DH10B',
                             'Escherichia coli BL21(DE3)']
                }

    # Provide the pickled database file name
    dbfilename = '../geneticsystems.db'

    # customtest = testsfm.analyze.ModelTest(transl_rate_models,datasets,dbfilename,nprocesses=1,verbose=True)
    customtest = testsfm.analyze.ModelTest(transl_rate_models,datasets,dbfilename,verbose=True)
    customtest.run(filename='model_calcs.db')
