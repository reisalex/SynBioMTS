''' Let's identify the most predictive anti-SD sequence for future use with
RBS Calculator v2.0. Currently the model uses the 9 3'-most nucleotides of
the 16S rRNA, however there are 13 total unbase-paired nucleotides, more or
fewer of which may be interacting with the mRNA ribosome binding site.'''

import testsfm
import cPickle as pickle
import numpy as np
import pandas as pd

# We're going to use shelve to store model predictions
# as a dictionary-like persistance object of pandas dataframes
import shelve

# Import RBS Calculator
import sys
sys.path.append('/home/alex/Private-Code')
from DNAc import *

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

def find_best_start(mRNA, start_pos, predictions):
    '''Finds the number of start codons, the most highly translated start codon and the start range associated with it.
        Disallows start codons that have stop codons following them.'''
    
    stop_codons = ["TAA", "TAG", "TGA", "UAA", "UAG", "UGA"]
    num_start_codons = 0
    dG_tot_list = []
    dG_best = None
    
    for index in range(len(predictions.RBS_list)):
    
        RBS = predictions.RBS_list[index]
        in_frame = (start_pos - RBS.start_position) % 3
        stop_codon_present = any([ (mRNA[pos:pos+3] in stop_codons) for pos in range(RBS.start_position+3, start_pos, 3)]) # stop codon present?

        if (in_frame == 0) and (not stop_codon_present): # if start is in frame and no downstream stop codon
            num_start_codons += 1
            dG_tot_list.append(RBS.dG_total)

        if (dG_best == None or RBS.dG_total < dG_best) and (in_frame == 0) and (not stop_codon_present):
            best_start_index = index
            dG_best = RBS.dG_total
            best_start_pos = RBS.start_position
    
    RBS = predictions.RBS_list[best_start_index]
    return (RBS,best_start_pos)

def test_model(rRNA,calcs,ref_error=None):

    # Useful lambda functions
    calc_dG_apparent = lambda K,beta,prot: np.log(prot/K)/beta
    calc_TIR = lambda K,beta,dG: K*np.exp(-beta*dG)

    all_dG_apparent = np.array([])
    all_error = np.array([])

    # get each subset, and fit proportionality constant
    for SG in calcs["SUBGROUP"].unique():
        dG_total = np.array(calcs["dG_total"][calcs["SUBGROUP"]==SG])
        protein  = np.array(database["PROT.MEAN"][database["SUBGROUP"]==SG])

        # print SG
        # print dG_total
        # print protein

        # determine outliers with initial fit
        (beta,K) = testsfm.stats.fit_linear_model(dG_total,np.log(protein),slope=-0.45)
        dG_apparent = calc_dG_apparent(K,beta,protein)
        ddG = dG_total - dG_apparent
        ddG_abs = np.absolute(ddG)
        outliers = testsfm.stats.find_outliers(ddG_abs)
        keepers = np.invert(outliers)

        # reapply fit with trimmed dataset & calculate error
        dG_total_trimmed = dG_total[keepers]
        protein_trimmed  = protein[keepers]
        (beta,K) = testsfm.stats.fit_linear_model(dG_total_trimmed,np.log(protein_trimmed),slope=-0.45)
        dG_apparent = calc_dG_apparent(K,beta,protein)
        TIR = calc_TIR(K,0.45,dG_total)
        TIR_error = protein/TIR

        # save dG_apparent & ddG
        all_dG_apparent = np.append(all_dG_apparent,dG_apparent)
        all_error = np.append(all_error,TIR_error)

        # print K
        # wait = input("value=")

    # calculate R^2
    (R,Pearsonp) = testsfm.stats.correlation(np.array(calcs["dG_total"]),all_dG_apparent)
    Rsqr = R**2

    # calculate F-test for equal variances against "AGGAGG"
    # ddG_ref are the ddG values computed for "AGGAGG"
    if ref_error == None: ref_error = all_error
    (h,F,Ftestp) = testsfm.stats.vartest2(all_error,ref_error,logNormal=True,alpha=0.05,test="F")

    # calculate Kullback-Leibler divergence
    (NKLdiv,KLdiv,KLdivmax) = testsfm.stats.normKLdiv(all_error,b=4)

    # return results
    results = {
    'rRNA': rRNA,
    'Rsqr': Rsqr,
    'Pearsonp': Pearsonp,
    'F': F,
    'Ftestp': Ftestp,
    'TIR_error': all_error,
    'error var': np.exp(np.var(np.log(all_error))),
    'NKLdiv': NKLdiv
    }

    return results


if __name__ == "__main__":

    # last 13 nt of Escherichia coli 16S rRNA sequence
    rRNA_16S_end = "GAUCACCUCCUUA"

    # Generate all possible anti-SD sequences of length 6-13 nt
    candidates = [(rRNA_16S_end[i:i+j],i,j) for i in xrange(10) for j in xrange(4,14) if i+j<(len(rRNA_16S_end)+1)]
    candidates = sorted(candidates, key=lambda x: (x[2], x[1]))

    # define database filters
    # we are only going to test on E. coli for now
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
    
    # add models to interface.Models
    # model names are the tested 16S rRNA sequence
    transl_rate_models = testsfm.interface.Models()
    for (rRNA,pos,length) in candidates:
        transl_rate_models.add(alias=rRNA,model=RBSCalc_v2,rRNA=rRNA)

    customtest = testsfm.analyze.ModelTest(transl_rate_models,dbfilename,filters,verbose=True)
    customtest.run(filename='test_16S_rRNA.db')

    # Now that the test run, let's calculate the proportionality constants,
    # and the statistics (including F-test for equal variance, and R^2 values)

    # load the database
    handle = open(dbfilename,'r')
    database = pickle.load(handle)
    handle.close()
    
    # load the model calculations
    modelcalcs = shelve.open('test_16S_rRNA.db')

    # calculate reference model first
    rRNA = "CCUCCU"
    calcs = modelcalcs[rRNA]
    results = test_model(rRNA,calcs)
    ref_error = results["TIR_error"]
    # calculate all models and export to excel
    all_results = [test_model(rRNA,modelcalcs[rRNA],ref_error) for (rRNA,pos,length) in candidates]

    output = pd.DataFrame(all_results)
    output.to_excel('test_Ecoli_aSD_seqs.xlsx',sheet_name='Results')

    modelcalcs.close()