
import testsfm
import cPickle as pickle

# We're going to use shelve to store model predictions
# as a dictionary-like persistance object of pandas dataframes
import shelve

# Import RBS Calculator
import sys
sys.path.append('/home/alex/test-sfm/models')
sys.path.append('/home/alex/Private-Code')
from DNAc import *

handle = open('/home/alex/test-sfm/models/rRNA_16S_3p_ends.p','r')
rRNA_16S_3p_ends = pickle.load(handle)
handle.close()

# Required function in RBS Calculators and UTR Designer
def get_rRNA(organism):
    return rRNA_16S_3p_ends[organism]

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


# Load and wrap models
import RBS_Calculator_v1_0
def RBSCalc_v1(sequence,startpos,organism,temp):
    start_range = [0,startpos+1]
    rRNA = get_rRNA(organism)
    model = RBS_Calculator_v1_0.RBS_Calculator(sequence,start_range,rRNA)
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
    'warnings' : RBS.warnings,
    '5pUTR' : sequence[:best_start_pos],
    'CDS' : sequence[best_start_pos:]
    }

    return results

import RBS_Calculator_v1_1
def RBSCalc_v1_1(sequence,startpos,organism,temp):
    start_range = [0,startpos+1]
    rRNA = get_rRNA(organism)
    model = RBS_Calculator_v1_1.RBS_Calculator(sequence,start_range,rRNA)
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
    'warnings' : RBS.warnings,
    '5pUTR' : sequence[:best_start_pos],
    'CDS' : sequence[best_start_pos:]
    }

    return results

import RBS_Calculator_v2_0
def RBSCalc_v2(sequence,organism,temp,startpos):
    start_range = [0,startpos+1]
    rRNA = get_rRNA(organism)
    model = RBS_Calculator_v2_0.RBS_Calculator(sequence,start_range,rRNA)
    model.temp = temp
    model.run()
    output = model.output()
    (RBS,best_start_pos) = find_best_start(sequence,startpos,output)
    
    results_dict = {
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
    'dot_bracket_structure' : RBS.bracket_string_mRNA,
    '5pNterminus_ssRNA' : RBS.ssRNA,
    'len_ssRNA' : RBS.ssRNA_len,
    'AU_count' : RBS.AU_count,
    'spacer_sequence' : RBS.spacer_seq,
    'mRNA_rRNA_hybrid_seq': RBS.mRNA_rRNA_hybrid_seq,
    'most_5p_paired_SD_mRNA' : RBS.most_5p_paired_SD_mRNA,
    'most_3p_paired_SD_mRNA' : RBS.most_3p_paired_SD_mRNA,
    'aligned_most_5p_paired_SD_mRNA' : RBS.aligned_most_5p_SD_border_mRNA,
    'aligned_most_3p_paired_SD_mRNA' : RBS.aligned_most_3p_SD_border_mRNA,
    'dG_mRNA_post_footprint' : RBS.dG_mRNA_post_footprint,
    'dG_footprint_subtracted' : RBS.dG_footprint_subtracted,
    '5pUTR' : sequence[:best_start_pos],
    'CDS' : sequence[best_start_pos:],
    }

    return results

import UTR_Designer
def wrap_UTR_Designer(sequence,startpos,organism,temp):
    start_range = [0,startpos+1]
    rRNA = get_rRNA(organism)
    model = UTR_Designer.RBS_Calculator(sequence,start_range,rRNA)
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
    'warnings' : RBS.warnings,
    '5pUTR' : sequence[:best_start_pos],
    'CDS' : sequence[best_start_pos:],
    }

    return results

handle = open('../models/RBSDesigner_ALL.p','r')
RBS_Designer = pickle.load(handle)
handle.close()
def wrap_RBS_Designer(sequence):
    # RBS['translation efficiency']
    # RBS['Expression']
    # RBS['exposure probability']
    # RBS['mRNA-ribosome probability']
    # RBS['SD']
    # RBS['dG_SD:aSD']
    # RBS['spacer length']
    return RBS_Designer[sequence]

from emopec._emopec import predict_spacing, get_expression
def EMOPEC(sequence,startpos):
    RBS = predict_spacing(sequence[:startpos])
    preseq,SD,spacer,Expression = RBS
    maxExpression = get_expression(sd_seq='AGGAGA',sd_dist=len(spacer))
    minExpression = get_expression(sd_seq='TTGGGC',sd_dist=len(spacer))
    Expression_percent = (Expression - minExpression)/(maxExpression - minExpression)

    results = {
    'preseq' : preseq,
    'SD_seq' : SD,
    'spacer_seq' : spacer,
    'spacer_length' : len(spacer),
    'Expression' : Expression,
    'Expression_percent': Expression_percent
    }

    return results

if __name__ == "__main__":

    # register models with the interface.Models object
    transl_rate_models = testsfm.interface.Models()
    # transl_rate_models.register("RBSCalc_v1",RBSCalc_v1)
    # transl_rate_models.register("RBSCalc_v1_1",RBSCalc_v1_1)
    # transl_rate_models.register("RBSCalc_v2",RBSCalc_v2)
    # transl_rate_models.register("UTRDesigner",wrap_UTR_Designer)
    # transl_rate_models.register("RBSDesigner",wrap_RBS_Designer)
    transl_rate_models.register("EMOPEC",EMOPEC)

    # Define datasets to run model calculations on
    # In this case, we're goign to specify the 856IC datasets
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

    # Provide the pickled database file name
    dbfileName = '../geneticsystems.db'

    customtest = testsfm.analyze.ModelTest(transl_rate_models,datasets,dbfileName,nprocesses=1,verbose=True)
    # customtest.run()