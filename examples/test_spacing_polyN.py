
import synbiomts
import cPickle as pickle

# We're going to use shelve to store model predictions
# as a dictionary-like persistance object of pandas dataframes
import shelve

# Import models
import sys
sys.path.append('../models')

# Import private Salis Lab code (latest versions of RBS Calculator)
sys.path.append('/home/alex/Private-Code')
from DNAc import *

from PyVRNA import PyVRNA
RNAEnergyModel = PyVRNA(dangles=0)

handle = open('../models/rRNA_16S_3p_ends.p','r')
rRNA_16S_3p_ends = pickle.load(handle)
handle.close()

# Required function in RBS Calculators
def get_rRNA(organism):
    # temporary until I update database:
    if organism=="Corynebacterium glutamicum B-2784": organism = 'Corynebacterium glutamicum R'
    elif organism=="Pseudomonas fluorescens A506":    return 'ACCTCCTTT'
    else: pass
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

import RBS_Calculator_v2_1
def RBSCalc_v2_1(sequence,organism,temp,startpos):
    start_range = [0,startpos+1]
    rRNA = get_rRNA(organism)
    model = RBS_Calculator_v2_1.RBS_Calculator(sequence,start_range,rRNA)
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
    'predicted_5pUTR' : sequence[:best_start_pos],
    'predicted_CDS' : sequence[best_start_pos:],
    }

    # Save dot-parentheses structure for initial state mRNA structure
    bpx = RBS.initial_structure['bp_x']
    bpy = RBS.initial_structure['bp_y']
    viennafld = RNAEnergyModel.bp2vienna(length=[len(RBS.sequence)],bpx=bpx, bpy=bpy)
    results['initial_structure'] = viennafld

    warning_dict = _calc_warning_flags(RBS)
    results.update(warning_dict)

    return results

if __name__ == "__main__":

    # add models to interface.Container
    models = synbiomts.interface.Container()
    models.add(RBSCalc_v2_1)
    models.setform(["RBSCalc_v2_1"], x="dG_total", y="PROT.MEAN", yScale='ln', a1=-0.45)

    # define database filters
    filters = { "DATASET": ['Egbert_PNAS_2012'] }

    # Provide the pickled database file name
    dbfilename = '../geneticsystems.db'

    # customtest = synbiomts.analyze.ModelTest(transl_rate_models,dbfilename,filters,nprocesses=1,verbose=True)
    customtest = synbiomts.analyze.ModelTest(models,dbfilename,filters,add_data=True,verbose=True)
    customtest.run()

    # Let's explore model error!
    # We noticed that poly-A, poly-U, poly-AU, and poly-AC had different slopes
    # Hypothesis: dG_stacking of the spacer sequence is a free energy
    # that the ribosome must overcome (initial state)
    # But we can't predict dG_stacking easily, so we can use this dataset
    # to estimate dG_stacking(NN)

    # Use scipy.optimize.minimize
    import scipy.optimize as optimize

    



    # Write model predictions and statistics to Excel
    # with open("labels/labels1.txt","r") as f:
    #     predictLabels = [x.strip('\n') for x in f.readlines()]
    # with open("labels/labels_stats.txt","r") as f:
    #     statsLabels = [x.strip('\n') for x in f.readlines()]
    # customtest.to_excel('Egbert_Output',predictLabels,statsLabels)