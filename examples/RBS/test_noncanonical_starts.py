
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
    elif organism=="Escherichia coli BL21(DE3)":      return 'ACCTCCTTA'
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
def RBSCalc_v2_1(sequence,organism,temp,startpos,start_energies):

    start_range = [0,startpos+1]
    rRNA = get_rRNA(organism)
    model = RBS_Calculator_v2_1.RBS_Calculator(sequence,start_range,rRNA)
    model.temp = temp

    model.start_codon_energies = start_energies
    model.start_codons = model.start_codon_energies.keys()

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

    return results

if __name__ == "__main__":

    # We're going to use the NanoLuc luciferase measurements on both plasmids,
    # with p15A and the oriS/BAC origin of replications to parameterize dG_start
    filters = { 'DATASET': ['Hecht_NAR_2017'],
                'PROTEIN': ['NanoLuc']
                }

    dbfilename = '../geneticsystems.db'

    # We calculated the RNA duplex Gibbs free energies at 37 C using Turner 2004
    # nearest neighbor free energy paramters (see Supp #)
    dG_start_INN = {
    "ATG": -2.76,
    "GTG": -3.47,
    "TTG": -2.11,
    "CTG": -2.11,
    "ATA": -0.67,
    "ATC": -0.67,
    "ATT": -0.67,
    "CAT": 0.0,
    "TGC": 0.0,
    "CGC": 0.0,
    "GGA": 0.0,
    "TAG": 0.0
    }

    # add models to interface.Container
    models = synbiomts.interface.Container()
    models.add(RBSCalc_v2_1,start_energies=dG_start_INN)
    modelName = 'RBSCalc_v2_1'
    models.setform([modelName], x='dG_total', y='PROT.MEAN', yScale='ln', a1=-0.45)

    # customtest = synbiomts.analyze.ModelTest(models,dbfilename,filters,nprocesses=1,verbose=True)
    test = synbiomts.analyze.ModelTest(models,dbfilename,filters,add_data=True,verbose=True)
    test.run()

    # Let's fit the NanoLuc data to determine the apparent dG_starts
    '''
    import numpy as np
    from synbiomts import dbms

    calcs = test.predictions[modelName]
    stats = test.statistics[modelName]
    K = np.exp(stats['intercept'])
    beta = 0.45

    x = []
    for sg,K in zip(stats['SUBGROUP'],np.exp(stats['intercept'])):
        kargs = {'SUBGROUP':[sg]}
        df = dbms.filter(calcs,kargs,False)
        dG_apparent = np.log(df['PROT.MEAN']/K)/-beta
        ddG = dG_apparent - df['dG_total']
        x.append(df['dG_start'] + ddG)

    start_codons = [CDS[0:3] for CDS in df['predicted_CDS']]

    a = np.array((x[0],x[1]))
    dG_start_avg = np.mean(a,axis=0)
    scale = dG_start_INN['ATG'] - dG_start_avg[2]
    dG_start = dG_start_avg + scale
    '''

    # Having executed the preceeding block, we identified the following
    # adjusted free energies. Here, the free energies are scaled relative to
    # AUG, and the very positive values indicate IF-discrimination
    start_codon_energies = {
    'ATG': -2.76,
    'GTG': -0.42,
    'TTG': 1.81,
    'CTG': 7.09,
    'ATC': 7.18,
    'ATA': 8.19,
    'ATT': 13.23,
    'CAT': 13.41,
    'GGA': 16.89,
    'TGC': 17.99,
    'CGC': 18.17,
    'TAG': 19.14
    }

    # Let's try out the new parameters on both the NanoLuc and the sfGFP data
    filters = {'DATASET':['Hecht_NAR_2017']}
    models.add(RBSCalc_v2_1,start_energies=start_codon_energies)
    modelName = 'RBSCalc_v2_1'
    models.setform([modelName], x='dG_total', y='PROT.MEAN', yScale='ln', a1=-0.45)
    test = synbiomts.analyze.ModelTest(models,dbfilename,filters,add_data=True,verbose=True)
    test.run()

    with open("labels/labels1.txt","r") as f:
        predictLabels = [x.strip('\n') for x in f.readlines()]
    with open("labels/labels_stats.txt","r") as f:
        statsLabels = [x.strip('\n') for x in f.readlines()]
    test.to_excel('Hecht_Output',predictLabels,statsLabels)