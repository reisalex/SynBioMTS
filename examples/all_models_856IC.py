
import testsfm


# Load and wrap models
transl_rate_models = testsfm.interface.Models()

import RBS_Calculator_v1_0
def RBSCalc_v1(sequence,startpos,organism,temp):
    start_range = [0,startpos+1]
    model = RBS_Calculator_v1_0.RBS_Calculator(sequence,start_range,rRNA)
    model.temp = temp
    model.run()
    output = model.output()
    (RBS,best_start_pos) = find_best_start(sequence,start_pos,output)
    
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
    '5pUTR' : sequence[:,best_start_pos],
    'CDS' : sequence[best_start_pos:],
    'RBSobj': RBS
    }

    return results

import RBS_Calculator_v1_1
def RBSCalc_v1_1(sequence,startpos,organism,temp):
    start_range = [0,startpos+1]
    model = RBS_Calculator_v1_1.RBS_Calculator(sequence,start_range,rRNA)
    model.temp = temp
    model.run()
    output = model.output()
    (RBS,best_start_pos) = find_best_start(sequence,start_pos,output)

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
    '5pUTR' : sequence[:,best_start_pos],
    'CDS' : sequence[best_start_pos:],
    'RBSobj': RBS
    }

    return results

import RBS_Calculator_v2_0
def RBSCalc_v2(sequence,organism,temp,startpos):
    start_range = [0,startpos+1]
    model = RBS_Calculator_v2_0.RBS_Calculator(sequence,start_range,rRNA)
    model.temp = temp
    model.run()
    output = model.output()
    (RBS,best_start_pos) = find_best_start(sequence,start_pos,output)
    
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
    'warnings' : RBS.warnings,
    
    'blank' : "  ",
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
    
    'num_standby_site_hairpins' : RBS.num_standby_site_hairpins,
    'standby_site_hairpin_seqs' : str([seq for seq,bp,length,dG in RBS.standby_site_hairpins]).translate(None, "[]'u"),
    'standby_site_hairpin_num_bp' : str([bp for seq,bp,length,dG in RBS.standby_site_hairpins]).translate(None, "[]'u"),
    'standby_site_hairpin_seq_len' : str([length for seq,bp,length,dG in RBS.standby_site_hairpins]).translate(None, "[]'u"),
    'dG_standby_site_hairpin' : str([dG for seq,bp,length,dG in RBS.standby_site_hairpins]).translate(None, "[]'u"),
    
    'SD_hairpin_seqs' : str([seq for seq,bp,length,dG in RBS.SD_hairpins]).translate(None, "[]'u"),
    'SD_hairpin_num_bp' : [bp for seq,bp,length,dG in RBS.SD_hairpins][0],
    'SD_hairpin_seq_len' : [length for seq,bp,length,dG in RBS.SD_hairpins][0],
    'dG_SD_hairpin' : [dG for seq,bp,length,dG in RBS.SD_hairpins][0],
    
    'SD_start_hairpin_seqs' : str([seq for seq,bp,length,dG in RBS.SD_start_hairpins]).translate(None, "[]'u"),
    'SD_start_hairpin_num_bp' : [bp for seq,bp,length,dG in RBS.SD_start_hairpins][0],
    'SD_start_hairpin_seq_len' : [length for seq,bp,length,dG in RBS.SD_start_hairpins][0],
    'dG_SD_start_hairpin' : [dG for seq,bp,length,dG in RBS.SD_start_hairpins][0],
    
    'start_hairpin_seqs' : str([seq for seq,bp,length,dG in RBS.start_hairpins]).translate(None, "[]'u"),
    'start_hairpin_num_bp' : [bp for seq,bp,length,dG in RBS.start_hairpins][0],
    'start_hairpin_seq_len' : [length for seq,bp,length,dG in RBS.start_hairpins][0],
    'dG_start_hairpin' : [dG for seq,bp,length,dG in RBS.start_hairpins][0],
    
    'footprint_hairpin_seqs' : str([seq for seq,bp,length,dG in RBS.footprint_hairpins]).translate(None, "[]'u"),
    'footprint_hairpin_num_bp' : [bp for seq,bp,length,dG in RBS.footprint_hairpins][0],
    'footprint_hairpin_seq_len' : [length for seq,bp,length,dG in RBS.footprint_hairpins][0],
    'dG_footprint_hairpin' : [dG for seq,bp,length,dG in RBS.footprint_hairpins][0],

    'post_footprint_hairpin_seqs' : str([seq for seq,bp,length,dG in RBS.post_footprint_hairpins]).translate(None, "[]'u"),
    'post_footprint_hairpin_num_bp' : [bp for seq,bp,length,dG in RBS.post_footprint_hairpins][0],
    'post_footprint_hairpin_seq_len' : [length for seq,bp,length,dG in RBS.post_footprint_hairpins][0],
    'dG_post_footprint_hairpin' : [dG for seq,bp,length,dG in RBS.post_footprint_hairpins][0],

    '5pUTR' : sequence[:,best_start_pos],
    'CDS' : sequence[best_start_pos:],
    'RBSobj': RBS
    }

    return results

import UTR_Designer

handle = open('models/RBSDesigner_ALL.p','r')

from emopec._emopec import predict_spacing, get_expression


        
    elif model == "EMOPEC":
    
        preseq,SD,spacer,Expression = RBS
        
        
        maxExpression = get_expression(sd_seq='AGGAGA',sd_dist=len(spacer))
        minExpression = get_expression(sd_seq='TTGGGC',sd_dist=len(spacer))
        
        Expression_percent = (Expression - minExpression)/(maxExpression - minExpression)
        
        results_dict = {
        'preseq' : preseq,
        'SD_seq' : SD,
        'spacer_seq' : spacer,
        'spacer_length' : len(spacer),
        'Expression' : Expression,
        'Expression_percent': Expression_percent,
        }
    
    elif model == "RBS Designer":
        pass
        #RBS['translation efficiency']
        #RBS['Expression']
        #RBS['exposure probability']
        #RBS['mRNA-ribosome probability']
        #RBS['SD']
        #RBS['dG_SD:aSD']
        #RBS['spacer length']
    
    #=== All models ===#
    
    results_dict['fivepUTR'] = mRNA[:start_pos]
    results_dict['CDS'] = mRNA[start_pos:]

    return results_dict




transl_rate_models.register("RBSCalc_v1",RBSCalc_v1)
transl_rate_models.register("RBSCalc_v1_1",RBSCalc_v1_1)
transl_rate_models.register("RBSCalc_v2",RBSCalc_v2)







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
fileName = '../geneticsystems.db'

customtest = testsfm.analyze.ModelTest(transl_rate_models,datasets,fileName)
customtest.run()