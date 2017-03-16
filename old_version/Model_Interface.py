import sys,os,re
import pickle

sys.path.append('models')
from RBS_Calculator_v1_0 import RBS_Calculator as RBS_Calculator_v1_0
from RBS_Calculator_v1_1 import RBS_Calculator as RBS_Calculator_v1_1
from RBS_Calculator_v2_0 import RBS_Calculator as RBS_Calculator_v2_0
from RBS_Calculator_v2_1 import RBS_Calculator as RBS_Calculator_v2_1
#from RBS_Calc_v2_UTR_Designer import RBS_Calculator as RBS_Calc_v2_UTR_Designer
from UTR_Designer import RBS_Calculator as UTR_Designer
from emopec._emopec import predict_spacing, get_expression

handle = open('models/RBSDesigner_ALL.p','r')
RBSDesigner = pickle.load(handle)
handle.close()

toolbox = {'RBS Calculator v1.0': RBS_Calculator_v1_0,
           'RBS Calculator v1.1': RBS_Calculator_v1_1,
           'RBS Calculator v2.0': RBS_Calculator_v2_0,
           'RBS Calculator v2.1': RBS_Calculator_v2_1,
           #'RBS_Calc_v2_UTR_Designer': RBS_Calc_v2_UTR_Designer,
           'UTR Designer':        UTR_Designer,
           'EMOPEC':              predict_spacing,
           'RBS Designer':        RBSDesigner}

def _call_model(data):
    
    paper       = data['paper']
    subset      = data['subset']
    indx        = data['indx']
    
    model       = data['model']
    mRNA        = data['mRNA']
    rRNA        = data['rRNA']
    
    start_pos   = data['start_pos']
    T = data['temperature']
    
    # print paper,"\t",subset,"\t",indx,"\t",model
    
    if model in ["RBS Calculator v1.0","UTR Designer","RBS Calculator v1.1","RBS Calculator v2.0","RBS Calculator v2.1","RBS_Calc_v2_UTR_Designer"]:
        
        start_range = [0,start_pos+1]
        m = toolbox[model](mRNA,start_range,rRNA)
        m.temp = T
        
        # Testing hypotheses
        #m.optimal_spacing = 8
        #m.footprint = 34
        #m.fixed_mRNA_post_cutoff=100
        #m.dangle = 'all'
        #m.Use_centroid = True
        #m.RNA_model = 'rna_andronescu2007'
        
        #print m.RNA_model
        
        m.run()
        output = m.output()
        
        (RBS,best_start_pos) = _find_best_start(mRNA,start_pos,output)
        results_dict = _get_output(model,mRNA,RBS,best_start_pos)
        
    elif model == "EMOPEC":
        
        RBS = toolbox[model](mRNA[0:start_pos])
        results_dict = _get_output(model,mRNA,RBS,start_pos)
        
    elif model == "RBS Designer":
        # Uses RBSDesigner results calculated in Windows
        
        if mRNA in toolbox[model].keys(): RBS = toolbox[model][mRNA]
        else:
            RBS = {'translation efficiency': -1,
                    'Expression': -1,
                    'exposure probability': -1,
                    'mRNA-ribosome probability': -1,
                    'SD': -1,
                    'dG_SD:aSD': -1,
                    'spacer length': -1}
        
        results_dict = RBS
        
    return results_dict


def _get_output(model,mRNA,RBS,start_pos):

    if model in ["RBS Calculator v1.0","RBS Calculator v1.1"]:
        
        results_dict = {
        'TIR' : RBS.tir,
        'dG_total' : RBS.dG_total,
        'dG_mRNA_rRNA' : RBS.dG_mRNA_rRNA,
        'dG_mRNA' : RBS.dG_mRNA,
        'dG_start' : RBS.dG_start,
        'dG_standby' : RBS.dG_standby,
        'dG_spacing' : RBS.dG_spacing,
        'used_mRNA_sequence' : RBS.sequence,
        'warnings' : RBS.warnings,
        'RBSobj': RBS
        }
        
    elif model == "UTR Designer":
        
        results_dict = {
        'TIR' : RBS.tir,
        'dG_total' : RBS.dG_total,
        'dG_SD': RBS.dG_SD_aSD,
        'dG_start' : RBS.dG_start,
        'dG_spacing' : RBS.dG_spacing,
        'dG_direct': RBS.dG_direct,
        'dG_indirect': RBS.dG_indirect,
        'used_mRNA_sequence' : RBS.sequence,
        'warnings' : RBS.warnings,
        'RBSobj': RBS
        }
        
    elif model in ["RBS Calculator v2.0","RBS Calculator v2.1"]:
        
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
        
        #'DB_struc': RBS.DB_struc,
        #'aDB_struc': RBS.aDB_struc,
        #'dG_DB_aDB': RBS.dG_DB_aDB,
        
        'RBSobj': RBS
        }
    
    elif model == "RBS_Calc_v2_UTR_Designer":
        
        results_dict = {
        'TIR' : RBS.tir,
        'dG_total' : RBS.dG_total,
        'dG_mRNA_rRNA': RBS.dG_mRNA_rRNA,
        'dG_standby' : RBS.dG_standby,
        'dG_start' : RBS.dG_start,
        'dG_spacing' : RBS.dG_spacing,
        'dG_direct': RBS.dG_direct,
        'dG_indirect': RBS.dG_indirect,
        'used_mRNA_sequence' : RBS.sequence,
        'warnings' : RBS.warnings,
        'RBSobj': RBS
        }        
        
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
    

def _find_best_start(mRNA, start_pos, predictions):
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
    
    #fivepUTR = mRNA[:best_start_pos]
    #CDS = mRNA[best_start_pos:] 
    
    return (RBS,best_start_pos)
