import os,sys
import xlrd
from math import *
import numpy as np
import random

rRNA_dic = {'Escherichia coli str. K-12 substr. DH10B' : 'ACCTCCTTA',
            'Escherichia coli str. K-12 substr. MG1655' : 'ACCTCCTTA',
            'Escherichia coli BL21(DE3)' : 'ACCTCCTTA',
            'Salmonella enterica subsp. enterica serovar Typhimurium str. LT2' : 'ACCTCCTTA',
            'Corynebacterium glutamicum B-2784' : 'ACCTCCTTT',
            'Pseudomonas fluorescens A506' : 'ACCTCCTTT',
            'Bacillus subtilis subsp. subtilis str. 168' : 'ACCTCCTTT',
            'Bacteroides thetaiotaomicron VPI-5482' : 'ACCTCCTTT'}

start_codons = ['ATG','GTG','CTG','TTG']

def _import(db,datasets):
    
    # ==================================================================
    paper = 'EspahBorujeni_NAR_2013'
    
    if paper in datasets:
        
        path = _get_filepath(paper)
        wb = xlrd.open_workbook(path,'rb')
        sheet = wb.sheet_by_index(0)
        db[paper] = {}
        
        db[paper]['authors'] = "Amin Espah Borujeni, Anirudh S. Channarasappa, & Howard M. Salis"
        db[paper]['title'] = "Translation rate is controlled by coupled trade-offs between site accessibility, select RNA unfolding and sliding at upstream standby sites"
        db[paper]['journal'] = "Nucleic Acids Research, 2014, Vol. 24, No. 4; doi: 10.1093/nar/gkt1139"
        db[paper]['abbrev'] = "StandbySite"
        
        db[paper]['xls_ID']    = sheet.col_values(colx=0, start_rowx=5, end_rowx=141)
        db[paper]['fivepUTR']  = sheet.col_values(colx=4, start_rowx=5, end_rowx=141)
        db[paper]['CDS']       = sheet.col_values(colx=5, start_rowx=5, end_rowx=141)
        db[paper]['fluo_mean'] = sheet.col_values(colx=7, start_rowx=5, end_rowx=141)
        db[paper]['fluo_std']  = sheet.col_values(colx=8, start_rowx=5, end_rowx=141)
    
        #=== Added extended dataset for StandbySite sequences ===#

        paper_ext = 'EspahBorujeni_NAR_2013_extended'
        path = _get_filepath(paper_ext)
        wb = xlrd.open_workbook(path,'rb')
        sheet = wb.sheet_by_index(0)
        
        db[paper]['xls_ID']    += sheet.col_values(colx=0, start_rowx=3, end_rowx=41)
        db[paper]['fivepUTR']  += sheet.col_values(colx=4, start_rowx=3, end_rowx=41)
        db[paper]['CDS']       += sheet.col_values(colx=5, start_rowx=3, end_rowx=41)
        db[paper]['fluo_mean'] += sheet.col_values(colx=9, start_rowx=3, end_rowx=41)
        db[paper]['fluo_std']  += sheet.col_values(colx=10, start_rowx=3, end_rowx=41)
        
        num_seqs = len(db[paper]['xls_ID'])
        
        db[paper]['num_seqs'] = num_seqs
        db[paper]['temp']     = [37.0]*num_seqs
        db[paper]['protein']  = ['RFP']*num_seqs
        db[paper]['organism'] = ["Escherichia coli str. K-12 substr. DH10B"]*num_seqs
        db[paper]['info']     = ["Escherichia coli str. K-12 substr. DH10B expressing RFP"]*num_seqs
        
        db[paper]['mRNA'] = [fivepUTR+CDS for fivepUTR,CDS in zip(db[paper]['fivepUTR'],db[paper]['CDS'])]
        db[paper]['rRNA'] = [rRNA_dic[organism] for organism in db[paper]['organism']]
        db[paper]['start_pos'] = [len(fivepUTR) for fivepUTR in db[paper]['fivepUTR']]
        
        db[paper]['Expression'] = db[paper]['fluo_mean']
        
        db = _format_ascii(db,paper)
        db = _subdivide_dataset(db,paper)
        
    # ==================================================================
    paper = 'EspahBorujeni_NAR_2015'    
    
    if paper in datasets:
        
        path = _get_filepath(paper)
        wb = xlrd.open_workbook(path,'rb')
        sheet = wb.sheet_by_index(0)
        db[paper] = {}
        
        db[paper]['authors'] = "Amin Espah Borujeni, Bennis M. Mishler, Jingzhi Wang, Walker Huso, & Howard M. Salis"
        db[paper]['title'] = "Automated Physics-Based Design of Synthetic Riboswitches from Diverse RNA Aptamers"
        db[paper]['journal'] = "Nucleic Acids Research, 2015, doi: 10.1093/nar/gkv1289"
        db[paper]['abbrev'] = "Riboswitch"
        
        db[paper]['xls_ID']    = sheet.col_values(colx=0, start_rowx=3, end_rowx=75)
        db[paper]['pre_aptamer']  = sheet.col_values(colx=3, start_rowx=3, end_rowx=75)
        db[paper]['aptamer']      = sheet.col_values(colx=4, start_rowx=3, end_rowx=75)
        db[paper]['post_aptamer'] = sheet.col_values(colx=5, start_rowx=3, end_rowx=75)
        
        db[paper]['fivepUTR']  = ["{}{}{}".format(pre,aptamer,post) for pre,aptamer,post \
                                  in zip(db[paper]['pre_aptamer'],db[paper]['aptamer'],db[paper]['post_aptamer'])]
        
        db[paper]['CDS']       = sheet.col_values(colx=6, start_rowx=3, end_rowx=75)
        db[paper]['fluo_mean'] = sheet.col_values(colx=16, start_rowx=3, end_rowx=75)
        db[paper]['fluo_std']  = sheet.col_values(colx=17, start_rowx=3, end_rowx=75)
        
        num_seqs = len(db[paper]['xls_ID'])
        
        db[paper]['num_seqs'] = num_seqs
        db[paper]['temp']     = [37.0]*num_seqs
        db[paper]['protein']  = sheet.col_values(colx=1, start_rowx=3, end_rowx=75)
        db[paper]['organism'] = ["Escherichia coli str. K-12 substr. DH10B"]*num_seqs
        db[paper]['info']     = ["{}_{}".format(host,protein) for host,protein \
                                 in zip(db[paper]['organism'],db[paper]['protein'])]

        db[paper]['mRNA'] = [fivepUTR+CDS for fivepUTR,CDS in zip(db[paper]['fivepUTR'],db[paper]['CDS'])]
        db[paper]['rRNA'] = [rRNA_dic[organism] for organism in db[paper]['organism']]
        db[paper]['start_pos'] = [len(fivepUTR) for fivepUTR in db[paper]['fivepUTR']]
        
        db[paper]['Expression'] = db[paper]['fluo_mean']
        
        db = _format_ascii(db,paper)
        db = _subdivide_dataset(db,paper)
        
    # ==================================================================
    paper = 'EspahBorujeni_JACS_2016'    
    
    if paper in datasets:
        
        path = _get_filepath(paper)
        wb = xlrd.open_workbook(path,'rb')
        sheet = wb.sheet_by_index(0)
        db[paper] = {}
        
        db[paper]['authors'] = "Amin Espah Borujeni, Howard M. Salis"
        db[paper]['title'] = "Translation Initiation is Controlled by RNA Folding Kinetics via a Ribosome Drafting Mechanism"
        db[paper]['journal'] = "JACS"
        db[paper]['abbrev'] = "RibosomeDrafting"
        
        db[paper]['xls_ID']    = sheet.col_values(colx=0, start_rowx=6, end_rowx=42)
        db[paper]['fivepUTR']  = sheet.col_values(colx=3, start_rowx=6, end_rowx=42)
        db[paper]['CDS']       = sheet.col_values(colx=4, start_rowx=6, end_rowx=42)
        db[paper]['fluo_mean'] = sheet.col_values(colx=9, start_rowx=6, end_rowx=42)
        db[paper]['fluo_std']  = sheet.col_values(colx=10, start_rowx=6, end_rowx=42)
        
        num_seqs = len(db[paper]['xls_ID'])
        
        db[paper]['num_seqs'] = num_seqs
        db[paper]['temp']     = [37.0]*num_seqs
        db[paper]['protein']  = ['RFP']*num_seqs
        db[paper]['organism'] = ["Escherichia coli str. K-12 substr. DH10B"]*num_seqs
        db[paper]['info']     = ["Escherichia coli str. K-12 substr. DH10B expressing RFP"]*num_seqs

        db[paper]['mRNA'] = [fivepUTR+CDS for fivepUTR,CDS in zip(db[paper]['fivepUTR'],db[paper]['CDS'])]
        db[paper]['rRNA'] = [rRNA_dic[organism] for organism in db[paper]['organism']]
        db[paper]['start_pos'] = [len(fivepUTR) for fivepUTR in db[paper]['fivepUTR']]
        
        db[paper]['Expression'] = db[paper]['fluo_mean']
        
        db = _format_ascii(db,paper)
        db = _subdivide_dataset(db,paper)    
    
    # ==================================================================
    paper = 'EspahBorujeni_Footprint'    
    
    if paper in datasets:
        
        path = _get_filepath(paper)
        wb = xlrd.open_workbook(path,'rb')
        sheet = wb.sheet_by_index(0)
        db[paper] = {}
        
        db[paper]['authors'] = "Amin Espah Borujeni, Daniel P. Cetnar, Howard M. Salis"
        db[paper]['title'] = "Precise Quantification of Translation Inhibition by RNA structures that Overlap with the Ribosome Footprint at N-terminal Coding Sections"
        db[paper]['journal'] = "Nucleic Acids Research"
        db[paper]['abbrev'] = "Footprint"
        
        db[paper]['xls_ID']    = sheet.col_values(colx=0, start_rowx=5, end_rowx=32)
        db[paper]['fivepUTR']  = sheet.col_values(colx=3, start_rowx=5, end_rowx=32)
        db[paper]['CDS']       = sheet.col_values(colx=4, start_rowx=5, end_rowx=32)
        db[paper]['fluo_mean'] = sheet.col_values(colx=10, start_rowx=5, end_rowx=32)
        db[paper]['fluo_std']  = sheet.col_values(colx=11, start_rowx=5, end_rowx=32)
        
        num_seqs = len(db[paper]['xls_ID'])
        
        db[paper]['num_seqs'] = num_seqs
        db[paper]['temp']     = [37.0]*num_seqs
        db[paper]['protein']  = ['RFP']*num_seqs
        db[paper]['organism'] = ["Escherichia coli str. K-12 substr. DH10B"]*num_seqs
        db[paper]['info']     = ["Escherichia coli str. K-12 substr. DH10B expressing RFP"]*num_seqs        
        
        db[paper]['mRNA'] = [fivepUTR+CDS for fivepUTR,CDS in zip(db[paper]['fivepUTR'],db[paper]['CDS'])]
        db[paper]['rRNA'] = [rRNA_dic[organism] for organism in db[paper]['organism']]
        db[paper]['start_pos'] = [len(fivepUTR) for fivepUTR in db[paper]['fivepUTR']]
        
        db[paper]['Expression'] = db[paper]['fluo_mean']
        
        db = _format_ascii(db,paper)
        db = _subdivide_dataset(db,paper)
    
    # ==================================================================
    paper = 'Salis_Nat_Biotech_2009'
    
    if paper in datasets:
    
        path = _get_filepath(paper)
        wb = xlrd.open_workbook(path,'rb')
        sheet = wb.sheet_by_index(0)
        db[paper] = {}
        
        db[paper]['authors'] = "Howard M. Salis, Ethan A. Mirsky, & Christopher A. Voigt"
        db[paper]['title'] = "Automated design of synthetic ribosome binding sites to control protein expression"
        db[paper]['journal'] = "Nature Biotechnology, 2009, Vol. 27, No. 10; doi: 10.1038/nbt.1568"
        db[paper]['abbrev'] = "Salis2009"
        
        db[paper]['xls_ID']    = sheet.col_values(colx=0, start_rowx=3, end_rowx=135)
        db[paper]['excel_seq']  = sheet.col_values(colx=3, start_rowx=3, end_rowx=135)
        
        SacI = 'GAGCTC'
        SacI_positions = [seq.find(SacI) for seq in db[paper]['excel_seq']]
        start_positions = [p - 5 for p in SacI_positions]
        
        # Find correct starts
        for i in range(len(db[paper]['excel_seq'])):
            seq = db[paper]['excel_seq'][i]
            start_pos = start_positions[i]
            start = seq[start_pos:start_pos+3]
            while start not in start_codons:
                start_pos -= 3
                start = seq[start_pos:start_pos+3]
            start_positions[i] = start_pos
        
        RFP1_CDS = sheet.cell_value(1,3)
        
        db[paper]['fivepUTR'] = [seq[0:start_pos] for seq,start_pos in \
                                zip(db[paper]['excel_seq'],start_positions)]
        
        db[paper]['CDS']      = [seq[start_pos:SacI_pos] + RFP1_CDS for seq,start_pos,SacI_pos in \
                                zip(db[paper]['excel_seq'],start_positions,SacI_positions)]        
        
        num_seqs = len(db[paper]['xls_ID'])
        
        db[paper]['num_seqs'] = num_seqs
        db[paper]['fluo_mean'] = sheet.col_values(colx=4, start_rowx=3, end_rowx=135)
        db[paper]['fluo_std']  = sheet.col_values(colx=5, start_rowx=3, end_rowx=135)
        db[paper]['temp']      = [37.0]*num_seqs
        db[paper]['protein']   = ['RFP']*num_seqs
        db[paper]['organism']  = ["Escherichia coli str. K-12 substr. DH10B"]*num_seqs
        db[paper]['info']      = ["Escherichia coli str. K-12 substr. DH10B expressing RFP"]*num_seqs         
        
        db[paper]['mRNA'] = [fivepUTR+CDS for fivepUTR,CDS in zip(db[paper]['fivepUTR'],db[paper]['CDS'])]
        db[paper]['rRNA'] = [rRNA_dic[organism] for organism in db[paper]['organism']]
        db[paper]['start_pos'] = [len(fivepUTR) for fivepUTR in db[paper]['fivepUTR']]
        
        db[paper]['Expression'] = db[paper]['fluo_mean']
        
        db = _format_ascii(db,paper)
        db = _subdivide_dataset(db,paper)
        
    # ==================================================================
    paper = 'Farasat_MSB_2014'
    
    if paper in datasets:
    
        path = _get_filepath(paper)
        wb = xlrd.open_workbook(path,'rb')
        sheet = wb.sheet_by_index(0)
        db[paper] = {}
        
        db[paper]['authors'] = "Iman Farasat, Manish Kushwaha, Jason Collens, Michael Easterbrook, Matthew Guido, & Howard M. Salis"
        db[paper]['title'] = "Efficient search, mapping, and optimization of multi-protein genetic systems in diverse bacteria"
        db[paper]['journal'] = "Molecular Systems Biology. 2014 Jun 21; 10:731. doi: 10.15252/msb.20134955"
        db[paper]['abbrev'] = "dRBSlibraries"
        
        db[paper]['xls_ID']    = sheet.col_values(colx=0, start_rowx=1, end_rowx=146)
        db[paper]['info']      = sheet.col_values(colx=1, start_rowx=1, end_rowx=146)
        db[paper]['organism']  = sheet.col_values(colx=2, start_rowx=1, end_rowx=146)
        db[paper]['protein']   = sheet.col_values(colx=3, start_rowx=1, end_rowx=146)
        db[paper]['preseq']    = sheet.col_values(colx=5, start_rowx=1, end_rowx=146)
        db[paper]['RBS']       = sheet.col_values(colx=6, start_rowx=1, end_rowx=146)
        db[paper]['CDS']       = sheet.col_values(colx=7, start_rowx=1, end_rowx=146)
        db[paper]['fluo_mean'] = sheet.col_values(colx=9, start_rowx=1, end_rowx=146)
        db[paper]['fluo_std']  = sheet.col_values(colx=10, start_rowx=1, end_rowx=146)
        
        db[paper]['fivepUTR']  = ["{}{}".format(preseq,RBS) for preseq,RBS \
                                  in zip(db[paper]['preseq'],db[paper]['RBS'])]

        db[paper]['num_seqs'] = len(db[paper]['xls_ID'])
        db[paper]['temp']     = [37.0]*db[paper]['num_seqs']
        
        db[paper]['mRNA'] = [fivepUTR+CDS for fivepUTR,CDS in zip(db[paper]['fivepUTR'],db[paper]['CDS'])]
        db[paper]['rRNA'] = [rRNA_dic[organism] for organism in db[paper]['organism']]
        db[paper]['start_pos'] = [len(fivepUTR) for fivepUTR in db[paper]['fivepUTR']]
        
        db[paper]['Expression'] = db[paper]['fluo_mean']
        
        db = _format_ascii(db,paper)
        db = _subdivide_dataset(db,paper)
    
    # ==================================================================
    paper = 'Tian_NAR_2015'
    
    if paper in datasets:
    
        path = _get_filepath(paper)
        wb = xlrd.open_workbook(path,'rb')
        sheet = wb.sheet_by_index(0)
        db[paper] = {}
        
        db[paper]['authors'] = "Tian Tian, & Howard M. Salis"
        db[paper]['title'] = "A Predictive Biophysical Model of Translational Coupling to Coordinate and Control Protein expression in Bacterial Operons"
        db[paper]['journal'] = "Nucleic Acids Research, 2015 doi: 10.1093/nar/gkv635"
        db[paper]['abbrev'] = "TranslationalCoupling"
        
        db[paper]['xls_ID']    = sheet.col_values(colx=0, start_rowx=2, end_rowx=26)
        db[paper]['info']      = sheet.col_values(colx=1, start_rowx=2, end_rowx=26)
        db[paper]['organism']  = sheet.col_values(colx=2, start_rowx=2, end_rowx=26)
        db[paper]['protein']   = sheet.col_values(colx=3, start_rowx=2, end_rowx=26)
        db[paper]['fivepUTR']  = sheet.col_values(colx=4, start_rowx=2, end_rowx=26)
        db[paper]['CDS']       = sheet.col_values(colx=5, start_rowx=2, end_rowx=26)
        db[paper]['fluo_mean'] = sheet.col_values(colx=12, start_rowx=2, end_rowx=26)
        db[paper]['fluo_std']  = sheet.col_values(colx=13, start_rowx=2, end_rowx=26)
        db[paper]['num_seqs'] = len(db[paper]['xls_ID'])
        db[paper]['temp']     = [37.0]*db[paper]['num_seqs']
        
        db[paper]['mRNA'] = [fivepUTR+CDS for fivepUTR,CDS in zip(db[paper]['fivepUTR'],db[paper]['CDS'])]
        db[paper]['rRNA'] = [rRNA_dic[organism] for organism in db[paper]['organism']]
        db[paper]['start_pos'] = [len(fivepUTR) for fivepUTR in db[paper]['fivepUTR']]
        
        db[paper]['Expression'] = db[paper]['fluo_mean']
        
        db = _format_ascii(db,paper)
        db = _subdivide_dataset(db,paper)
    
    # ==================================================================
    paper = 'Mimee_Cell_Sys_2015'
    
    if paper in datasets:
    
        path = _get_filepath(paper)
        wb = xlrd.open_workbook(path,'rb')
        sheet = wb.sheet_by_index(0)
        db[paper] = {}
        
        db[paper]['authors'] = "Mark Mimee, Alex C. Tucker, Christopher A. Voigt, and Timothy K. Lu"
        db[paper]['title'] = "Programming a Human Commensal Bacterium, Bacteroides thetaiotaomicron, to Sense and Respond to Stimuli in the Murine Gut Microbiota"
        db[paper]['journal'] = "Cell Systems 1, 62-71, July 29, 2015"
        db[paper]['abbrev'] = "Bthetaiotaomicron"
        
        db[paper]['xls_ID']    = sheet.col_values(colx=3, start_rowx=3, end_rowx=146)
        db[paper]['info']      = sheet.col_values(colx=0, start_rowx=3, end_rowx=146)
        db[paper]['organism']  = sheet.col_values(colx=1, start_rowx=3, end_rowx=146)
        db[paper]['protein']   = sheet.col_values(colx=2, start_rowx=3, end_rowx=146)
        db[paper]['fivepUTR']  = sheet.col_values(colx=6, start_rowx=3, end_rowx=146)
        db[paper]['CDS']       = sheet.col_values(colx=7, start_rowx=3, end_rowx=146)
        db[paper]['fluo_mean'] = sheet.col_values(colx=8, start_rowx=3, end_rowx=146)
        db[paper]['fluo_std']  = sheet.col_values(colx=9, start_rowx=3, end_rowx=146)
        db[paper]['num_seqs'] = len(db[paper]['xls_ID'])
        db[paper]['temp']     = [37.0]*db[paper]['num_seqs']
        
        db[paper]['mRNA'] = [fivepUTR+CDS for fivepUTR,CDS in zip(db[paper]['fivepUTR'],db[paper]['CDS'])]
        db[paper]['rRNA'] = [rRNA_dic[organism] for organism in db[paper]['organism']]
        db[paper]['start_pos'] = [len(fivepUTR) for fivepUTR in db[paper]['fivepUTR']]
        
        db[paper]['Expression'] = db[paper]['fluo_mean']
        
        db = _format_ascii(db,paper)
        db = _subdivide_dataset(db,paper)
    
    # ==================================================================
    paper = 'Bonde_NatMethods_IC_2016'
    # Only including the IC data for now from Bonde et al., Nat Methods
    
    if paper in datasets:
    
        path = _get_filepath(paper)
        wb = xlrd.open_workbook(path,'rb')
        sheet = wb.sheet_by_index(0)
        db[paper] = {}
        
        db[paper]['authors'] = "Mads T Bonde, Margit Pederse, Michael S Klausen, Sheila I Jensen, Tune Wulff, Scott Harrison, Alex T Nielsen, Markus J Herrgard, Morten O A Sommer"
        db[paper]['title'] = "Predictable tuning of protein expression in bacteria"
        db[paper]['journal'] = "Nature Methods, 13: 3, 233-236, 2016"
        db[paper]['abbrev'] = "EMOPEC"
        
        db[paper]['xls_ID']    = sheet.col_values(colx=0, start_rowx=1, end_rowx=107)
        db[paper]['info']      = sheet.col_values(colx=1, start_rowx=1, end_rowx=107)
        db[paper]['organism']  = sheet.col_values(colx=2, start_rowx=1, end_rowx=107)
        db[paper]['protein']   = sheet.col_values(colx=3, start_rowx=1, end_rowx=107)
        db[paper]['fivepUTR']  = sheet.col_values(colx=5, start_rowx=1, end_rowx=107)
        db[paper]['CDS']       = sheet.col_values(colx=6, start_rowx=1, end_rowx=107)
        db[paper]['fluo_mean'] = sheet.col_values(colx=11, start_rowx=1, end_rowx=107)
        db[paper]['num_seqs'] = len(db[paper]['xls_ID'])
        db[paper]['fluo_std']  = [0.0]*db[paper]['num_seqs']
        db[paper]['temp']     = [37.0]*db[paper]['num_seqs']
        
        db[paper]['mRNA'] = [fivepUTR+CDS for fivepUTR,CDS in zip(db[paper]['fivepUTR'],db[paper]['CDS'])]
        db[paper]['rRNA'] = [rRNA_dic[organism] for organism in db[paper]['organism']]
        db[paper]['start_pos'] = [len(fivepUTR) for fivepUTR in db[paper]['fivepUTR']]
        
        db[paper]['Expression'] = db[paper]['fluo_mean']
        
        db = _format_ascii(db,paper)
        db = _subdivide_dataset(db,paper)
    
    # ==================================================================
    paper = 'Kosuri_PNAS_2013'
    
    if paper in datasets:
        
        path = _get_filepath(paper)
        wb = xlrd.open_workbook(path,'rb')
        data = {}
        
        data['authors'] = "Sriram Kosuri, Daniel B. Goodman, George M. Church"
        data['title'] = "Composability of regulatory sequences controlling transcription and translation in Escherichia coli"
        data['journal'] = "Proc Natl Acad Sci USA, 2013, Vol. 110 no. 34"
        data['abbrev'] = "Kosuri2013"        
        
        sheet = wb.sheet_by_name("Promoters")
        promoter_dict = {}
        for i in range(1,113):
            promoter = sheet.cell_value(i,0)
            TSS = int(sheet.cell_value(i,7))
            promoter_dict[promoter] = {}
            promoter_dict[promoter]['sequence'] = sheet.cell_value(i,9).replace(' ','').replace('"','')
            promoter_dict[promoter]["TSS"] = TSS
        
        sheet = wb.sheet_by_name("RBSs")
        RBS_dict = {}
        for i in range(1,112):
            RBS = sheet.cell_value(i,0)
            RBS_dict[RBS] = (sheet.cell_value(i,9).replace(' ','').replace('"',''))[:-3]
        
        sheet = wb.sheet_by_index(0)
        
        data['xls_ID'] = sheet.col_values(colx=2, start_rowx=1, end_rowx=9337)
        
        data['promoter_ID'] = sheet.col_values(colx=0, start_rowx=1, end_rowx=9337)
        data['promoter']    = [promoter_dict[ID]['sequence'] for ID in data['promoter_ID']]
        data['TSS']         = [promoter_dict[ID]['TSS'] for ID in data['promoter_ID']]
        data['RBS_ID']      = sheet.col_values(colx=1, start_rowx=1, end_rowx=9337)
        data['RBS']         = [RBS_dict[ID] for ID in data['RBS_ID']]
        data['CDS']         = sheet.col_values(colx=85, start_rowx=1, end_rowx=9337)
        
        data['fivepUTR']  = [p[TSS:]+RBS for p,TSS,RBS \
                                  in zip(data['promoter'],data['TSS'],data['RBS'])]
        
        num_seqs = len(data['xls_ID'])
        
        data['num_seqs'] = num_seqs
        data['temp']     = [30.0]*num_seqs        
        data['gene']  = ['sfGFP']*num_seqs
        data['organism'] = ['Escherichia coli str. K-12 substr. MG1655']*num_seqs
        data['info']     = ["{} expressing {}".format(host,protein) for host,protein \
                                 in zip(data['organism'],data['gene'])]
        
        data['mRNA'] = [fivepUTR+CDS for fivepUTR,CDS in zip(data['fivepUTR'],data['CDS'])]
        data['rRNA'] = [rRNA_dic[organism] for organism in data['organism']]
        data['start_pos'] = [len(fivepUTR) for fivepUTR in data['fivepUTR']]
        
        x = np.array(sheet.col_values(colx=3, start_rowx=1, end_rowx=9337))
        data['Count.Protein'] = x.astype(float)
        
        data['Count.A.RNA'] = np.array(sheet.col_values(colx=31, start_rowx=1, end_rowx=9337))
        data['Count.B.RNA'] = np.array(sheet.col_values(colx=32, start_rowx=1, end_rowx=9337))
        data['Count.RNA']   = np.array(sheet.col_values(colx=33, start_rowx=1, end_rowx=9337))
        data['Count.A.DNA'] = np.array(sheet.col_values(colx=34, start_rowx=1, end_rowx=9337))
        data['Count.B.DNA'] = np.array(sheet.col_values(colx=35, start_rowx=1, end_rowx=9337))
        data['Count.DNA']   = np.array(sheet.col_values(colx=36, start_rowx=1, end_rowx=9337))
        
        sheet = wb.sheet_by_name("NGS counts")
        data['total.rna.a'] = sheet.cell_value(0,1)
        data['total.rna.b'] = sheet.cell_value(1,1)
        data['total.dna.a'] = sheet.cell_value(2,1)
        data['total.dna.b'] = sheet.cell_value(3,1)
        
        bin_list = ['bin'+str(num) for num in range(1,13)]
        data['bin_list'] = bin_list    
        
        sheet = wb.sheet_by_name("bin_vals")
        flow_seq = {b: {} for b in bin_list}
        for b,col in zip(bin_list,range(1,13)):
            flow_seq[b]['cell_fraction'] = sheet.cell_value(2,col)
            flow_seq[b]['fluo'] = sheet.cell_value(4,col)
            #left_edge = sheet.cell_value(3,col)
            #right_edge = sheet.cell_value(3,col+1)
            #flow_seq[b]['fluo'] = np.exp((np.log(left_edge) + np.log(right_edge))/2.0)
        data['bins'] = flow_seq
        
        # read contig counts from Excel
        sheet = wb.sheet_by_index(0)
        for b,col in zip(bin_list,range(4,16)):
            x = np.array(sheet.col_values(colx=col, start_rowx=1, end_rowx=9337))
            data[b] = x.astype(float)   
        
        data = run_FS_calcs(data)        
        data.pop('bins')
        data.pop('bin_list')
        
        db[paper] = {}
        db[paper] = data        
        db = _format_ascii(db,paper)
        db = _subdivide_dataset(db,paper)
    
    # ==================================================================
    paper = 'Goodman_Science_2013'
    
    if paper in datasets:
    
        path = _get_filepath(paper)
        wb = xlrd.open_workbook(path,'rb')
        
        data = {}
        data['authors'] = "Daniel B. Goodman, George M. Church, Sriram Kosuri"
        data['title'] = "Causes and Effects of N-Terminal Codon Bias in Bacterial Genes"
        data['journal'] = "Science, 2013, Vol. 342"
        data['abbrev'] = "Goodman2013"      

        sheet = wb.sheet_by_name("Goodman_data")
        
        data['xls_ID']   = sheet.col_values(colx=0, start_rowx=1, end_rowx=6598)
        
        data['promoter'] = sheet.col_values(colx=27, start_rowx=1, end_rowx=6598)
        data['TSS']      = [int(round(x)) for x in sheet.col_values(colx=24, start_rowx=1, end_rowx=6598)]
        data['RBS']      = sheet.col_values(colx=28, start_rowx=1, end_rowx=6598)
        data['CDS']      = sheet.col_values(colx=76, start_rowx=1, end_rowx=6598)
        data['N_terminal_CDS'] = sheet.col_values(colx=1, start_rowx=1, end_rowx=6598)
        
        data['fivepUTR'] = [p[TSS:]+RBS+N+"CAT" for p,TSS,RBS,N in \
                                 zip(data['promoter'],data['TSS'],data['RBS'],data['N_terminal_CDS'])]
        
        num_seqs = len(data['xls_ID'])
        data['num_seqs'] = num_seqs
        data['temp']     = [30.0]*num_seqs        
        data['gene']  = ['sfGFP']*num_seqs
        data['organism'] = ['Escherichia coli str. K-12 substr. MG1655']*num_seqs        
        data['info']     = ["{} expressing {}".format(host,protein) for host,protein \
                                 in zip(data['organism'],data['gene'])]
        
        data['mRNA'] = [fivepUTR+CDS for fivepUTR,CDS in zip(data['fivepUTR'],data['CDS'])]
        data['rRNA'] = [rRNA_dic[organism] for organism in data['organism']]
        data['start_pos'] = [len(fivepUTR) for fivepUTR in data['fivepUTR']]
        
        x = np.array(sheet.col_values(colx=23, start_rowx=1, end_rowx=6598))
        data['Count.Protein'] = x.astype(float)
        
        data['Count.DNA']   = np.array(sheet.col_values(colx=3, start_rowx=1, end_rowx=6598))
        data['Count.RNA']   = np.array(sheet.col_values(colx=4, start_rowx=1, end_rowx=6598))
        data['Count.A.DNA'] = np.array(sheet.col_values(colx=5, start_rowx=1, end_rowx=6598))
        data['Count.A.RNA'] = np.array(sheet.col_values(colx=6, start_rowx=1, end_rowx=6598))
        data['Count.B.DNA'] = np.array(sheet.col_values(colx=7, start_rowx=1, end_rowx=6598))
        data['Count.B.RNA'] = np.array(sheet.col_values(colx=8, start_rowx=1, end_rowx=6598))
        
        sheet = wb.sheet_by_name("NGS counts")
        data['total.rna.a'] = sheet.cell_value(0,1)
        data['total.rna.b'] = sheet.cell_value(1,1)
        data['total.dna.a'] = sheet.cell_value(2,1)
        data['total.dna.b'] = sheet.cell_value(3,1)
        
        bin_list = ['bin'+str(num) for num in range(1,13)]
        data['bin_list'] = bin_list    

        sheet = wb.sheet_by_name("bin_vals")
        flow_seq = {b: {} for b in bin_list}
        for b,col in zip(bin_list,range(1,13)):
            flow_seq[b]['cell_fraction'] = sheet.cell_value(2,col)
            flow_seq[b]['fluo'] = sheet.cell_value(4,col)
            #left_edge = sheet.cell_value(3,col)
            #right_edge = sheet.cell_value(3,col+1)
            #flow_seq[b]['fluo'] = np.exp((np.log(left_edge) + np.log(right_edge))/2.0)            
        data['bins'] = flow_seq
        
        # read contig counts from Excel
        sheet = wb.sheet_by_index(0)
        for b,col in zip(bin_list,range(11,23)):
            x = np.array(sheet.col_values(colx=col, start_rowx=1, end_rowx=6598))
            data[b] = x.astype(float)   
            
        data = run_FS_calcs(data)
        data.pop('bins')
        data.pop('bin_list')
        
        db[paper] = {}
        db[paper] = data
        db = _format_ascii(db,paper)
        db = _subdivide_dataset(db,paper)
        
    return db

'''
def MC_sample(weights,vals)
    cum = np.cumsum(weights)
    rng = [random.random() for i in range(1000)]
    return

# CODE BLOCK THAT WAS IN process_FS_data
#=== Approximate the variance of the ratio b/w Protein and RNA ==+#
# Monte Carlo sampling of protein distribution
# breaking down from using vector math

m_prot = [] # bin values
for b in bin_list:
    m_prot.append(flow_seq[b]['fluo'])

data['sampled.Prot.dist'] = []

for i in range(len(data['num_seqs'])):
    weights = []
    for b in bin_list:
        weights.append(a[b][i])
    sampled_prot_dist = MC_sample(weights,m_prot)
    data['sampled.Prot.dist'].append(sampled_prot_dist)

# Now sample RNA distribution (two states)
weights = [0.5,0.5]
for i in range(len(data['num_seqs'])):
    vals = [data['RNA.A'],data['RNA.B']]
    sampled_RNA_dist = MC_sample(weights,vals)
    data['sampled.RNA.dist'].append(sampled_RNA_dist)
'''

def run_FS_calcs(data):
    
    # Calculate RNA levels from NGS data
    # Two replicates, A & B
    data['RNA.A'] = (data['Count.A.RNA']/data['total.rna.a']) \
                   /(data['Count.A.DNA']/data['total.dna.a'])
    
    data['RNA.B'] = (data['Count.B.RNA']/data['total.rna.a']) \
                   /(data['Count.B.DNA']/data['total.dna.b'])
    
    data['RNA'] = (data['RNA.A'] + data['RNA.B'])/2.0
    
    data['RNA.var'] = ((data['RNA.A']-data['RNA'])**2.0 + \
                       (data['RNA.B']-data['RNA'])**2.0) / 2.0
    
    data['RNA.std'] = np.sqrt(data['RNA.var'])

    # Calculate the reconstructed protein levels
    # calculate the normalized fractional contribution of each bin j per sequence i (a_ij)
    a = {}
    denominator = np.zeros(data['num_seqs'])
    for b in data['bin_list']:
        denominator += data['bins'][b]['cell_fraction']*data[b]/data['Count.Protein']
    
    # Geometric mean to calculate protein level
    data['Protein'] = np.ones(data['num_seqs'])
    for b in data['bin_list']:
        a[b] = (data['bins'][b]['cell_fraction']*data[b]/data['Count.Protein'])/denominator
        data['Protein'] *= np.exp( a[b]*np.log(data['bins'][b]['fluo']) )
    
    data['a_ij'] = a
    
    # Arithmetic mean of protein level
    data['Protein.arithmetic'] = np.zeros(data['num_seqs'])
    for b in data['bin_list']:
        data['Protein.arithmetic'] += a[b]*data['bins'][b]['fluo']
    
    # Variance & standard deviation of linear-scaled protein data
    var = 0.0
    for b in data['bin_list']:
        var += a[b]*( data['bins'][b]['fluo'] - data['Protein.arithmetic'] )**2.0        
    data['Protein.var'] = var
    data['Protein.std'] = np.sqrt(data['Protein.var'])    
    
    # Calculate apparent translation rate
    data['transl.rate'] = data['Protein']/data['RNA']
    data['transl.rate.var'] = approx_var_ratio(data['Protein'],data['Protein.var'],data['RNA'],data['RNA.var'])
    data['transl.rate.std'] = np.sqrt(data['transl.rate.var'])
    
    # Final steps
    for key,val in data.iteritems():
        if isinstance(val,np.ndarray):
            data[key] = list(val)

    data['Expression'] = data['transl.rate']
    
    return data

def approx_var_ratio(mu1,var1,mu2,var2):
    # assumed cov(x,y) = 0
    var = (mu2**2.0 * var1**2.0 + mu1**2.0 * var2**2.0) / mu2**2.0
    return var

def _get_filepath(ds):
    return "spreadsheets/{}.xls".format(ds)

def _format_ascii(db,paper):
    for key,val in db[paper].iteritems():
        if isinstance(val,list) and isinstance(val[0],basestring):
            db[paper][key] = map(lambda x: x.encode('ascii'), val)
    return db

def _subdivide_dataset(db,paper):
    
    set_labels = list(set(db[paper]['info']))
    unique_sets = {l: [] for l in set_labels}
    for indx,val in zip(range(db[paper]['num_seqs']),db[paper]['info']):
        unique_sets[val].append(indx)

    temp = {}
    for key,val in db[paper].iteritems():
        if isinstance(val,list):
            subdvd_list = []
            for l in set_labels:
                sub = []
                for i in unique_sets[l]:
                    sub.append(val[i])
                subdvd_list.append(sub)
            temp[key] = subdvd_list
        else:
            temp[key] = val
    
    temp['subset_size'] = []
    for subset in temp['mRNA']:
        temp['subset_size'].append(len(subset))
    
    db[paper] = temp
    db[paper]['num_subsets'] = len(set_labels)
    
    return db

if __name__ == "__main__":
    
    datasets = ['EspahBorujeni_NAR_2013',   # StandbySite
                'EspahBorujeni_NAR_2015',   # Riboswitch
                'EspahBorujeni_JACS_2016',  # RibosomeDrafting
                'EspahBorujeni_Footprint',  # Footprint
                'Salis_Nat_Biotech_2009',   # RBSCalculator
                'Farasat_MSB_2014',         # dRBS
                'Tian_NAR_2015',            # TranslationCoupling
                'Mimee_Cell_Sys_2015',      # Bthetaiotaomicron
                'Bonde_NatMethods_IC_2016', # EMOPEC
                'Kosuri_PNAS_2013',         # Kosuri
                'Goodman_Science_2013'      # Goodman
                ]
    
    mRNA_database = {}
    mRNA_database = _import(mRNA_database,datasets)
    
    import cPickle
    handle = open('mRNA_database.p','w')
    cPickle.dump(mRNA_database, handle)
    handle.close()
    
