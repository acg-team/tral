# (C) 2013 Elke Schaper

from os.path import join
import re
import logging
import numpy as np
from collections import defaultdict

#from . import repeat_info

logger = logging.getLogger(__name__)

'''
Functions to check whether the exon structure and the repeat unit structure are correlated.

'''


################################### Repeat exon structure #########################################

def get_exon_structure(my_TR, exon_structure):

    '''
    Exon structure parameters:
    # exon_structure = [{'Exon Chr End (bp)': '43478424', 'Exon Chr Start (bp)': '43477252', 'CDS End': '272', 'CDS Start': '1'}, {'Exon Chr End (bp)': '43476658', 'Exon Chr Start (bp)': '43476498', 'CDS End': '433', 'CDS Start': '273'}, {'Exon Chr End (bp)': '43476158', 'Exon Chr Start (bp)': '43476036', 'CDS End': '556', 'CDS Start': '434'}, {'Exon Chr End (bp)': '43475664', 'Exon Chr Start (bp)': '43475564', 'CDS End': '657', 'CDS Start': '557'}, {'Exon Chr End (bp)': '43475416', 'Exon Chr Start (bp)': '43475194', 'CDS End': '880', 'CDS Start': '658'}, {'Exon Chr End (bp)': '43475046', 'Exon Chr Start (bp)': '43474707', 'CDS End': '951', 'CDS Start': '881'}]
    '''

    # Extract the exon structure in protein coordinates
    exon_start_aa = [int(i['CDS Start']) for i in exon_structure if 'CDS Start' in i and i['CDS Start'] != '']
    exon_end_aa = [int(i['CDS End']) for i in exon_structure if 'CDS End' in i and i['CDS End'] != '']
    
    if not exon_end_aa == []:
        protein_length = exon_end_aa[-1]
    else:
        protein_length = None
        
    # Extract the TR structure in protein coordinates
    index = my_TR.begin
    tr_start_aa = []
    for iMSA in my_TR.msa:
        unit_length = len(iMSA.replace("-",""))
        tr_start_aa.append(index)
        index += unit_length
        
    return exon_start_aa, protein_length, tr_start_aa, index


def get_exon_measures(my_TR, exon_structure):
    
    '''
    Calculate diverse exon measures
    
    Exon structure parameters:
    # exon_structure = [{'Exon Chr End (bp)': '43478424', 'Exon Chr Start (bp)': '43477252', 'CDS End': '272', 'CDS Start': '1'}, {'Exon Chr End (bp)': '43476658', 'Exon Chr Start (bp)': '43476498', 'CDS End': '433', 'CDS Start': '273'}, {'Exon Chr End (bp)': '43476158', 'Exon Chr Start (bp)': '43476036', 'CDS End': '556', 'CDS Start': '434'}, {'Exon Chr End (bp)': '43475664', 'Exon Chr Start (bp)': '43475564', 'CDS End': '657', 'CDS Start': '557'}, {'Exon Chr End (bp)': '43475416', 'Exon Chr Start (bp)': '43475194', 'CDS End': '880', 'CDS Start': '658'}, {'Exon Chr End (bp)': '43475046', 'Exon Chr Start (bp)': '43474707', 'CDS End': '951', 'CDS Start': '881'}]
    '''
    
    exon_start_aa, protein_length, tr_start_aa, exon_stop_aa = get_exon_structure(my_TR, exon_structure)
    if exon_start_aa == [] or tr_start_aa == []:
        return None
    
    nExon = len([i for i in exon_start_aa if i > tr_start_aa[0] and i < tr_start_aa[-1]]) + 1
    
    nTR_in_Exon_Max = 0
    for iExon_start, iExon_end in zip(exon_start_aa[:-1],exon_start_aa[1:]):
        lTR_Start_in_Exon = [iTR_start for iTR_start in tr_start_aa if iTR_start >= iExon_start and iTR_start < iExon_end]
        if len(lTR_Start_in_Exon) > nTR_in_Exon_Max:
            nTR_in_Exon_Max = len(lTR_Start_in_Exon)
            
    return nExon, nTR_in_Exon_Max


def print_exon_structure(my_TR, exon_structure):
    
    '''
    Exon structure parameters:
    # exon_structure = [{'Exon Chr End (bp)': '43478424', 'Exon Chr Start (bp)': '43477252', 'CDS End': '272', 'CDS Start': '1'}, {'Exon Chr End (bp)': '43476658', 'Exon Chr Start (bp)': '43476498', 'CDS End': '433', 'CDS Start': '273'}, {'Exon Chr End (bp)': '43476158', 'Exon Chr Start (bp)': '43476036', 'CDS End': '556', 'CDS Start': '434'}, {'Exon Chr End (bp)': '43475664', 'Exon Chr Start (bp)': '43475564', 'CDS End': '657', 'CDS Start': '557'}, {'Exon Chr End (bp)': '43475416', 'Exon Chr Start (bp)': '43475194', 'CDS End': '880', 'CDS Start': '658'}, {'Exon Chr End (bp)': '43475046', 'Exon Chr Start (bp)': '43474707', 'CDS End': '951', 'CDS Start': '881'}]
    '''
    
    exon_start_aa, tr_start_aa, protein_length = get_exon_structure(my_TR, exon_structure)
    
    if exon_start_aa == [] or tr_start_aa == []:
        return ""
        
    protein_structure = ""   
    for i in range(1,protein_length+1):
        if i in exon_start_aa:
            if i in tr_start_aa:
                protein_structure += "B"
            else:
                protein_structure += "E"
        elif i in tr_start_aa:
            protein_structure += "x"
        else:
            protein_structure += "_"
            
    return protein_structure
        
    


def main():

    exon_structure = [{'Exon Chr End (bp)': '43478424', 'Exon Chr Start (bp)': '43477252', 'CDS End': '272', 'CDS Start': '1'}, {'Exon Chr End (bp)': '43476658', 'Exon Chr Start (bp)': '43476498', 'CDS End': '433', 'CDS Start': '273'}, {'Exon Chr End (bp)': '43476158', 'Exon Chr Start (bp)': '43476036', 'CDS End': '556', 'CDS Start': '434'}, {'Exon Chr End (bp)': '43475664', 'Exon Chr Start (bp)': '43475564', 'CDS End': '657', 'CDS Start': '557'}, {'Exon Chr End (bp)': '43475416', 'Exon Chr Start (bp)': '43475194', 'CDS End': '880', 'CDS Start': '658'}, {'Exon Chr End (bp)': '43475046', 'Exon Chr Start (bp)': '43474707', 'CDS End': '951', 'CDS Start': '881'}]

    #[{'Exon Chr End (bp)': '120068186', 'Exon Chr Start (bp)': '120067591', 'CDS End': '500', 'CDS Start': '1'}, {'Exon Chr End (bp)': '120054800', 'Exon Chr Start (bp)': '120054672', 'CDS End': '629', 'CDS Start': '501'}, {'Exon Chr End (bp)': '120053986', 'Exon Chr Start (bp)': '120053709', 'CDS End': '907', 'CDS Start': '630'}, {'Exon Chr End (bp)': '120050255', 'Exon Chr Start (bp)': '120043356', 'CDS End': '1116', 'CDS Start': '908'}]
    exon_structure = my_gene_tree['leafs']['ENSP00000295628']['exon']
    my_TR = my_TR_data['homologs']['ENSP00000295628']['TR']['HMM']
    
    exon_structure_aa = [[i['CDS Start'],i['CDS End']] for i in exon_structure]



    


