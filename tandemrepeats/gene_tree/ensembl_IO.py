# Copyright (C) 2012 by Elke Schaper (elke@inf.ethz.ch) 

"""
    IO functions for Ensembl Files 
    Current emphasis: Ensembl Xref Tables
""" 

import logging
import re
import glob
import os

import logging
logging.basicConfig(level=logging.DEBUG)

def clean_sequence_file(infile_name):
    """ remove the following type of lines from Ensembl Biomart files:
>
Sequence unavailable
    """
    
    with open(infile_name, 'r') as f:
        lines = f.readlines()

    # patterns
    pat_header = '>'
    pat_empty = '>$'
    
    # States
    # 1: searching empty >
    # 2: anything after empty >

    state = 1
    with open(infile_name, 'w') as f:
        for i, line in enumerate(lines):
    
            if 2 == state:
                if re.match(pat_header, line):
                    state = 1
                    logging.debug(" * (2->1) Found first header after empty header")
            
            if 1 == state:
                if re.match(pat_empty, line):
                    state = 2
                    logging.debug(" * (1->2) Found empty header")
                else:
                    f.write(line)
                
        

    



def get_data(infile, dict_ID=None, all=False):
    """Read Data from Ensembl Table for <dict_ID>: (ID_tag, ID)
     At current load file each time you need an entry. Alternative: Read whole file and save as dict
    
    Arguments:
    infile: File stream
    
    Returns a generator function yielding one gene tree and all its sequences per request.
    
    Postcondition: infile points to EOF.
    """
    
    ## patterns
    # No patterns are necessary, as we are dealing with a simple tab separated file
        
    # Our possible parser states:
    #
    # 1: searching for header
    # 2: searching for appropiate data
    
    state = 1
    for i, line in enumerate(infile):
        #logging.debug("Line: {0}: {1}".format(str(i), line[0:-1]))
        if 1 == state:
            header = line.split("\t")
            if len(header) > 0:
                state = 2
                result = []
                header = [iHeader.strip() for iHeader in header]
                logging.debug(" * (1->2) Found header: {0}".format("*".join(header)))
                try:
                    index_ID = header.index(dict_ID[0])
                    logging.debug("index_ID: {0}".format(index_ID))
                except:
                    logging.debug("Column name not in header: {0}".format(dict_ID[0]))
                    return None
                continue
        
        if 2 == state:
            data = line.rstrip().split("\t")
            if len(data) > index_ID:
                if data[index_ID] == dict_ID[1]:
                    result.append({iHeader: iData for iHeader,iData in zip(header, data)})
                    if all == False:
                        return result[0]
    else:
        if result == []:
            logging.debug(" No data was found for: {0}".format(dict_ID[1]))
            return None
        else:
            return result
                    
                
if __name__ == "__main__":

    clean_sequence_file("/Users/elkeschaper/Downloads/petromyzon_marinus_unspliced_gene.fasta")
    clean_sequence_file("/Users/elkeschaper/Downloads/pongo_abelii_unspliced_gene.fasta")
    clean_sequence_file("/Users/elkeschaper/Downloads/pteropus_vampyrus_unspliced_gene.fasta")
    clean_sequence_file("/Users/elkeschaper/Downloads/petromyzon_marinus_cds.fasta")
    #os.chdir("/Users/elkeschaper/Downloads/")
    #for iFile in glob.glob("*.fasta"):
    #    print(iFile)
    #    clean_sequence_file(iFile)