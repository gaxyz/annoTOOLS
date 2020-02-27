#!/usr/bin/env python3
"""
This script builds proteomes from a species list

The proteomes are composed of longest CDS for proper orthology inference
"""

import os
import sys
import pickle
from Bio import SeqIO


print( "> Build longest CDS proteomes" )

sys.path.insert( 0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from annoTOOLS.classes import *
from annoTOOLS.functions import *


sp_list = ["SMAN",
            "SJAP",
            "TREG",
            "FHEP",
            "FGIG",
            "FBUS",
            "CSIN",
            "OVIV",
            "SMED"
            ]
pickle_dir = "/export/home/grijo/projects/platy/data/pickle/geneDicts"
protein_root_dir = "/export/home/grijo/projects/platy/data/species_raw"
output_dir = "/export/home/grijo/projects/platy/data/proteins"

print( "\n>> Looping through species..." )

for sp in sp_list:
    
    print( "\n>>> READING {0} SPECIES... <<<".format(sp) )
    faa_file = protein_root_dir +  "/" + sp + "/" + sp + ".faa"
    gene_pickle = pickle_dir + "/" + sp + ".genedict.pkl"

    # Read proteome into dict
    print( ">>>>>> Loading proteome..." )
    proteome = {}
    for seq in SeqIO.parse( faa_file, "fasta" ):
        proteome[seq.id] = seq.seq

    # Read pickle
    print( ">>>>>> Loading annotation (pickle)..." )
    pickle_in = open(gene_pickle, 'rb')
    gene_dict = pickle.load(pickle_in)
    
    
    # Get longest CDSs   
    print( ">>>>>> Getting longest CDS..." )
    namedict = {}
    for gene in gene_dict:
        names = gene_dict[gene].get_longest( "CDS" )
        if names is not None:
            
            newname = names[0]
            oldname = names[1]
            # FGIG and FBUS are not annotated correctly
            if sp in ["FBUS", "FGIG"]:
                newname = newname.split(".")[0]
                oldname = oldname.split(".")[0]

            namedict[ oldname ]  = newname


    print( ">>>>>> Writing longest CDS proteome..." )
    with open( output_dir + "/" + sp + ".faa", 'w' ) as handle:
        # Write only longest proteins and change names
        for entry in proteome:
            try:
                header = ">" + namedict[ entry ] + "\n"
                sequence  = str(proteome[entry]) + "\n"
                handle.write( header )
                handle.write( sequence )
            except KeyError:
                pass
            
               
                





    print( ">>> SPECIES {0} DONE! <<<".format(sp))
    pickle_in.close()
    

