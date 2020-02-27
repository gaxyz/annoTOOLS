#!/usr/bin/env python3

import os
import sys

sys.path.insert( 0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from annoTOOLS.classes import *
from annoTOOLS.functions import *

root_dir = "/export/home/grijo/projects/platy/data/species_raw"
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
outfile = "/export/home/grijo/projects/platy/data/crossref.tab"

with open( outfile , 'w' ) as handle:

    for species in sp_list:

       sp_name = species
       sp_gff = root_dir + "/" + species +  "/" + species + ".gff3"
       print( "Processing {0} annotation from {1} file...".format(sp_name,sp_gff ) ) 
       # First, read gff3 entry into variable
       print(">>>>> Reading gff")
       sp = read_gff3( sp_gff )


       # Now, update nesting structure (gene,protein,mRNA,exon)
       print(">>>>> Updating nesting structure")
       update_nesting( sp_gff, sp )

       # Convert to gene class
       print(">>>>> Converting to gene class")
       gene_dict = gff3_to_gene( sp, sp_name )

       # Print table of oldname -> newname

       print(">>>>> Writing crossreference table")
       for i in gene_dict:

           handle.write( gene_dict[ i ].get_rename_table() )

       print("Done! <<<<<<<<<<<<<<<<<<") 
