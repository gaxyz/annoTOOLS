#!/usr/bin/env python3

import os
import sys

sys.path.insert( 0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from annoTOOLS.classes import *
from annoTOOLS.functions import *


# First, read gff3 entry into variable

sman_gff = "../../sandbox/SMAN.gff3"
sman = read_gff3( sman_gff )


# Now, update nesting structure (gene,protein,mRNA,exon)

update_nesting( sman_gff, sman )

# Convert to gene class

gene_dict = gff3_to_gene( sman, "SMAN" )


# Print table of oldname -> newname


for i in gene_dict:

    gene_dict[ i ].print_rename_table()


