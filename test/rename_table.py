#!/usr/bin/env python3

import os
import sys

sys.path.insert( 0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from annoTOOLS.classes import *
from annoTOOLS.functions import *


####
#sp_gff = sys.argv[1]
#sp_gff = sys.argv[2]

####
sp_gff = "../../sandbox/SMAN.gff3"
sp_name = "SMAN"
##


# First, read gff3 entry into variable
sp = read_gff3( sp_gff )


# Now, update nesting structure (gene,protein,mRNA,exon)

update_nesting( sp_gff, sp )

# Convert to gene class

gene_dict = gff3_to_gene( sp, sp_name )



# Print table of oldname -> newname


for i in gene_dict:

    gene_dict[ i ].print_rename_table()


