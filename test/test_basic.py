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

# Print to test

for i in sman:
    print(sman[i].get_entry())
