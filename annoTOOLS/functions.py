"""
This script creates a new database from a gff3 file

It will produce two outputs:

    + a new gff3 with a coherent nesting structure.
    + a crossreference of old <-> new ids for tractable data.
"""

from annoTOOLS.classes import gff3_entry
from annoTOOLS.classes import gene

gff3 = "../sandbox/SMAN.gff3" # For testing purposes

def read_gff3( gff3 ):
    """
    Convert gff file into gff3 dictionary with ID as dict entry.
    """

    gff3_dict = {}
    with open( gff3, 'r' ) as handle:
        for line in handle:
            if not line.startswith("#"):
                
                entry = gff3_entry( line )

                if entry.unsupported:
                    pass
                else:
                    # Both id and feature must be specified b/c some CDS and mRNA share ids
                    gff3_dict[ (entry.id, entry.feature) ] = entry

    return gff3_dict


def update_nesting( gff3_file, gff3_dictionary ):

    """
    Update gene nesting structure for exon and CDS.
    """

    with open( gff3_file , 'r' ) as handle:
        for line in handle:
            if not line.startswith("#"):
                entry = gff3_entry(line)
                # If feature is not supported, skip
                if entry.unsupported:
                    pass
                else:
                    # Check if its exon or CDS
                    if entry.feature == "exon" or entry.feature == "CDS":
                        # Try if mRNA is present, if its not, do nothing
                        try:
                            gene = gff3_dictionary[ (entry.mRNA, "mRNA") ].gene
                            gff3_dictionary[entry.id, entry.feature ].add_gene_feature( entry.feature, gene )
                        except KeyError:
                            pass


                    else:
                        pass


def gff3_to_gene( gff3_dict , species ):
    
    """
    Turn a gff3_entry dictionary into a gene dictionary

    number is a counter number for id generation
    """

    gene_dict = {}

    # Do first traversal to add genes
    counter = 0
    for key in gff3_dict:
        
        entry = gff3_dict[ key ]
        if entry.feature == "gene":

            g=gene( species,
                    counter,
                    entry.chromosome,
                    entry.start,
                    entry.end,
                    entry.strand,
                    entry.id )
            gene_dict[ entry.id ] = g
            counter+=1

            
    # Do second traversal to add nested features
    for key in gff3_dict:

        entry = gff3_dict[key]

        if entry.feature == "gene":
            pass
        elif entry.feature == "exon":
            gene_dict[ entry.gene ].add_exon( entry.start, entry.end, entry.id )
        elif entry.feature == "mRNA":
            gene_dict[ entry.gene ].add_mrna( entry.start, entry.end, entry.id )
        elif entry.feature == "CDS":
            gene_dict[ entry.gene ].add_cds( entry.start, entry.end, entry.id  )
        else:
            pass
    

        
    return gene_dict    



















