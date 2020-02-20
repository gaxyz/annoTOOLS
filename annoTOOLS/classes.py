#!/usr/bin/env python3



##########################
#    CLASSES             #
##########################

class gene:
    """Class made for storing exons, mRNAs and CDSs within genes"""

    def __init__(self, species, number, chromosome, start, end, strand, original_id):
    
        # Define gene attributes
        self.species = str(species)
        self.number = str(number)
        self.chromosome = str(chromosome)
        self.start = int(start)
        self.end = int(end)
        self.strand = str(strand)
        self.original_id = str(original_id)

        # Initialize attributes for exons

        self.exons = dict()

        # Initialize attributes for mRNAs

        self.mRNAs = dict()

        # Initialize attributes for CDSs    
    
        self.CDSs = dict()

        self.newnumber = ("0"*7)[0:len(self.number) - 1] + self.number 

        self.name = self.species + "GEN" + self.newnumber

    
    def add_exon(self, start, end , oldname):

        # Assign number
        #### NUMBER ASSIGNMENT NOT IN GENOMIC ORDER!!!!!!
        
        n = len(self.exons) + 1
        
        # Name it
        
        name = self.species + "EXO" + self.newnumber + "." + str(n)

        # Convert to within-gene coordinate

        start = int(start) - self.start
        end = int(end) - self.start

        # Add to dict

        self.exons[name] = [int(start), int(end), oldname ]

    def add_mrna(self, start, end , oldname ):

        # Assign number
        n = len(self.mRNAs) + 1 
        # Name it

        name = self.species + "MRN" + self.newnumber + "." + str(n)

        # Convert to within-gene coordinate

        start = int(start) - self.start
        end = int(end) - self.start

        # Add to dict

        self.mRNAs[name] = [int(start), int(end), oldname ]

    def add_cds(self, start, end, oldname ):

        # Assign number
        n = len(self.CDSs) + 1
        # Name it

        name = self.species + "CDS" + self.newnumber + "." + str(n)

        # Convert to within-gene coordinate

        start = int(start) - self.start
        end = int(end) - self.start

        # Add to dict


        self.CDSs[name] = [int(start) , int(end), oldname ]

    # Method to get longest (specified) feature ID

    def get_longest(self, feature):

        dictionary = {"CDS": self.CDSs ,"MRNA": self.MRNAs , "EXON": self.exons}
        longest_len = 0
        longest_id = ""
        for key in dictionary[feature]:

            length = dictionary[feature][key][1] - dictionary[feature][key][0]
            if length > longest_len:
                longest_id = key

        return longest_id
            
    # Method to print OLDNAME -> NEWNAME table

    def print_rename_table( self ):

        # Format:
        # OLDNAME FATURE NEWNAME

        # gene
        g = "{0} gene {1}".format( self.original_id, self.name )
        print( g )
        # exons
        for entry in self.exons:
            e = "{0} exon {1}".format( self.exons[entry][2], entry )
            print(e)

        # mrnas
        for entry in self.mRNAs:
            e = "{0} mRNA {1}".format( self.mRNAs[entry][2], entry )
            print(e)

        # cdss
        for entry in self.CDSs:
            e = "{0} CDS {1}".format( self.CDSs[entry][2], entry )
            print(e)




class gff3_entry:
    """Class for containing gff information. Only CDS, mRNA, exon and gene entries are kept"""

    def __init__(self, line):

        sline = line.split()
        self.chromosome = sline[0]
        self.source = sline[1]
        self.feature = sline[2]
        self.start = sline[3]
        self.end = sline[4]
        
        self.score = sline[5]
        self.strand = sline[6]
        self.phase = sline[7]
       

        # Define accepted features
        self.unsupported = False 
        accepted_features = ["gene", "mRNA", "CDS", "exon"]
        if self.feature not in accepted_features:
            self.unsupported = True
         
            
        if not self.unsupported:
        
            attributes = sline[8].split(";")
            # Make attributes explicit

            att_dict = {}
            for i in attributes:
     
                attribute = i.split("=")[0]
                # We lose some information here
                if attribute == "ID" or attribute == "Parent":
                    att_dict[attribute] = i.split("=")[1].split(":")[1]

            # Handle attributes to allow for proper nesting structure ( genes(mRNA(CDS, exon) )
            # Second level(CDS, exon) has to be initialized with none and then added after reading gff
            self.id = att_dict["ID"]

            if self.feature == "exon":
                # Initinalize
                self.gene = None
                self.exon = self.id
                self.mRNA = att_dict["Parent"]
                self.CDS = None

            elif self.feature == "mRNA":
                # MRNA
                self.gene = att_dict["Parent"]
                self.exon = None
                self.mRNA = self.id
                self.CDS = None

            elif self.feature == "CDS":
                self.gene = None
                self.exon = None
                self.mRNA = att_dict["Parent"]
                self.CDS = self.id

            elif self.feature == "gene":
                self.gene = self.id
                self.exon = None
                self.mRNA = None
                self.CDS = None

    def get_entry(self):

        """
        Method for returning information in gff3 format.
        General annotation except IDs is not returned (it is lost in reading the entry)
        """
        # Collapse attributes
        # Some information has been lost when coercing into gff3_entry class
        attributes = "ID={0};gene={1};mRNA={2};exon={3};CDS={4}".format(self.id,
                                                                        self.gene,
                                                                        self.mRNA,
                                                                        self.exon,
                                                                        self.CDS)
        # Gather info
        entry = "{0} {1} {2} {3} {4} {5} {6} {7} {8}".format(self.chromosome,
                                                     self.source,
                                                     self.feature,
                                                     self.start,
                                                     self.end,
                                                     self.score,
                                                     self.strand,
                                                     self.phase,
                                                     attributes)
        # Return entry
        return entry



    def add_gene_feature( self, feature , gene_id):

        """
        Add a gene feature to gff3 entry.
        
        This is intended for exons and CDS as
        they lack an explicit link to their 
        gene.

        Parameters:

        feature  For safety purposes (dont wanna modify genes or mRNAs). 
                 It is either exon or CDS.

                 
        gene_id  The gene id to be associated with either the exin or the CDS.
        """


        # This function is intended for modifying exons and CDS as they lack a gene-link
        if feature == "exon" or feature == "CDS":
            self.gene = gene_id

        else:
            print( "Error: Only exon and CDS entries are modifiable with add_gene_feature() method." )

