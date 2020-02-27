#!/usr/bin/env python3 
import  sys

ingff3 = sys.argv[1]
outgff3 = sys.argv[2]

with open( ingff3 , 'r' ) as infile:
    with open( outgff3 , 'w' ) as outfile:

        for line in infile:

            
            if not line.startswith("#"):
                
                feature = line.split()[2]


                if feature == "exon":
                    linelist = line.split()
                    attributes = linelist[8]

                    attributes = attributes.replace( "exon:Transcript", "exon" )
                    attributes = attributes.replace( ":" ,"." , 2)
                    attributes = attributes.replace(".", ":" , 1 )

                    linelist[8] = attributes
                    newline = "\t".join(linelist) + "\n"

                    outfile.write(newline)

                elif feature == "CDS":

                    linelist = line.split()
                    attributes = linelist[8]

                    attributes = attributes.replace( "CDS:Transcript", "CDS" )
                    attributes = attributes.replace( ":" ,"." , 2)
                    attributes = attributes.replace(".", ":" , 1 )

                    linelist[8] = attributes
                    newline = "\t".join(linelist) + "\n"

                    outfile.write(newline )


                else:
                    newline = line
                    outfile.write(newline)
                

            else:
                outfile.write(line)


