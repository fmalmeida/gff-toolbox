## Def overview help
usage_overview="""
gff-toolbox:

            Overview

This command uses BCBio library to let you get the gist of your GFF file

usage:
    gff-toolbox overview [ -h|--help ] [ --input <gff> ]

options:
    -h --help                                               Show this screen
    -i, --input=<gff>                                       Input GFF file. GFF file must not contain sequences with it

example:
    gff-toolbox overview -i input.gff
"""

##################################
### Loading Necessary Packages ###
##################################
import Bio.SeqIO
from pprintpp import pprint
from BCBio import GFF
from BCBio.GFF import GFFExaminer

##################################################
### Function for checking available qualifiers ###
##################################################
def check_gff(infile):

    ## File
    in_file = infile

    ## Module
    examiner = GFFExaminer()

    ## Open connection
    gff = GFF.parse(infile)
    in_handle = open(in_file)
    summary = examiner.available_limits(in_handle)
    in_handle.close()

    ## Print
    print(f"""
Note: This command lets you get the gist of your gff file. It summarises as a dictionary the available
limits of your gff, which means, that it counts the ids, sources and feature types that are found in the
input gff. Also, we check the 9th column to get an overview of the available attributes/qualifiers.

Input: {infile} gff file.

Overview:
    """)
    pprint(summary)
    print("""
Attributes:

Note: Showing the available qualifiers and attributes found in the the first record of the GFF. These
qualifiers can be used in the other commands of the toolbox for filtering, plotting, etc. The qualifiers
are the fields separated by ';' found in the GFF attributes column (9th column).
    """)

    # Check qualifiers
    for rec in gff:
        print(rec.features[0])
        exit()
