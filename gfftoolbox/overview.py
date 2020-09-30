## Def overview help
usage_overview="""
gff-toolbox:

            Overview

This command uses BCBio library to let you get the gist of your GFF file

usage:
    gff-toolbox overview [ -h|--help ] [ --input <gff> ]

options:
    -h, --help                                               Show this screen
    -i, --input=<gff>                                       Input GFF file. GFF file must not contain sequences with it. [Default: stdin]

example:

    ## Getting the overview of a generic GFF file
    
$ gff-toolbox overview -i input.gff
"""

##################################
### Loading Necessary Packages ###
##################################
import Bio.SeqIO
from pprintpp import pprint
from BCBio import GFF
from BCBio.GFF import GFFExaminer
import tempfile
import sys

##################################################
### Function for checking available qualifiers ###
##################################################
def check_gff(infile):

    # Checking for stdin
    if infile == "stdin":
        tmp = tempfile.NamedTemporaryFile(mode = "w+t") # Create tmp file to work as input
        temp_file = open(tmp.name, 'w')
        for line in sys.stdin:
            temp_file.writelines(f"{line}")

        temp_file.seek(0)
        infile = tmp.name

    ## Module
    examiner = GFFExaminer()

    ## Open connection
    in_handle = open(infile)
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
#     print("""
# Attributes:
#
# Note: Showing the available qualifiers and attributes found in the the first record of the GFF. These
# qualifiers can be used in the other commands of the toolbox for filtering, plotting, etc.
#     """)
#
#     # Check qualifiers
#     for rec in GFF.parse(infile):
#         print(rec.features[0])
#         exit()
