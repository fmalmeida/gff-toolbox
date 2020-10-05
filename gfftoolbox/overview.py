## Def overview help
usage_overview="""
gff-toolbox:

            Overview

This command uses BCBio library to let you get the gist of your GFF file

usage:
    gff-toolbox overview [ -h|--help ] [ --input <gff> ]

options:
    -h, --help                                               Show this screen
    -i, --input=<gff>                                        Input GFF file. GFF file must not contain sequences with it. [Default: stdin]

example:

    ## Getting the overview of a generic GFF file

$ gff-toolbox overview -i Athaliana_ref.gff.gz

$ gff-toolbox overview -i Kp_ref.gff
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
import binascii
import gzip

###################
### Stdin Check ###
###################
def stdin_checker(input):
    # Checking for stdin
    if input == "stdin":
        tmp = tempfile.NamedTemporaryFile(mode = "w+t", delete = False) # Create tmp file to work as input
        temp_file = open(tmp.name, 'w')
        for line in sys.stdin:
            temp_file.writelines(f"{line}")

        temp_file.seek(0)
        input = tmp.name

    else:
        pass

    return input

###################
### Gzip opener ###
###################
def gzip_opener(input, mode_in):

    if  binascii.hexlify(open(input, 'rb').read(2)) == b"1f8b" or str(input).endswith(".gz"):
        return gzip.open(input, mode=mode_in)
    else:
        return open(input, mode=mode_in)

##################################################
### Function for checking available qualifiers ###
##################################################
def overview(infile):

    ## Module
    examiner = GFFExaminer()

    ## Open connection
    summary = examiner.available_limits(gzip_opener(stdin_checker(infile), "rt"))

    ## Print
    print(f"""
Note: This command lets you get the gist of your gff file. It summarises as a dictionary the available
limits of your gff, which means, that it counts the ids, sources and feature types that are found in the
input gff. Also, we check the 9th column to get an overview of the available attributes/qualifiers.

Input: {infile} gff file.

Overview:
    """)
    pprint(summary)
