## Def filter help
usage_filter="""
gff-toolbox:

            Filter

This command uses Biopython library to parse and filter your GFF file as you wish. It targets the
attributes column of the GFF.

usage:
    gff-toolbox filter [-h|--help ]
    gff-toolbox filter [ --mode loose ] [ --input <gff> ] [ --pattern <string> --column <int> --sort --header ]
    gff-toolbox filter [ --mode exact ] [ --input <gff> ] [ --chr <chr_limits> ]

options:

                                                Generic parameters

    -h --help                                               Show this screen.

    -i, --input=<gff>                                       Input GFF file. GFF file must not contain sequences, only features.

    -m, --mode=<search_mode>                                In which mode to search for patterns: loose or exact? In loose mode,
                                                            the GFF is scanned in a grep-like manner via pandas dataframes in which
                                                            user must specify a pattern and a column to do the search. This mode is
                                                            recommended for simple searches with simple GFFs (not nested). The exact
                                                            mode scans the GFF with Biopython and BCBio packages, treating it as
                                                            python dictionary. It is recommended for more complex searches and complex
                                                            GFFs, such as nested GFFs. [Default: exact]

                                                Loose search mode parameters

    -c, --column=<int>                                      Apply pattern search in which GFF columns? [Default: 9]

    -p, --pattern=<string>                                  Pattern to search in the GFF file. Can be a list of patterns separated by commas.

    --sort                                                  Sort the GFF by the contig and start position.

    --header                                                Print GFF header (##gff-version 3)? Some programs require this header.

                                                Exact search mode parameters

    --chr=<chr_limits>                                      Apply a filter based on the chr/contig/sequence ids (Column 1). Can be a list of patterns separated by commas.


example:

    ## Simple filter in any column: wheter a line contain a pattern in a specific column (like grep)

$ gff-toolbox filter --mode loose --sort --header -i test/input.gff -c 2 -p "barrnap:0.9"

    ## Explain

$ gff-toolbox filter -i test/input.gff --chr contig_7_segment0
"""

##################################
### Loading Necessary Packages ###
##################################
import sys
import pandas as pd
import re
from BCBio import GFF

####################################
### Function to import gff as df ###
####################################
def read_gff_df(input):
    return pd.read_csv(input, sep = "\t", comment = "#",
                        names=['1', '2', '3', '4', '5', '6', '7', '8', '9'])

# Guide: ['Chr', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'

######################################
### Function to import gff as dict ###
######################################
def read_gff_dict(input, chr_limits):

    # Check for the limits imposed by the user
    limit_info = {}

    # Chr limit?
    if chr_limits != None:
        chr_list = list(chr_limits.split(','))
        limit_info['gff_id'] = chr_list

    in_handle = open(input)
    return GFF.parse(in_handle, limit_info=limit_info)
    in_handle.close()


###################################
### Filter df by column pattern ###
###################################
def df_col_filter(df, column, pattern):

    # Split csv
    pat_list = list(pattern.split(','))

    # Filter
    return df[
        df[str(column)].str.contains('|'.join(pat_list))
    ].to_csv(sep='\t', index=False, header=False).strip()

######################################################
### Function for simple filter with single pattern ###
######################################################
def filter_loose_mode(input_gff, column, pattern, sort, header):

    # Read GFF file
    df = read_gff_df(input_gff)

    # Sort
    if sort == True:
        df = df.sort_values(by=['1', '4'])

    # header
    if header == True:
        print("##gff-version 3")

    # Filter
    print(df_col_filter(df=df, column=str(column), pattern=str(pattern)))

#######################################################
### Function for complex filter with single pattern ###
#######################################################
def filter_exact_mode(input_gff, chr_limits):

    # Parse fields
    gff_dict = read_gff_dict(input=input_gff, chr_limits=chr_limits)

    # Check
    for rec in gff_dict:
        GFF.write([rec], sys.stdout)

################
### Def main ###
################
def filter(input_gff, column, pattern, sort, header, mode, chr_limits):

    # Simple filter
    if mode == "loose":
        filter_loose_mode(input_gff=input_gff, column=column, pattern=pattern, sort=sort, header=header)

    # Complex filter
    elif mode == "exact":
        filter_exact_mode(input_gff=input_gff, chr_limits=chr_limits)

    # Error
    else:
        print(f"""
Error: --mode must be either 'loose' or 'exact'. {mode} is incorrect.
        """)
