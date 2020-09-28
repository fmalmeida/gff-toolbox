## Def filter help
usage_filter="""
gff-toolbox:

            Filter

This command uses Biopython library to parse and filter your GFF file as you wish. It targets the
attributes column of the GFF.

usage:
    gff-toolbox filter [-h|--help ] [ --input <gff> --mode <search_mode> ] [ --pattern <string> ] [ --column <int> --sort --header ]

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

    --sort                                                  Sort the GFF by the contig and start position.

    --header                                                Print GFF header (##gff-version 3)? Some programs require this header.

                                                Loose search mode parameters

    -c, --column=<int>                                      Apply pattern search in which GFF columns? [Default: 9]

    -p, --pattern=<string>                                  Pattern to search in the GFF file. Can be a list of patterns separated by commas.

                                                Exact search mode parameters


example:

    ## Simple filter in any column: wheter a line contain a pattern in a specific column (like grep)

$ gff-toolbox filter -i test/input.gff -c 2 -p "barrnap:0.9" --sort

    ## Complex filter, based on the fields found in a given column. In the example we search for lines
    ## where the product field in the attributes column (9th col) are "16S ribosomal RNA".
    ## In the example GFF, attributes are separated by ";" and its values are marked by "="

$ gff-toolbox filter -i test/input.gff -c 9 -p "16S ribosomal RNA" --sort -f ";~~~product~~~="
"""

##################################
### Loading Necessary Packages ###
##################################
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
def read_gff_dict(input):

    # limit_info = dict(gff_id=["chr1"], gff_source=["Coding_transcript"])
    # GFF.parse(in_handle, limit_info=limit_info)

    in_handle = open(input)
    return GFF.parse(in_handle)
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
def filter_exact_mode(input_gff):

    # Parse fields
    gff_dict = read_gff_dict(input_gff)

    # Check
    for rec in gff_dict:
        print(rec.features[0])

################
### Def main ###
################
def filter(input_gff, column, pattern, sort, header, mode):

    # Simple filter
    if mode == "loose":
        filter_loose_mode(input_gff=input_gff, column=column, pattern=pattern, sort=sort, header=header)

    # Complex filter
    elif mode == "exact":
        filter_exact_mode(input_gff=input_gff)

    # Error
    else:
        print(f"""
Error: --mode must be either 'loose' or 'exact'. {mode} is incorrect.
        """)
