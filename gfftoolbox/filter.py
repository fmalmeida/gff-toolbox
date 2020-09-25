## Def filter help
usage_filter="""
gff-toolbox:

            Filter

This command uses Biopython library to parse and filter your GFF file as you wish. It targets the
attributes column of the GFF.

usage:
    gff-toolbox filter [ -h|--help ] [ --input <gff> ] [ --pattern <string> | --list <file> ] [ --field <string> --column <int> --sort ]

options:
    -h --help                                               Show this screen

    -i, --input=<gff>                                       Input GFF file. GFF file must not contain sequences with it

    -c, --column=<int>                                      Apply pattern search in which GFF columns? [Default: 9]

    -p, --pattern=<string>                                  Pattern to search in the GFF file

    -l, --list=<file>                                       List of patterns to search in the GFF file

    -f, --field=<string>                                    Useful when you want to check for a pattern in a specific field of a given column.
                                                            Format: -f delimiter~~~field~~~val-sign. Remember to use quotation. Read the example below

    --sort                                                  Sort the GFF by the contig and start position

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

####################################
### Function to import gff as df ###
####################################
def read_gff(input):
    return pd.read_csv(input, sep = "\t", comment = "#",
                        names=['1', '2', '3', '4', '5', '6', '7', '8', '9'])

# Guide: ['Chr', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'

###################################
### Filter df by column pattern ###
###################################
def df_col_filter(df, column, pattern):
    return df[
        df[str(column)].str.contains(f"{pattern}")
        ].to_csv(sep='\t', index=False, header=False).strip()

##############################################
### Function for filtering 9th column list ###
##############################################
def Filter(string, substr):
    return [str for str in string if
             any(sub in str for sub in substr)]

######################################################
### Function for simple filter with single pattern ###
######################################################
def simple_filter_gff_pattern(input_gff, column, pattern, sort):

    # Read GFF file
    df = read_gff(input_gff)

    if int(column) == 9:
        print("wait")

    else:

        # Sort
        if sort == True:
            df = df.sort_values(by=['1', '4'])

        # header
        print("##gff-version 3")

        # Filter
        print(df_col_filter(df, str(column), str(pattern)))

#######################################################
### Function for complex filter with single pattern ###
#######################################################
def complex_filter_gff_pattern(input_gff, column, pattern, sort, field):

    # Parse fields
    delimiter,field,sign = str(field).split('~~~')

    # Check
    print(f"""
Delimiter is: {delimiter}

Field is: {field}

Sign is: {sign}
    """)

### Def main
def filter(input_gff, column, pattern, sort, field):

    # Simple filter
    if field == None:
        simple_filter_gff_pattern(input_gff, column, pattern, sort)

    # Complex filter
    elif field != None:
        complex_filter_gff_pattern(input_gff, column, pattern, sort, field)
