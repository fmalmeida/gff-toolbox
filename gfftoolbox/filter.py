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

########################################
### Function for execution of parser ###
########################################
def filter_gff_pattern(input_gff, column, pattern, sort):

    # Read GFF file
    df = read_gff(input_gff)

    if int(column) == 9:
        print("wait")

    else:

        # Sort
        if sort == True:
            df = df.sort_values(by=['4'])

        # header
        print("##gff-version 3")

        # Filter
        print(df_col_filter(df, str(column), str(pattern)))
