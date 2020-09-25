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
                        names=['Chr', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'])

##############################################
### Function for filtering 9th column list ###
##############################################
def Filter(string, substr):
    return [str for str in string if
             any(sub in str for sub in substr)]

########################################
### Function for execution of parser ###
########################################
def filter_gff_pattern(input_gff):

    # Read GFF file
    my_df = read_gff(input_gff)

    # Check
    print(my_df.head())
