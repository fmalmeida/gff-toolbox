## Def filter help
usage_filter="""
gff-toolbox:

            Filter

This command uses Biopython library to parse and filter your GFF file as you wish. It targets the
attributes column of the GFF.

usage:
    gff-toolbox filter [-h|--help ]
    gff-toolbox filter [ --mode loose ] [ --input <gff> ] [ --pattern <string> --column <int> --start <start_position> --end <end_position> --strand <strand> --sort --header ]
    gff-toolbox filter [ --mode exact ] [ --input <gff> ] [ --chr <chr_limits> --source <source_limits> --type <type_limits> --start <start_position> --end <end_position> --strand <strand> --attributes <file_with_attributes> ]

options:

                                                Generic parameters

    -h --help                                               Show this screen.

    -i, --input=<gff>                                       Input GFF file. GFF file must not contain sequences, only features. [Default: stdin]

    -m, --mode=<search_mode>                                In which mode to search for patterns: loose or exact?
                                                            The loose mode, scans the GFF in a grep-like manner via pandas dataframes in which the user must specify
                                                            a pattern and a column to search it. Recommended for simple searches were nest structure is not a must.
                                                            The exact mode scans the GFF with Biopython and BCBio packages, treating it as python dictionary. It is
                                                            recommended for more complex searches and complex GFFs, such as nested GFFs. [Default: exact]

    --strand=<strand>                                       Apply a filter based on the strand of the feature. Options: plus or minus. By default, everything is given.
                                                            In exact mode, this filter is applied in the parent feature, if it passes, it's children are also printed.
                                                            The contrary is also true. In the loose mode it is applied directly to all features, nested or not.

    --start=<start_position>                                Apply a filter to select features starting from this position. In exact mode, this filter is applied in the
                                                            parent feature, if it passes, it's children are also printed. The contrary is also true. In the loose mode
                                                            it is applied directly to all features, nested or not.

    --end=<end_position>                                    Apply a filter to select features until this position. In exact mode, this filter is applied in the parent
                                                            feature, if it passes, it's children are also printed. The contrary is also true.

                                                Loose search mode parameters (Handy in general cases)

    -c, --column=<int>                                      Apply pattern search in which GFF columns?. [Default: 9]

    -p, --pattern=<string>                                  Pattern to search in the GFF file. Can be a list of patterns separated by commas.

    --sort                                                  Sort the GFF by the contig and start position. Be aware, it can disorganize nested gffs.

    --header                                                Print GFF header (##gff-version 3)? Some programs require this header.

                                                Exact search mode parameters (Very useful for nested GFFs)

    --chr=<chr_limits>                                      Apply a filter based on the chr/contig/sequence ids (Column 1). Can be a list of patterns separated by commas.
                                                            This step only works using the complete string for full-matches (it does not work with partial-matches based
                                                            substrings of the desired pattern).

    --source=<source_limits>                                Apply a filter based on the source column (Column 2). Can be a list of patterns separated by commas.
                                                            This step works using the complete string (with full-matches) or substrings of the desired pattern,
                                                            working with partial-matches.

    --type=<type_limits>                                    Apply a filter based on the type column (Column 3). Can be a list of patterns separated by commas.
                                                            This step works using the complete string (with full-matches) or substrings of the desired pattern,
                                                            working with partial-matches. In the loose mode it is applied directly to all features, nested or not.

    --attributes=<file_with_attributes>                     Pass a file containing the desired key/value tuple to search in the 9th column. The header of the file is the
                                                            attribute key in which to search for the values given in the following it. Since it maintains the nest and
                                                            organization of the file, it is useful for filtering nested GFFs based on a list of genes, parents or products.
                                                            The maintainence of the nest structure would be difficult to have with simpler commands such as `grep -f filep`
                                                            since children and parents seldom have the same attribute keys.

                                                            This file must a header starting with '##', whithout space and its values following it. E.g.:

                                                                    ##ID
                                                                    desired gene id 1
                                                                    desired gene id 2
                                                                    ...


example:

    ## Simple filter in any column: wheter a line contain a pattern in a specific column (like grep)
    ## Check the features that have the word "putative" in their attributes.

$ gff-toolbox filter --mode loose --sort --header -i Kp_ref.gff -p "putative"

    ## In the example below, we filter the GFF in a more complex manner:
    ## All the CDS(s) found in the sequence named NC_016845.1 that
    ## have the word "transcriptional regulator" in their attributes.
    ##
    ## It works in both ways:

$ gff-toolbox filter -i Kp_ref.gff --chr NC_016845.1 --type CDS | gff-toolbox filter --mode loose -p "transcriptional regulator"
$ gff-toolbox filter --mode loose -i Kp_ref.gff -p "transcriptional regulator" | gff-toolbox filter --chr NC_016845.1 --type CDS

    ## Filtering a set of genes and its childs using a file containing the desired attributes.
    ## K. pneumoniae annotation.

$ gff-toolbox filter -i Kp_ref.gff --attributes atts.txt

    ## Filtering a set of genes and its childs using a file containing the desired attributes.
    ## A. thaliana annotation. Also give a custom start position for features to be printed?

$ gff-toolbox filter -i Athaliana_ref.gff --attributes atts2.txt --start 5900
"""

##################################
### Loading Necessary Packages ###
##################################
import sys
import pandas as pd
import re
from BCBio import GFF
from BCBio.GFF import GFFExaminer
import subprocess
import tempfile
import binascii
import gzip

#######################
### Def error class ###
#######################
class Error(Exception):
    """Base class for other exceptions"""
    pass

class noChild(Error):
    """it has no childs"""
    pass

class goChild(Error):
    """let's check the child"""
    pass

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
    if  binascii.hexlify(open(input, 'rb').read(2)) == b"1f8b" or input.endswith(".gz"):
        return gzip.open(input, mode=mode_in)
    else:
        open(input, mode=mode_in)

####################################
### Function to import gff as df ###
####################################
def read_gff_df(input, strand, start_pos, end_pos):

    # Guide: ['Chr', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'
    df = pd.read_csv(input, sep = "\t", comment = "#", names=['1', '2', '3', '4', '5', '6', '7', '8', '9'])

    ## Wants plus strand
    if strand == "plus":
        df = df[ df['7'] != "-" ]

    ## Wants minus strand
    elif strand == "minus":
        df = df[ df['7'] != "+" ]

    # Remove features based on position
    ## Min (start)
    if start_pos != None:
        df = df [ df['4'] >= int(start_pos) ]

    ## Max (end)
    if end_pos != None:
        df = df[ df['5'] <= int(end_pos) ]

    return df

#####################################################
### Function to check common member between lists ###
#####################################################
def common_member(a, b):

    result = not set(a).isdisjoint(set(b))
    return result

######################################
### Function to import gff as dict ###
######################################
def read_gff_dict(input, chr_limits, source_limits, type_limits, strand, start_pos, end_pos, att_file):

    # Check for the limits imposed by the user
    limit_info = {}

    # Chr limit?
    if chr_limits != None:
        chr_list = list(chr_limits.split(','))
        limit_info['gff_id'] = chr_list

    # GFF limits
    ## Module
    examiner = GFFExaminer()

    ## Open connections
    summary = examiner.available_limits(open(input))
    gff_sources = []
    gff_types = []

    ## Get sources
    for keys in summary['gff_source'].keys():
        key = str(keys).replace("('", "").replace("',)", "")
        gff_sources.append(key)

    ## Get types
    for keys in summary['gff_type'].keys():
        key = str(keys).replace("('", "").replace("',)", "")
        gff_types.append(key)

    # Source limits?
    if source_limits != None:
        src_definitive  = [] # The real patterns in the GFF
        src_list = list(source_limits.split(',')) # Patterns given by user
        for source in gff_sources:
            if bool(re.search('|'.join(src_list), str(source))): # Check wheter the pattern given by user is present in the GFF
                src_definitive.append(str(source)) # Select the real pattern that have the pattern given by user, for biopython dict

            limit_info['gff_source'] = list(set(src_definitive))

    # Type limits?
    if type_limits != None:
        type_definitive = [] # The real patterns in the GFF
        type_list = list(type_limits.split(',')) # Patterns given by user
        for type in gff_types:
            if bool(re.search('|'.join(type_list), str(type))): # Check wheter the pattern given by user is present in the GFF
                type_definitive.append(str(type)) # Select the real pattern that have the pattern given by user, for biopython dict

        limit_info['gff_type'] = list(set(type_definitive))

    # Open GFF for more customisable / exact filters
    records = []

    for rec in GFF.parse(open(input), limit_info=limit_info):

        ###############################################
        ### Complex search: loop in attributes nest ###
        ### Remove features based on attributes it  ###
        ### tries to maintain the nestness of the   ###
        ### GFF, which means, it searchs the value  ###
        ### in the first, second and third level of ###
        ### the nest, respectively. If found in     ###
        ### first, all its childs are given; If     ###
        ### found in second, its parent and its     ###
        ### childs are given; If found in the third,###
        ### its parents and grand-parents are given ###
        ###############################################

        if att_file != None:

            indexes = [] # Indexes to be saved by this filter -- only applies if given by user

            att_filter = {} # Store data from file
            for line in open(att_file, "r"):
                line = line.strip()
                if line.startswith("#"): # Is a header?
                    key = line.replace("##", "")
                    att_filter[key] = []
                else:
                    att_filter[key].append(line)

            # Parse GFF based on file data
            for field in att_filter:

                for index, f in enumerate(rec.features):

                    #############################
                    ### PARENT / MAIN feature ###
                    #############################

                    # Checking the value of the attribute key in the main feature
                    try:
                        if len(f.qualifiers[field]) > 0:

                            if common_member(f.qualifiers[field], list(set(att_filter[field]))):
                                indexes.append(int(index)) # Save this feature and its childs.

                            else:
                                # Despite having the same attribute-key, the value is not what the
                                # user wants. Maybe it refers to the same key in a child feature?
                                raise goChild

                    except (goChild, KeyError):

                        ###############################
                        ### 1st level child feature ###
                        ###############################

                        sub_index_list = [] # Create a list for saved indexes

                        # Checking the value of the attribute key in the 1st level child feature
                        for sub_index, sub in enumerate(f.sub_features):

                            try:

                                if len(sub.qualifiers[field]) > 0:

                                    if common_member(sub.qualifiers[field], list(set(att_filter[field]))):
                                        indexes.append(int(index)) # Save the main feature
                                        sub_index_list.append(int(sub_index)) # Save the subfeature

                                    else:
                                        # Despite having the same attribute-key, the value is not what the
                                        # user wants. Maybe it refers to the same key in a child feature?
                                        raise goChild

                            except(goChild, KeyError):

                                ###############################
                                ### 2nd level child feature ###
                                ###############################

                                subsub_index_list = [] # Create a list for saved indexes

                                # Checking the value of the attribute key in the 2nd level child feature
                                for subsub_index, subsub in enumerate(sub.sub_features):

                                    # Inside the 2nd level child

                                    try:

                                        if len(subsub.qualifiers[field]) > 0:

                                            if common_member(subsub.qualifiers[field], list(set(att_filter[field]))) == True:
                                                indexes.append(int(index)) # Save the main feature
                                                sub_index_list.append(int(sub_index)) # Save the subfeature
                                                subsub_index_list.append(int(subsub_index)) # Save this subsubfeature

                                            else:
                                                pass

                                    except:
                                        pass

                                # Finish the parsing for subsubfeatures
                                # Remove childs that were not selected -- after the loop in all childs
                                # outside the for loop. Change the subfeatures of this parent in the main object
                                rec.features[index].sub_features[sub_index].sub_features = [ i for j, i in enumerate(sub.sub_features) if j in list(set(subsub_index_list)) ]

                        # Finish the parsing for subfeatures
                        # Remove childs that were not selected -- after the loop in all childs
                        # outside the for loop. Change the subfeatures of this parent in the main object
                        rec.features[index].sub_features = [ i for j, i in enumerate(f.sub_features) if j in list(set(sub_index_list)) ]

            # Filter out the features that do not pass the filter
            rec.features = [i for j, i in enumerate(rec.features) if j in list(set(indexes))]

        ##############################################################################
        ### START of simpler filters, it is better to put it after the complex one ###
        ### based on the loop in attributes of parents and childs since it depends ###
        ### on it                                                                  ###
        ##############################################################################


        # Simpler filters: start, end, strand
        simpler_filtered_out_indexes = [] # Indexes that must be removed by the simpler filters

        # Remove features based on strand
        ## Wants plus strand
        if strand == "plus":
            for index, f in enumerate(rec.features):
                try:
                    if int(f.location.strand) == -1:
                        # Strand equals to the minus strand, thus remove
                        simpler_filtered_out_indexes.append(int(index))
                except:
                    # There is a problem ... the strand columns is not correct
                    # Thus we let it stay since we cannot assure is procedence
                    # And the user can check it.
                    pass

        ## Wants minus strand
        elif strand == "minus":
            for index, f in enumerate(rec.features):
                try:
                    if int(f.location.strand) == 1:
                        # Strand equals to the plus strand, thus remove
                        simpler_filtered_out_indexes.append(int(index))
                except:
                    # There is a problem ... the strand columns is not correct
                    # Thus we let it stay since we cannot assure is procedence
                    # And the user can check it.
                    pass

        # Remove features based on position
        ## Min (start)
        if start_pos != None:
            for index, f in enumerate(rec.features):
                if int(f.location.start) + 1 < int(start_pos): # Biopython is zero-based
                    simpler_filtered_out_indexes.append(int(index))

        ## Max (end)
        if end_pos != None:
            for index, f in enumerate(rec.features):
                if int(f.location.end) > int(end_pos):
                    simpler_filtered_out_indexes.append(int(index))

        # Filter out the features that do not passed the simpler filters
        rec.features = [i for j, i in enumerate(rec.features) if j not in list(set(simpler_filtered_out_indexes))]

        # Save record by appending
        records.append(rec)

    # Done
    return records


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
def filter_loose_mode(input_gff, column, pattern, sort, header, strand, start_pos, end_pos):

    # Read GFF file
    df = read_gff_df(input=input_gff, strand=strand, start_pos=start_pos, end_pos=end_pos)

    # Sort
    if sort == True:
        df = df.sort_values(by=['1', '4'])

    # header
    if header == True:
        print("##gff-version 3")

    # Filter
    if pattern != None:
        print(df_col_filter(df=df, column=str(column), pattern=str(pattern)))
    else:
        print(df.to_csv(sep='\t', index=False, header=False).strip())

#######################################################
### Function for complex filter with single pattern ###
#######################################################
def filter_exact_mode(input_gff, chr_limits, source_limits, type_limits, start_pos, end_pos, strand, att_file):

    # Parse fields
    gff_dict = read_gff_dict(input=input_gff, chr_limits=chr_limits, source_limits=source_limits, type_limits=type_limits,
                             strand=strand, start_pos=start_pos, end_pos=end_pos, att_file=att_file)

    # Print the records filtered (each record is a sequence)
    for record in gff_dict:
        if len(record.features) > 0:
            GFF.write([record], sys.stdout) # Write filtered gff

################
### Def main ###
################
def filter(input_gff, column, pattern, sort, header, mode, chr_limits, source_limits, type_limits, start_pos, end_pos, strand, att_file):

    # Checking for stdin
    input_gff = stdin_checker(input_gff)

    # Gzip?
    input_gff = gzip_opener(input_gff)

    # Simple filter
    if mode == "loose":
        filter_loose_mode(input_gff=input_gff, column=column, pattern=pattern, sort=sort, header=header,
                          start_pos=start_pos, end_pos=end_pos, strand=strand)

    # Complex filter
    elif mode == "exact":
        filter_exact_mode(input_gff=input_gff, chr_limits=chr_limits, source_limits=source_limits, type_limits=type_limits,
                          start_pos=start_pos, end_pos=end_pos, strand=strand, att_file=att_file)

    # Error
    else:
        print(f"""
Error: --mode must be either 'loose' or 'exact'. {mode} is incorrect.
        """)
