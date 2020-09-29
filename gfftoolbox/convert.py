## Def convert help
usage_convert="""
gff-toolbox:

            Convert

This command uses several python libraries to provide an easy way to convert your GFF data

usage:
    gff-toolbox convert [ -h|--help ]
    gff-toolbox convert [ --input <gff> --format <out_format> --id <field_identifier> ]

options:
    -h --help                                               Show this screen
    -i, --input=<gff>                                       Input GFF file. GFF file must not contain sequences with it.
    -f, --format=<out_format>                               Convert to which format? Options: json
    --id=<field_identifier>                                 Field name of the attribute in the 9th column to get feature IDs/names from. [Default: ID]

example:

$ gff-toolbox convert -i test/input.gff -f json
"""

##################################
### Loading Necessary Packages ###
##################################
import sys
import time
from collections import namedtuple
import gzip
import urllib.request, urllib.parse, urllib.error
import json

######################################
### GFF columns names -- immutable ###
######################################
gff_cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

#############################################
### Parse attribute columns as dictionary ###
#############################################
def att_to_dict(attributes):

    final = {} # Initialize dictionary
    for atts in attributes.split(";"):
        try:
            key, value = atts.split("=")
            final[key] = value.strip()
        except ValueError:
            print(f"Can't split attributes. Problem in {atts}")
            break

    return final

#######################
### Convert to JSON ###
#######################
def gff2json(filename, identifier):

    gff_dict = {} # Initialize gff as dict

    openFunc = gzip.open if filename.endswith(".gz") else open # Parse with transparent decompression

    with openFunc(filename) as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            # If this fails, the file format is not standard-compatible
            assert len(parts) == len(gff_cols)
            # Separate Values
            seq        = parts[0]
            source     = parts[1]
            feature    = parts[2]
            start      = parts[3]
            end        = parts[4]
            score      = parts[5]
            strand     = parts[6]
            phase      = parts[7]
            attributes = att_to_dict(parts[8])
            gff_dict[attributes[identifier]] = {
                "seqid"     : seq,
                "source"    : source,
                "type"      : feature,
                "start"     : start,
                "end"       : end,
                "score"     : score,
                "strand"    : strand,
                "phase"     : phase,
                "attributes": attributes
            }

    json.dump(gff_dict, sys.stdout, indent=4)

################
### Def main ###
################
def convert(filename, format, identifier):

    if format == "json" :
        gff2json(filename=filename, identifier=identifier)

    else:
        print(f"""
Error: I can't understand {format} format ... Please check the help.
        """)
