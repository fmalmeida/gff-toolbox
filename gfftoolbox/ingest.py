## Def convert help
usage_ingest="""
gff-toolbox:

            Mongo-ingest

This command uses several python libraries to provide an easy way to convert your GFF data

usage:
    gff-toolbox mongo-ingest [ -h|--help ] [ --input <tsv> ] [--gff_feature gene --db_name <db_name> --genome_name <genome_name> --mongo_path <mongo_path> ]

options:
    -h, --help                                              Show this screen
    -i, --input=<tsv>                                       Annotation file in TSV (tab-separated values) format describing, for each line, the gff feature
                                                            (default: gene, can be changed in parameter --gff_feature) id and
                                                            corresponding annotations that should be included in mongodb database.
                                                            The annotation file must contain four columns (#locusName\tId\tIdType\tdescription). [Default: stdin].
    -l, --gff_feature=<feature_type>                        Which GFF feature type must be annotated. [Default: gene].
    -d, --db_name=<db_name>                                 Name of existing mongodb database to update with anotations.
                                                            If database doesnt exist, create it using gff-toolbox convert module [Default: annotation_db].
    -n, --genome_name=<genome_name>                         When loading the mongodb this will be used as collection name. [Default: Genome].
    -p, --mongo_path=<mongo_path>                           Where to load your mongoDB? [Default: ./mongodb].
                                                            If you insert a path that already have a mongoDB in it will include (append)
                                                            the GFF as new collection (<genome_name>) in a new or existing DB (<db_name>).


example:

    ## Converting a GFF to a mongoDB and updating annotation for genes listed in gene.functions.txt file
    ## This will be added to a mongo db called <db_name> in a collection named <genome_name>
    ## DBs are writen in the localhost 27027 mongo db connection of mongo shell

$ gff-toolbox convert --format mongodb -i Kp_ref.gff -n Kp

$ gff-toolbox mongo-ingest -i gene.functions.txt -n Kp

"""

##################################
### Loading Necessary Packages ###
##################################
import sys
import os
import re
import time
import tempfile
from collections import namedtuple
import gzip
import urllib.request, urllib.parse, urllib.error
import json
import pymongo
from pymongo import MongoClient
from pymongo import collection
from BCBio import GFF
from Bio import SeqIO
from io import StringIO
import binascii
import pathlib

######################################
### GFF columns names -- immutable ###
######################################
gff_cols = ["recid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

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

########################
### Useful functions ###
########################

# Check for stdin
def stdin_checker(input):
    # Checking for stdin
    if input == "stdin":
        tmp = tempfile.NamedTemporaryFile(mode = "w+t", delete = False) # Create tmp file to work as input
        temp_file = open(tmp.name, 'w+t')
        for line in sys.stdin:
            temp_file.writelines(f"{line}")

        temp_file.seek(0)
        input = tmp.name

    else:
        pass

    return input

# Print fasta format
def fasta_printer(id, seq):
    print(f">{id}\n{seq}")

# Get feature identifier
def tag_getter(record, seq):
    if record.id != "":
        tag = record.id
    else:
        tag = f"{record.type}_{seq.id}:{record.location.start}-{record.location.end}"

    return tag

# Open gzipped files
def gzip_opener(input, mode_in):
    if  binascii.hexlify(open(input, 'rb').read(2)) == b"1f8b" or input.endswith(".gz"):
        return gzip.open(input, mode=mode_in)
    else:
        return open(input, mode=mode_in)

# Remove nest from GFF
def _flatten_features(rec):
        """Make sub_features in an input rec flat for output.

        GenBank does not handle nested features, so we want to make
        everything top level.
        """
        out = []
        for f in rec.features:
            cur = [f]
            while len(cur) > 0:
                nextf = []
                for curf in cur:
                    out.append(curf)
                    if len(curf.sub_features) > 0:
                        nextf.extend(curf.sub_features)
                cur = nextf
        rec.features = out
        return rec


def ingestAttributes(filename, db_name, collection_name, mongo_path, feature_type, has_header):

    # Start mongo
    current = pathlib.Path().absolute()
    # os.system(f"mkdir -p {mongo_path} {mongo_path}/{db_name}")
    os.system(f"cd {mongo_path} && mongod --dbpath {db_name} --logpath mongo_{collection_name}.log &")

    # Create connection
    client = MongoClient()

    # Create Database
    db = client[db_name]

    # Create Collection
    collection = db[collection_name]

    # Create a compound index to further accelerate the update of collection entries
    # requires pymongo
    collection.create_index( [ ("type", 1),  ("attributes.ID", 1) ])


    # Open
    contents = gzip_opener(stdin_checker(filename), "rt")



    for idxline, line in enumerate(contents):
        
        if idxline == 0 and has_header: continue
        
        try:
                gene_name, iddb, idType, description = line.strip("\n").split("\t")
        except ValueError:
                print("Error: Could not unpack the values in line number {}: {}.\n \
                The file with the set of attributes to ingest into GFF must have 4 \
                 columns and be tab-delimited.".format(str(idxline), line.strip("\n").split("\t")))
        
        if (idxline % 1000) == 0:
            print("Update Col 9 N=", str(idxline))
        if idType != "GO":
            if len(description) == 0:
                obj = {"DBTAG": idType, "ID": iddb}
            else:
                obj = {"DBTAG": idType, "ID": iddb, "Description": description}
            
            collection.update_one(
                # query                 
                {'type': 'gene', 'attributes.ID': gene_name},
                # update
                {"$addToSet": {"attributes.Dbxref": obj}}
                # update options (Optional)
            )
        else:
            if len(description) == 0:
                obj = {"DBTAG": "GO", "ID": iddb}
            else:
                obj = {"DBTAG": "GO", "ID": iddb, "Description": description}
            collection.update_one(
                {'type': 'gene', 'attributes.ID': gene_name},
                {"$addToSet": {"attributes.Ontology_term": obj}}
            )

    # Close
    client.close()
    os.system(f"cd {current}")



################
### Def main ###
################
# (filename, db_name, collection_name, mongo_path, feature_type, has_header)
def ingest(filename, feature_type, db_name, genome_name, mongo_path):

    ingestAttributes(filename=filename, db_name=db_name, collection_name=genome_name, mongo_path=mongo_path, feature_type=feature_type, has_header=True)
