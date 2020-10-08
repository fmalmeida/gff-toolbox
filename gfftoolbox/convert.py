## Def convert help
usage_convert="""
gff-toolbox:

            Convert

This command uses several python libraries to provide an easy way to convert your GFF data

usage:
    gff-toolbox convert [ -h|--help ]
    gff-toolbox convert [ --input <gff> --format json --genome_name <genome_name> ]
    gff-toolbox convert [ --input <gff> --format mongodb --db_name <db_name> --genome_name <genome_name> --mongo_path <mongo_path> ]
    gff-toolbox convert [ --input <gff> --format fasta --fasta <genome_file> --fasta_features <feature_types> --translation_table <int> ]
    gff-toolbox convert [ --input <gff> --format genbank --fasta <genome_file> --translation_table <int> ]

options:
                                                      Mandatory

    -h --help                                               Show this screen
    -i, --input=<gff>                                       Input GFF file. GFF file must not contain recuences with it. [Default: stdin].
    -f, --format=<out_format>                               Convert to which format? Options: json, mongodb, fasta-nt, fasta-aa, genbank. [Default: genbank]

                                                Converting to JSON and/or mongoDB
                        Obs: converted dbs are automatically added to a db in localhost 27027 ... dependent of mongo shell"

    -d, --db_name=<db_name>                                 Name of mongodb database to save results. Only for mongoDB. [Default: annotation_db].
    -n, --genome_name=<genome_name>                         Genome name. This will be used as main dictionary key to store the information of JSON.
                                                            When converting to mongodb this will be used as collection name. [Default: Genome].
    -p, --mongo_path=<mongo_path>                           Where to save your mongoDB? [Default: ./mongodb].
                                                            If you insert a path that already have a mongoDB in it will include (append)
                                                            the GFF as new collection (<genome_name>) in a new or existing DB (<db_name>).

                                                Converting to FASTA or Genbank

    --fasta=<genome_file>                                   Genomic fasta file to extract features from.
    -l, --fasta_features=<feature_types>                    A comma separated list of feature types to extract. Must be identical as written in the GFF file. [Default: CDS]
    -t, --translation_table=<int>                           NCBI's translation table number. For converting nucleotide sequences to protein.
                                                            Read more at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. [Default: 1]

example:

    ## Converting a GFF to JSON

$ gff-toolbox convert -i Kp_ref.gff -f json

    ## Converting a GFF to a mongoDB
    ## This will be added to a mongo db called <db_name> in a collection named <genome_name>
    ## DBs are writen in the localhost 27027 mongo db connection of mongo shell

$ gff-toolbox convert --format mongodb -i Kp_ref.gff

    ## Get CDS sequences from GFF to protein fasta

$ gff-toolbox convert -i Kp_ref.gff -f fasta-aa --fasta Kpneumoniae_genome.fasta -t 11

    ## Get mRNA sequences from GFF, in nucleotide fasta

$ gff-toolbox convert -i Athaliana_ref.gff.gz --fasta Athaliana_genome.fasta.gz --format fasta-nt --fasta_features mRNA

    ## Converting a GFF to GenBank

$ gff-toolbox convert -i Athaliana_ref.gff.gz --fasta Athaliana_genome.fasta.gz

Obs: The smaller the input, the fastest the convertion.
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

#######################
### Convert to JSON ###
#######################
def gff2json(filename):

    gff_dict = {} # Initialize gff as dict
    final    = [] # Final list for JSON

    # Open
    contents = gzip_opener(stdin_checker(filename), "rt")

    for line in contents:
        if line.startswith("#"): continue
        parts = line.strip().split("\t")
        #If this fails, the file format is not standard-compatible
        assert len(parts) == len(gff_cols)
        #Separate Values
        rec        = parts[0]
        source     = parts[1]
        feature    = parts[2]
        start      = parts[3]
        end        = parts[4]
        score      = parts[5]
        strand     = parts[6]
        phase      = parts[7]
        attributes = att_to_dict(parts[8])
        gff_dict = {
            # "CDS" : attributes[0],
            "recid"     : rec,
            "source"    : source,
            "type"      : feature,
            "start"     : start,
            "end"       : end,
            "score"     : score,
            "strand"    : strand,
            "phase"     : phase,
            "attributes": attributes
        }

        final.append(gff_dict)

    # json.dump(final, sys.stdout, indent=4)
    final = json.dumps(final, indent=4) # Converts to string
    return final

##########################
### Convert to mongoDB ###
##########################
def gff2mongo(filename, db_name, collection_name, mongo_path):

    # Start mongo
    current = pathlib.Path().absolute()
    os.system(f"mkdir -p {mongo_path} {mongo_path}/{db_name}")
    os.system(f"cd {mongo_path} && mongod --dbpath {db_name} --logpath mongo_{collection_name}.log &")

    # Create connection
    client = MongoClient()

    # Create Database
    db = client[db_name]

    # Create Collection
    collection = db[collection_name]

    # Save JSON as tmp.file
    tmp = tempfile.NamedTemporaryFile(mode = "w+t", delete = False) # Create tmp file to work as input
    temp_file = open(tmp.name, 'w+t')
    for line in gff2json(filename).strip():
        temp_file.writelines(f"{line}")

    temp_file.seek(0)
    input = tmp.name

    # Add document to collection


    contents = json.load(open(input))
    collection.insert_many(contents)

    # Close
    client.close()
    os.system(f"cd {current}")

########################
### Convert to FASTA ###
########################
def gff2fasta(input, fasta, features, format, translation_table):

    # Parse fasta
    rec_dict = SeqIO.to_dict(SeqIO.parse(gzip_opener(fasta, "rt"), "fasta"))

    # Read the gff
    for seq in GFF.parse(gzip_opener(stdin_checker(input), "rt"), base_dict=rec_dict):

        genome = seq.seq

        for rec in _flatten_features(seq).features:

            if bool(re.search("|".join(list(features.split(','))).lower(), str(rec.type.lower()))):
                tag = tag_getter(rec, seq=seq)
                if format == "fasta-aa":
                    fasta_printer(tag, genome[rec.location.start:rec.location.end].translate(table=translation_table))
                else:
                    fasta_printer(tag, genome[rec.location.start:rec.location.end])
            else:
                pass

##########################
### Convert to Genbank ###
##########################
def gff2gbk(input, fasta, translation_table):

    def _seq_getter(rec, fasta, translation_table):

        sequence = fasta[rec.location.start:rec.location.end]

        rec.seq = sequence

        if not bool(re.search("region", str(rec.type.lower()))):

            if bool(re.search("cds|protein", str(rec.type.lower()))):
                rec.qualifiers['translation']  = sequence.translate(table=translation_table)
            else:
                rec.qualifiers['sequence'] = sequence

        return rec

    genome = SeqIO.to_dict(SeqIO.parse(gzip_opener(fasta, "rt"), "fasta"))

    for seq in GFF.parse(gzip_opener(stdin_checker(input), "rt"), base_dict=genome):

        seq.annotations["molecule_type"] = "DNA"

        out = []

        for f in _flatten_features(seq).features:
            record = _seq_getter(f, seq.seq, translation_table)
            record.name = seq.name
            out.append(record)
        seq.features = out
        SeqIO.write(seq, sys.stdout, "genbank")

################
### Def main ###
################
def convert(filename, format, fasta, fasta_features, translation_table, db_name, genome_name, mongo_path):

    if format == "json" :
        print(gff2json(filename=filename).replace("[", "{").replace("]", "}"))

    elif format == "mongodb" :
        gff2mongo(filename=filename, db_name=db_name, collection_name=genome_name, mongo_path=mongo_path)

    elif format == "fasta-nt" or format == "fasta-aa":
        gff2fasta(input=filename, fasta=fasta, features=fasta_features,
                  format=format, translation_table=translation_table)
    elif format == "genbank":
        gff2gbk(input=filename, fasta=fasta, translation_table=translation_table)

    else:
        print(f"""
Error: I can't understand {format} format ... Please check the help.
        """)
