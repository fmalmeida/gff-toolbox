## Def convert help
usage_convert="""
gff-toolbox:

            Convert

This command uses several python libraries to provide an easy way to convert your GFF data

usage:
    gff-toolbox convert [ -h|--help ]
    gff-toolbox convert [ --input <gff> --format json ]
    gff-toolbox convert [ --input <gff> --format fasta --fasta <genome_file> --fasta_features <feature_types> --translation_table <int>]
    gff-toolbox convert [ --input <gff> --format genbank --fasta <genome_file> --translation_table <int>]

options:
    -h --help                                               Show this screen
    -i, --input=<gff>                                       Input GFF file. GFF file must not contain recuences with it. [Default: stdin].
    -f, --format=<out_format>                               Convert to which format? Options: json,fasta-nt, fasta-aa, genbank. [Default: genbank]

                                                Converting to FASTA or Genbank

    --fasta=<genome_file>                                   Genomic fasta file to extract features from.
    -l, --fasta_features=<feature_types>                        A comma separated list of feature types to extract. Must be identical as written
                                                            in the GFF file. [Default: CDS]
    -t, --translation_table=<int>                           NCBI's translation table number. For converting nucleotide sequences to protein.
                                                            Read more at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. [Default: 1]

example:

    ## Converting a GFF to JSON

$ gff-toolbox convert -i Kp_ref.gff -f json

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
import time
from collections import namedtuple
import gzip
import urllib.request, urllib.parse, urllib.error
import json
from BCBio import GFF
from Bio import SeqIO
from io import StringIO
import binascii

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
        temp_file = open(tmp.name, 'w')
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
        open(input, mode=mode_in)

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

    # Check if is stdin
    filename = stdin_checker(filename)

    # Open
    contents = gzip_opener(filename, "w+t")

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
    final = final.replace("[", "{").replace("]", "}")
    print(final)

########################
### Convert to FASTA ###
########################
def gff2fasta(input, fasta, features, format, translation_table):

    # Checking for stdin
    input = stdin_checker(input)

    # Subset GFF based on chr and feature type
    limit_info = dict(
            gff_type = list(features.split(','))
    )
    features = limit_info['gff_type']
    # Parse fasta
    rec_dict = SeqIO.to_dict(SeqIO.parse(gzip_opener(fasta, "rt"), "fasta"))

    # Read the gff
    for seq in GFF.parse(gzip_opener(input, "rt"), base_dict=rec_dict, limit_info=limit_info):

        genome = seq.seq

        for rec in _flatten_features(seq).features:

            if rec.type in features:
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

    # Checking for stdin
    input = stdin_checker(input)

    def _seq_getter(rec, fasta, translation_table):

        sequence = fasta[rec.location.start:rec.location.end]

        rec.seq = sequence

        if rec.type.lower() != "region":
            if str(rec.type).lower() == "cds" or str(rec.type).lower() == "protein":
                rec.qualifiers['translation']  = sequence.translate(table=translation_table)
            else:
                rec.qualifiers['sequence'] = sequence

        return rec

    genome = SeqIO.to_dict(SeqIO.parse(gzip_opener(fasta, "rt"), "fasta"))

    for seq in GFF.parse(gzip_opener(input, "rt"), base_dict=genome):

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
def convert(filename, format, fasta, fasta_features, translation_table):

    if format == "json" :
        gff2json(filename=filename)

    elif format == "fasta-nt" or format == "fasta-aa":
        gff2fasta(input=filename, fasta=fasta, features=fasta_features,
                  format=format, translation_table=translation_table)
    elif format == "genbank":
        gff2gbk(input=filename, fasta=fasta, translation_table=translation_table)

    else:
        print(f"""
Error: I can't understand {format} format ... Please check the help.
        """)
