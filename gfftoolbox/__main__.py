#!/usr/bin/env python3
license="""
Copyright 2020 Felipe Almeida (almeidafmarques@gmail.com)
https://github.com/fmalmeida/gff-toolbox

This file is part of gff-toolbox. gff-toolbox is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. gff-toolbox is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with gff-toolbox.
If not, see <http://www.gnu.org/licenses/>.
"""

## Def main help
usage="""
gff-toolbox: a compendium of commands to facilitate the work with GFF files

Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

usage:
    gff-toolbox [ -h|--help ] [ -v|--version ] [ --license ]
    gff-toolbox <command> [ -h|--help ] [ <args>... ]

options:
    -h --help                                               Show this screen
    -v --version                                            Show version information
    --license                                               Show LEGAL LICENSE information

commands:
    overview                                                Useful to get the gist of your GFF file
    filter                                                  It lets you filter your GFF based on one or multiple patterns
    convert                                                 Converts a GFF file into another formats
    plot                                                    Useful command to plot genomic regions from a GFF file
    mongo-ingest                                            Useful to include new or update annotations in existing GFF mongo databases

Use: `gff-toolbox <commmand> -h` to get more help and see examples.
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
from subprocess import call
import sys

########################
### Import functions ###
########################
from .overview import *
from .filter import *
from .version import *
from .convert import *
from .plot import *
from .ingest import *

## Defining main
def main():
    # Parse docopt
    __version__ = get_version()
    arguments = docopt(usage, version=__version__, help=False, options_first=True)

    ############################
    ### GFF overview command ###
    ############################
    if arguments['<command>'] == 'overview':
        
        # Parse docopt
        args_overview = docopt(usage_overview, version=__version__, help=False)
        
        if args_overview['overview'] and args_overview['--help']:
            print(usage_overview.strip())

        elif args_overview['overview'] and args_overview['--input']:
            overview(args_overview['--input'])

        else:
            print(usage_overview.strip())

    ##########################
    ### GFF filter command ###
    ##########################
    elif arguments['<command>'] == 'filter':

        # Parse docopt
        args_filter = docopt(usage_filter, version=__version__, help=False)

        if args_filter['filter'] and args_filter['--help']:
            print(usage_filter.strip())

        ## Run it
        elif args_filter['filter'] and args_filter['--input'] and not args_filter['--help']:
            filter(input_gff=args_filter['--input'], column=args_filter['--column'],
                   pattern=args_filter['--pattern'], sort=args_filter['--sort'],
                   header=args_filter['--header'], mode=args_filter['--mode'],
                   chr_limits=args_filter['--chr'], source_limits=args_filter['--source'],
                   type_limits=args_filter['--type'], start_pos=args_filter['--start'],
                   end_pos=args_filter['--end'], strand=args_filter['--strand'],
                   att_file=args_filter['--attributes'])

        else:
            print(usage_filter.strip())

    ###########################
    ### GFF convert command ###
    ###########################
    elif arguments['<command>'] == 'convert':

        # Parse docopt
        args_convert = docopt(usage_convert, version=__version__, help=False)

        if args_convert['convert'] and args_convert['--help']:
            print(usage_convert.strip())

        ## Run it
        elif args_convert['convert'] and args_convert['--input'] and args_convert['--format'] and not args_convert['--help']:
            convert(filename=args_convert['--input'], format=args_convert['--format'], fasta=args_convert['--fasta'],
                    fasta_features=args_convert['--fasta_features'], translation_table=args_convert['--translation_table'],
                    db_name=args_convert['--db_name'], genome_name=args_convert['--genome_name'], mongo_path=args_convert['--mongo_path'])

        else:
            print(usage_convert.strip())

    ###########################
    ### GFF ingest command ###
    ###########################
    elif arguments['<command>'] == 'mongo-ingest':

        # Parse docopt
        args_ingest = docopt(usage_ingest, version=__version__, help=False)

        if args_ingest['mongo-ingest'] and args_ingest['--help']:
            print(usage_ingest.strip())

        ## Run it
        elif args_ingest['mongo-ingest'] and args_ingest['--input'] and args_ingest['--db_name'] and not args_ingest['--help']:
            
            ingestAttributes(filename=args_ingest['--input'], feature_type=args_ingest['--gff_feature'],
             db_name=args_ingest['--db_name'], collection_name=args_ingest['--genome_name'],
              mongo_path=args_ingest['--mongo_path'], has_header=True)

        else:
            print(usage_ingest.strip())

    ########################
    ### GFF plot command ###
    ########################
    elif arguments['<command>'] == 'plot':

        # Parse docopt
        args_plot = docopt(usage_plot, version=__version__, help=False)

        ## Check GFF
        if args_plot['check-gff'] and args_plot['--input']:
            check_gff(args_plot['--input'])

        ## Single GFF
        if args_plot['--input'] and args_plot['--start'] and args_plot['--end'] and args_plot['--contig'] and not args_plot['--fofn']:
            print("Executing the pipeline for a single GFF input")
            single_gff(infile=args_plot['--input'], start=args_plot['--start'], end=args_plot['--end'],
                       contig=args_plot['--contig'], feature=args_plot['--feature'], coloring=args_plot['--color'],
                       custom_label=args_plot['--label'], outfile=args_plot['--output'], plot_title=args_plot['--title'],
                       qualifier=args_plot['--identification'], plot_width=args_plot['--width'], plot_height=args_plot['--height'])
            print("Done, checkout the results in {}".format(args_plot['--output']))

        ## Multiple GFFs
        elif args_plot['--fofn'] and args_plot['--start'] and args_plot['--end'] and args_plot['--contig']:
            print("Executing the pipeline for multiple GFF inputs")
            multiple_gff(input_fofn=args_plot['--fofn'], start=args_plot['--start'], end=args_plot['--end'],
                         contig=args_plot['--contig'], feature=args_plot['--feature'], outfile=args_plot['--output'],
                         qualifier=args_plot['--identification'], plot_title=args_plot['--title'],
                         plot_width=args_plot['--width'], plot_height=args_plot['--height'])
            print("Done, checkout the results in {}".format(args_plot['--output']))

        ## None
        else:
            print(usage_plot.strip())

    #####################
    ### Check license ###
    #####################
    elif arguments['--license']:
        print(license.strip())

    #######################################
    ### Without commands nor parameters ###
    #######################################
    else:
        print(usage.strip())

## Calling main
if __name__ == '__main__':
    main()
