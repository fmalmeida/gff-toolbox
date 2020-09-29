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
    overview                                                This command is very useful to get the gist of your GFF file
    filter                                                  This command lets you filter your GFF based on one or multiple patterns

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
            check_gff(args_overview['--input'])

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
                   end_pos=args_filter['--end'], strand=args_filter['--strand'])

        else:
            print(usage_filter.strip())

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
