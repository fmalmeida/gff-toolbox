#!/usr/bin/env python3
"""
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
    gff-toolbox [ -h|--help ] [ -v|--version ]
    gff-toolbox overview [ -h|--help ] [ -i <gff> ]
    gff-toolbox filter   [ -h|--help ] [ --input <gff> ] [ --patern <string> | --list <file> ] [ --field <string ]

options:
    -h --help                                               Show this screen
    -v --version                                            Show version information
    -i, --input=<gff>                                       Input GFF file. GFF file must not contain sequences with it
    -p, --pattern=<string>                                  Pattern to search in the GFF file
    -l, --list=<file>                                       List of patterns to search in the GFF file
    -c, --column=<int>                                      Apply pattern search in which GFF columns? [Default: 9]
    -f, --field=<string>                                    Apply pattern search in which column field? [Default: ID]

commands:
    overview                                                This command is very useful to get the gist of your GFF file
    filter                                                  This command lets you filter your GFF based on one or multiple patterns

Use: `gff-toolbox <commmand> -h` to get more help and see examples.
"""

## Def overview help
usage_overview="""
gff-toolbox:

            Overview

This command uses BCBio library to let you get the gist of your GFF file

usage:
    gff-toolbox overview [ -h|--help ] [ --input <gff> ]

options:
    -h --help                                               Show this screen
    -i, --input=<gff>                                       Input GFF file. GFF file must not contain sequences with it

example:
    gff-toolbox overview -i input.gff
"""

## Def overview help
usage_filter="""
gff-toolbox:

            Filter

This command uses Biopython library to parse and filter your GFF file as you wish. It targets the
attributes column of the GFF.

usage:
    gff-toolbox filter [ -h|--help ] [ --input <gff> ] [ --patern <string> | --list <file> ] [ --field <string ]

options:
    -h --help                                               Show this screen
    -i, --input=<gff>                                       Input GFF file. GFF file must not contain sequences with it
    -p, --pattern=<string>                                  Pattern to search in the GFF file
    -l, --list=<file>                                       List of patterns to search in the GFF file
    -c, --column=<int>                                      Apply pattern search in which GFF columns? [Default: 9]
    -f, --field=<string>                                    Only useful if applying the pattern search in the 9th column,
                                                            the attributes column. This parameter sets in which of the fields
                                                            of the 9th column to apply the pattern search [Default: ID]

example:
    gff-toolbox filter -i input.gff
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

## Defining main
def main():
    arguments = docopt(usage, version='v0.1 by Felipe Almeida', help=False)

    ## GFF overview
    if arguments['overview']:

        if arguments['overview'] and arguments['--help']:
            print(usage_overview.strip())

        elif arguments['overview'] and arguments['--input']:
            check_gff(arguments['--input'])

        else:
            print(usage_overview.strip())

    ## GFF filter
    if arguments['filter']:

        if arguments['filter'] and arguments['--help']:
            print(usage_filter.strip())

        elif arguments['filter'] and arguments['--input']:
            filter_gff_pattern(arguments['--input'])

        else:
            print(usage_filter.strip())

    ## None
    else:
        print(usage.strip())

## Calling main
if __name__ == '__main__':
    main()
