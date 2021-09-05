.. _overview:

GFF-toolbox overview
====================

This command uses `BCBio package <https://biopython.org/wiki/GFF_Parsing>`_ to get the gist of an input GFF file.

Help
----

.. code-block:: bash

    # Trigger help
    gff-toolbox overview -h

    # Help
    gff-toolbox:

                Overview

    This command uses BCBio library to let you get the gist of your GFF file

    usage:
        gff-toolbox overview [ -h|--help ] [ --input <gff> ]

    options:
        -h, --help                        Show this screen
        -i, --input=<gff>                 Input GFF file. GFF file must not contain sequences with it. [Default: stdin]

    example:

        ## Getting the overview of a generic GFF file

    $ gff-toolbox overview -i Athaliana_ref.gff.gz

    $ gff-toolbox overview -i Kp_ref.gff

Example
-------

We can use as an example `the GFF <https://github.com/fmalmeida/gff-toolbox/raw/master/test/Kp_ref.gff>`_ provided in the test directory.

.. code-block:: bash

    # Example
    gff-toolbox overview -i Kp_ref.gff

    # Output
    Input: Kp_ref.gff gff file.

    Overview:

    {
        'gff_id': {
            ('NC_016838.1', ): 290,
            ('NC_016839.1', ): 305,
            ('NC_016840.1', ): 9,
            ('NC_016841.1', ): 3,
            ('NC_016845.1', ): 10982,
            ('NC_016846.1', ): 321,
            ('NC_016847.1', ): 11,
        },
        'gff_source': {('RefSeq', ): 11921},
        'gff_source_type': {
            ('RefSeq', 'CDS'): 5779,
            ('RefSeq', 'exon'): 88,
            ('RefSeq', 'gene'): 5867,
            ('RefSeq', 'rRNA'): 25,
            ('RefSeq', 'region'): 7,
            ('RefSeq', 'sequence_feature'): 92,
            ('RefSeq', 'tRNA'): 62,
            ('RefSeq', 'tmRNA'): 1,
        },
        'gff_type': {
            ('CDS', ): 5779,
            ('exon', ): 88,
            ('gene', ): 5867,
            ('rRNA', ): 25,
            ('region', ): 7,
            ('sequence_feature', ): 92,
            ('tRNA', ): 62,
            ('tmRNA', ): 1,
        },
    }


