.. _convert:

GFF-toolbox convert
===================

This command uses several python libraries to provide an easy and simple way to explore and convert your annotations in GFF into other formats.

    Whenever using this command to extract protein sequences from the genome based on a GFF file, please understand that this command uses Biopython to do this extraction and convertion, not taking into account different reading frames ... Therefore, users may observe differences in the protein sequences if they already have a protein fasta file. If, you already have a fasta file of gene sequences and proteins, it is better to subset the data with `seqtk subseq <https://github.com/lh3/seqtk>`_.

Help
----

.. code-block:: bash

    # Trigger help
    gff-toolbox convert -h

    # Help
    gff-toolbox:

                Convert

    This command uses several python libraries to provide an easy way to convert your GFF data

    usage:
        gff-toolbox convert [ -h|--help ]
        gff-toolbox convert [ --input <gff> --fasta <genome_file> --translation_table <int> ] [ --format json|genbank ]
        gff-toolbox convert [ --input <gff> --fasta <genome_file> --translation_table <int> ] [ --format fasta --fasta_features <feature_types> ]
        gff-toolbox convert [ --input <gff> --fasta <genome_file> --translation_table <int> ] [--format mongodb --db_name <db_name> --genome_name <genome_name> --mongo_path <mongo_path> ]

    options:
                                                        General

        -h --help                                               Show this screen
        -i, --input=<gff>                                       Input GFF file. GFF file must not contain recuences with it. [Default: stdin].
        -f, --format=<out_format>                               Convert to which format? Options: json, mongodb, fasta-nt, fasta-aa, genbank. [Default: genbank]
        --fasta=<genome_file>                                   Genomic fasta file to extract features from.
        -t, --translation_table=<int>                           NCBI\'s translation table number. For converting nucleotide sequences to protein.
                                                                Read more at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. [Default: 1]

                                                    Converting to mongoDB
                            Obs: converted dbs are automatically added to a db in localhost 27027 ... dependent of mongo shell

        -d, --db_name=<db_name>                                 Name of mongodb database to save results. Only for mongoDB. [Default: annotation_db].
        -n, --genome_name=<genome_name>                         When converting to mongodb this will be used as collection name. [Default: Genome].
        -p, --mongo_path=<mongo_path>                           Where to save your mongoDB? [Default: ./mongodb].
                                                                If you insert a path that already have a mongoDB in it will include (append)
                                                                the GFF as new collection (<genome_name>) in a new or existing DB (<db_name>).

                                                    Converting to FASTA

        -l, --fasta_features=<feature_types>                    A comma separated list of feature types to extract. Must be identical as written in the GFF file. [Default: CDS]

Example
-------

We can use as an example some of the files stored in the provided `test directory <https://github.com/fmalmeida/gff-toolbox/tree/master/test>`_.

Get mRNA sequences from GFF, in nucleotide fasta
""""""""""""""""""""""""""""""""""""""""""""""""

.. code-block:: bash

    # Example
    ## Get mRNA sequences from GFF, in nucleotide fasta
    gff-toolbox convert -i Athaliana_ref.gff.gz --fasta Athaliana_genome.fasta.gz --format fasta-nt --fasta_features mRNA

    ## Output
    >rna-NM_099983.2
    AAATTATTAGATATACCAAaccagagaaaacaaatacataatCGGAGAAATACAGATTAcagagagcgagagagatcGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGTTTAGTGAAAATGGAGGATCAAGTTGGGTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTTGGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATCAGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTGTAAGTTCCGAATTTTCTGAATTTCATTTGCAAGTAATCGAtttaggtttttgattttagggtttttttttgttttgaacagtCCAGTCAAAGTACAAATCGAGAGATGCTATGTGGTACTTCTTCTCTCgtagagaaaacaacaaaggGAATCGACAGAGCAGGACAACGGTTTCTGGTAAATGGAAGCTTACCGGAGAATCTGTTGAGGTCAAGGACCAGTGGGGATTTTGTAGTGAGGGCTTTCGTGGTAAGATTGGTCATAAAAGGGTTTTGGTGTTCCTCGATGGAAGATACCCtgacaaaaccaaatctgATTGGGTTATCCACGAGTTCCACTACGACCTCTTACCAGAACATCAGGTTTTCttctattcatatatatatatatatatatatgtggatatatatatatgtggttTCTGCTGATTCATAGTTAGAATTTGAGTTATGCAAATTAGAAACTATGTAATGTAACTCTATTTAGGTTCAGCAGCTATTTTAGGCTTAGCTTACTCTCACCAATGTTTTATACTGATGAACTTATGTGCTTACCTCCGGAAATTTTACAGAGGACATATGTCATCTGCAGACTTGAGTACAAGGGTGATGATGCGGACATTCTATCTGCTTATGCAATAGATCCCACTCCCGCTTTTGTCCCCAATATGACTAGTAGTGCAGGTTCTGTGGTGAGTCTTTCTCCATATACACTTAGCTTTGAGTAGGCAGATCAAAAAAGAGCTTGTGTCTACTGATTTGATGTTTTCCTAAACTGTTGATTCGTTTCAGGTCAACCAATCACGTCAACGAAATTCAGGATCTTACAACACTTACTCTGAGTATGATTCAGCAAATCATGGCCAGCAGTTTAATGAAAACTCTAACATTATGCAGCAGCAACCACTTCAAGGATCATTCAACCCTCTCCTTGAGTATGATTTTGCAAATCACGGCGGTCAGTGGCTGAGTGACTATATCGACCTGCAACAGCAAGTTCCTTACTTGGCACCTTATGAAAATGAGTCGGAGATGATTTGGAAGCATGTgattgaagaaaattttgagtttttggtaGATGAAAGGACATCTATGCAACAGCATTACAGTGATCACCGGCCCAAAAAACCTGTGTCTGGGGTTTTGCCTGATGATAGCAGTGATACTGAAACTGGATCAATGGTAAGCTTTTTTTACTCATATATAATCACAACCTATATCGCTTCTATATCTCACACGCTGAATTTTGGCTTTTAACAGATTTTCGAAGACACTTCGAGCTCCACTGATAGTGTTGGTAGTTCAGATGAACCGGGCCATACTCGTATAGATGATATTCCATCATTGAACATTATTGAGCCTTTGCACAATTATAAGGCACAAGAGCAACCAAAGCAGCAGAGCAAAGAAAAGGTTTAACACTCTCACTGAGAAACATGACTTTGATACGAAATCTGAATcaacatttcatcaaaaagaTTTAGTCAAATGACCTCTAAATTATGAGCTATGGGTCTGCTTTCAGGTGATAAGTTCGCAGAAAAGCGAATGCGAGTGGAAAATGGCTGAAGACTCGATCAAGATACCTCCATCCACCAACACGGTGAAGCAGAGCTGGATTGTTTTGGAGAATGCACAGTGGAACTATCTCAAGAACATGATCATTGGTGTCTTGTTGTTCATCTCCGTCATTAGTTGGATCATTCTTGTTGGTTAAGAGGTCAAATCGGATTCTTGCtcaaaatttgtatttcttagaatgtgtgtttttttttgtttttttttctttgctctgttttctcgCTCCGGAAAAGTTTgaagttatattttattagtatgtaaagaagagaaaaagggggaaagaagagagaagaaaaatgcagaaaatcatatatatgaattggaaaaaagtatatgtaataataattagtGCATCGTTTTGTGGTGtagtttatataaataaagtgATATATAGTCTTGTATAAG
    >rna-NM_001331242.1
    ttagtaAGGTCTAATTCAATTTTTGGTGGCGATAATATTTGGCTTAGTCATAAAATACAGTATGGTATAATAATGTAAAGGTTTCTCTTATCTTCAAACCAAAAGACTATACTGGAAGCTGATGGGATCATACGATTCTGAAAAAATAAGACATATATTGCAACAGAGATCCAATTTGTATCAAAAATATTGTCGGCTCAAAAATCTGACCCACCAAGAATCTAATCAAGTGCGCGATTAAGCATACGGCTATGCATCTGGTCATTGTTGATTCAGTCATCACTGGTTTAAAGACAAACTTGCATTGTGAGATTCcaaaataacaacaacaaaaaacaatttgcatTGAGAACATTTTGAAGGTCTGACCTTTAAGAGCCATGGAGTTTGATGTTAAGAGAAGTAtatcgacaaaaaaaatcactgaCATTGGGAATTCCCATACCTGTAATAACAAAGATTCTCTATTTTTGAGCAAAGAGACATAACACCATGTTTAATCAATACAAGAAACTTTAGTGCATGATTTATGGAATGCTTAAGAAGTTTGGAACTTCAATATTAGGAATTAAGTGAGAGGTAAAGCTACAACATACCAACATCGCAAGCAGAAATATCTTGAAGTAACTAGAGATGAATATCCCCAAcataatctctcttcttctgcaaagtttttaaaaaaaattatacataaaCAATCTTCAGGTGACAGAAGTCTGAGATCTTTGATGAAAAACTCATAATAAgaactaagaagaagaagaatcagtaATCACCTGGAAACTTCATTTAGCAAACCCTTAGTCGCAATGGCAAAAGAGATGATAAATGCAGCGTTTGCAGATAAGACACCAATCAGAACCTGTTTCAAATGCGAAATTATTACCCTTTCTAAACAATCTCAATGACTTAAATCATTTAAAccttaaaggaaaaaaaatctaattaagtCCATTAAAAAGAAACGATCTAACCTTTATAGATAGAAGAACAGGGCTATCAGAAAAGCTCGATTCTTCATCACTTTTTCTCAGTAGCAAGCTTCTATCTATGATGAGTTCAGGGCTTTTCAAATAAGTTTGCCGATGAATCTCACCAACTACACATCTGCTAGCTACACTTTGATAGTAAAagattataaaacaaaaggataCAACAGTCTAGAAGAAGATAGGCGAAGACCAACTTCCACAACAGATGctgcacacacacacacaaaaaaaaaaaaagaacccaACAATTCTTATTGGATCAGAGACTACTCAATATCCCCAAACTTGGAAATTAGTTTGTTGCTTGAGGTCTAAGATACTTCTATATATGgaaaaagattttcaaagCCAGATATTTCCACAAGTTTGTAATATCAATTCAAGATAAGAGAGCTAGAATCAGACAGGAACTAGCAATGCTTGAAATCAAGAACTTGAATTGAAATAGTTTTTTACCTGAATATTGACAGTTGCTGGATTAATTGCATTGTAGAGGACGTGTCTATATACCTTTGGTCTGTGAAGGATTAAATCGATGAAAATAATCTgccaaagaaaacaattaaagaaccaaaaaccaaaattggaAAGAAATAGGGAAACACccaaaaagggaaagaaagtGATTAAAACAGACCATGCGTTCACACTCGATGTACTCATCTGCTACTTCCTTGCAATTTCCCTAAATATAACAATATGATCAAAGATGGAAACTTTGAAGAAATTTAATAGAGAATcttataaaccctaattggGTCAAAGAAGATCcattaatacaaaaatcttACGCATTTCATGAGACGAATGTTACCCGGAGAGTATTGAATGAACAATGACTTTACCCTAAAACCACATCCCACGCATCTGTGTTCACTCGCCGCcattgctctctctctctctctctctctctctctctcaagagaagaagaatacgGAGCAATTAGAGTCCGGGTCTGGGCTACTGTTTTAACCCTAAATGGGCTTATTCATGGGCCAAGTTTTTGAAGTCTtaactttaaatttgttaGGCCCACTTTTGCTCTAAGCCGGGGTATTTGTAccccaaaatttaaaaatcatatacacgttgtaatttataaatagttCAATTTGGATCAAAATCTTGTCCATATGACATAGCATTTTAAAATGCGTAGGTTCATGAATGAAACATATTATAGGCCTCAGATAAAGATATACATATTAAGTCTAAATTATTTAGTCTTCAGAATTTACCACACTTACTGAAAAGTCTAGTGGTTCACTAATATTATTACTGTCGTGTTACTTTCTATATATAGTTCATGACTTGTGAGTTGTGATGGATaagtttataagaaaataaatt
    ...

Converting a GFF to JSON, without grabbing sequences from genome
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. code-block:: bash

    # Example
    ## Converting a GFF to JSON, without grabbing sequences from genome
    gff-toolbox convert -i Kp_ref.gff -f json

    # Output
    {
        ...,
        {
            "recid": "NC_016845.1",
            "source": "RefSeq",
            "type": "gene",
            "start": "382",
            "end": "822",
            "score": ".",
            "strand": "-",
            "phase": ".",
            "attributes": {
                "ID": "gene-KPHS_00010",
                "Dbxref": "GeneID:11849782",
                "Name": "KPHS_00010",
                "gbkey": "Gene",
                "gene_biotype": "protein_coding",
                "locus_tag": "KPHS_00010"
            }
        },
        ...
    }

Converting a GFF to JSON, grabbing sequences from fasta
"""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. code-block:: bash

    # Example
    ## Converting a GFF to JSON, grabbing sequences from fasta
    gff-toolbox convert -i Kp_ref.gff -f json --fasta Kpneumoniae_genome.fasta

    # Output
    {
        ...,
        {
            "recid": "NC_016845.1",
            "source": "RefSeq",
            "type": "gene",
            "start": "382",
            "end": "822",
            "score": ".",
            "strand": "-",
            "phase": ".",
            "attributes": {
                "ID": "gene-KPHS_00010",
                "Dbxref": "GeneID:11849782",
                "Name": "KPHS_00010",
                "gbkey": "Gene",
                "gene_biotype": "protein_coding",
                "locus_tag": "KPHS_00010"
            },
            "sequence": "TCAGAGCAAATTTTTCCAGGATCCCAGCCAGATTTCCGCTGGATCCTCCGGAATTTCGTGTTCAAGGACGTTGATCTTCAGCGTCTCACCGATCTGTTTTGCCCCGCAGGCCTTCAGCGCCGCTTCTAATTTATCGATCGCGCCGCAAAAGGTGTCGTATTCGCGGCTGCCGATCCCTATAGCGCCAAAGCGCACGGCGGAAAGATCGGCCTGCTTCTCCTGCAGCGCGTCATAGAAGGGAACCAGGTTGTCGGGAACGTCGCCGGCGCCATGGGTGGAAGAGATAATTAACCAGATCCCTTGCGCCTGCAGATCGTCCACCAGCGGACCGTGCTGGATTTCAGTGCTATAGCCGGCTTCTTCCAGCTTTTCCGCCAGATGCTCTGCTACATATTCCGCACTGCCCAGGGTGCTGCCGCTGATAAGTGTAATGTCTGCCAT"
        },
        ...
    }





