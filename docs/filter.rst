.. _filter:

GFF-toolbox filter
==================

This command uses Biopython library to parse and filter your GFF file as you wish. It targets the attributes column of the GFF. Many options are possible since it has the options to search for patterns in a grep-like manner where users must specify the pattern and the column to search it or, it can use biopython dictionary structure **(default)** that is more recommended for nested GFFs.

Help
----

.. code-block:: bash

    # Trigger help
    gff-toolbox filter -h

    # Help
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

        -i, --input=<gff>                                       Input GFF file. GFF file must not contain sequences, only features [Default: stdin].

        -m, --mode=<search_mode>                                In which mode to search for patterns: loose or exact?
                                                                The loose mode, scans the GFF in a grep-like manner via pandas dataframes in which the user must specify
                                                                a pattern and a column to search it. Recommended for simple searches were nest structure is not a must.
                                                                The exact mode scans the GFF with Biopython and BCBio packages, treating it as python dictionary. It is
                                                                recommended for more complex searches and complex GFFs, such as nested GFFs. [Default: exact]

        --strand=<strand>                                       Apply a filter based on the strand of the feature. Options: plus or minus. By default, everything is given.
                                                                In exact mode, this filter is applied in the parent feature, if it passes, it\'s children are also printed.
                                                                The contrary is also true. In the loose mode it is applied directly to all features, nested or not.

        --start=<start_position>                                Apply a filter to select features starting from this position. In exact mode, this filter is applied in the
                                                                parent feature, if it passes, it\'s children are also printed. The contrary is also true. In the loose mode
                                                                it is applied directly to all features, nested or not.

        --end=<end_position>                                    Apply a filter to select features until this position. In exact mode, this filter is applied in the parent
                                                                feature, if it passes, it\'s children are also printed. The contrary is also true.

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

Differing loose and exact filter modes
--------------------------------------

Oi

Example
-------

We can use as an example some of the files stored in the provided `test directory <https://github.com/fmalmeida/gff-toolbox/tree/master/test>`_.

Filter per contig and pattern in attributes
""""""""""""""""""""""""""""""""""""""""""""

.. code-block:: bash

    # Example
    ## All the CDS(s) found in the sequence named NC_016845.1 that have the word "transcriptional regulator" in their attributes.
    gff-toolbox filter -i Kp_ref.gff --chr NC_016845.1 --type CDS | gff-toolbox filter --mode loose -p "transcriptional regulator"

    # Output
    NC_016845.1	RefSeq	CDS	922	1380	.	-	0	Dbxref=Genbank:YP_005224302.1,GeneID:11849790;ID=cds-YP_005224302.1;Name=YP_005224302.1;Parent=gene-KPHS_00020;gbkey=CDS;locus_tag=KPHS_00020;product=DNA-binding transcriptional regulator AsnC;protein_id=YP_005224302.1;transl_table=11
    NC_016845.1	RefSeq	CDS	15004	15705	.	-	0	Dbxref=Genbank:YP_005224314.1,GeneID:11844989;ID=cds-YP_005224314.1;Name=YP_005224314.1;Parent=gene-KPHS_00140;gbkey=CDS;locus_tag=KPHS_00140;product=putative transcriptional regulator;protein_id=YP_005224314.1;transl_table=11
    NC_016845.1	RefSeq	CDS	43194	43577	.	+	0	Dbxref=Genbank:YP_005224336.1,GeneID:11845014;ID=cds-YP_005224336.1;Name=YP_005224336.1;Parent=gene-KPHS_00360;gbkey=CDS;locus_tag=KPHS_00360;product=putative 2-component transcriptional regulator;protein_id=YP_005224336.1;transl_table=11
    NC_016845.1	RefSeq	CDS	78374	79072	.	-	0	Dbxref=Genbank:YP_005224372.1,GeneID:11845050;ID=cds-YP_005224372.1;Name=YP_005224372.1;Parent=gene-KPHS_00720;gbkey=CDS;locus_tag=KPHS_00720;product=DNA-binding transcriptional regulator CpxR;protein_id=YP_005224372.1;transl_table=11
    ...

Filter per start position and custom attributes
""""""""""""""""""""""""""""""""""""""""""""""""

.. code-block:: bash

    ## Filtering a set of genes and its childs using a file containing the desired attributes.
    ## A. thaliana annotation. Also give a custom start position for features to be printed?
    gff-toolbox filter -i Athaliana_ref.gff.gz --attributes atts2.txt --start 5900

    # Output
    ##gff-version 3
    ##sequence-region NC_003070.9 1 30427671
    NC_003070.9	annotation	remark	1	30427671	.	.	.	gff-version=3;sequence-region=%28%27NC_003070.9%27%2C 0%2C 30427671%29,%28%27NC_003071.7%27%2C 0%2C 19698289%29,%28%27NC_003074.8%27%2C 0%2C 23459830%29,%28%27NC_003075.7%27%2C 0%2C 18585056%29,%28%27NC_003076.8%27%2C 0%2C 26975502%29,%28%27NC_037304.1%27%2C 0%2C 367808%29,%28%27NC_000932.1%27%2C 0%2C 154478%29;species=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi%3Fid%3D3702
    NC_003070.9	RefSeq	gene	6788	9130	.	-	.	Dbxref=Araport:AT1G01020,TAIR:AT1G01020,GeneID:839569;ID=gene-AT1G01020;Name=ARV1;gbkey=Gene;gene=ARV1;gene_biotype=protein_coding;gene_synonym=T25K16.2,T25K16_2;locus_tag=AT1G01020
    NC_003070.9	RefSeq	mRNA	6788	9130	.	-	.	Dbxref=Araport:AT1G01020,GeneID:839569,Genbank:NM_001331242.1,TAIR:AT1G01020;ID=rna-NM_001331242.1;Name=NM_001331242.1;Parent=gene-AT1G01020;gbkey=mRNA;gene=ARV1;locus_tag=AT1G01020;orig_protein_id=gnl%7CJCVI%7CAT1G01020.5;orig_transcript_id=gnl%7CJCVI%7CmRNA.AT1G01020.5;product=ARV1 family protein;transcript_id=NM_001331242.1
    NC_003070.9	RefSeq	exon	8594	9130	.	-	.	Dbxref=Araport:AT1G01020,GeneID:839569,Genbank:NM_001331242.1,TAIR:AT1G01020;ID=exon-NM_001331242.1-1;Parent=rna-NM_001331242.1;gbkey=mRNA;gene=ARV1;locus_tag=AT1G01020;orig_protein_id=gnl%7CJCVI%7CAT1G01020.5;orig_transcript_id=gnl%7CJCVI%7CmRNA.AT1G01020.5;product=ARV1 family protein;transcript_id=NM_001331242.1
    [...]

Filter per pattern in grep-like manner
""""""""""""""""""""""""""""""""""""""

.. code-block:: bash

    ## Simple filter in any column: wheter a line contain a pattern in a specific column (like grep)
    ## Check the features that have the word "putative" in their attributes.
    gff-toolbox filter --mode loose --sort --header -i Kp_ref.gff -p "putative"

    # Output
    ##gff-version 3
    NC_016838.1	RefSeq	CDS	3127	4149	.	-	0	ID=cds-YP_005220812.1;Parent=gene-KPHS_p100050;Dbxref=Genbank:YP_005220812.1,GeneID:11818042;Name=YP_005220812.1;gbkey=CDS;locus_tag=KPHS_p100050;product=putative recombinase;protein_id=YP_005220812.1;transl_table=11
    NC_016838.1	RefSeq	CDS	5479	6486	.	-	0	ID=cds-YP_005220815.1;Parent=gene-KPHS_p100080;Dbxref=Genbank:YP_005220815.1,GeneID:11818045;Name=YP_005220815.1;gbkey=CDS;locus_tag=KPHS_p100080;product=putative regulator;protein_id=YP_005220815.1;transl_table=11
    NC_016838.1	RefSeq	CDS	9046	12090	.	-	0	ID=cds-YP_005220821.1;Parent=gene-KPHS_p100140;Dbxref=Genbank:YP_005220821.1,GeneID:11817910;Name=YP_005220821.1;gbkey=CDS;locus_tag=KPHS_p100140;product=putative DNA polymerase III alpha subunit;protein_id=YP_005220821.1;transl_table=11
    [...]

A more complex example
""""""""""""""""""""""

.. code-block:: bash

    ## In the example below, we filter the GFF in a more complex manner:
    ## All the CDS(s) found in the sequence named NC_016845.1 that
    ## have the word "transcriptional regulator" in their attributes.
    ##
    ## It works in both ways:
    gff-toolbox filter -i Kp_ref.gff --chr NC_016845.1 --type CDS | gff-toolbox filter --mode loose -p "transcriptional regulator"
    gff-toolbox filter --mode loose -i Kp_ref.gff -p "transcriptional regulator" | gff-toolbox filter --chr NC_016845.1 --type CDS

    # Output
    NC_016845.1	RefSeq	CDS	922	1380	.	-	0	Dbxref=Genbank:YP_005224302.1,GeneID:11849790;ID=cds-YP_005224302.1;Name=YP_005224302.1;Parent=gene-KPHS_00020;gbkey=CDS;locus_tag=KPHS_00020;product=DNA-binding transcriptional regulator AsnC;protein_id=YP_005224302.1;transl_table=11
    NC_016845.1	RefSeq	CDS	15004	15705	.	-	0	Dbxref=Genbank:YP_005224314.1,GeneID:11844989;ID=cds-YP_005224314.1;Name=YP_005224314.1;Parent=gene-KPHS_00140;gbkey=CDS;locus_tag=KPHS_00140;product=putative transcriptional regulator;protein_id=YP_005224314.1;transl_table=11
    NC_016845.1	RefSeq	CDS	43194	43577	.	+	0	Dbxref=Genbank:YP_005224336.1,GeneID:11845014;ID=cds-YP_005224336.1;Name=YP_005224336.1;Parent=gene-KPHS_00360;gbkey=CDS;locus_tag=KPHS_00360;product=putative 2-component transcriptional regulator;protein_id=YP_005224336.1;transl_table=11
    NC_016845.1	RefSeq	CDS	78374	79072	.	-	0	Dbxref=Genbank:YP_005224372.1,GeneID:11845050;ID=cds-YP_005224372.1;Name=YP_005224372.1;Parent=gene-KPHS_00720;gbkey=CDS;locus_tag=KPHS_00720;product=DNA-binding transcriptional regulator CpxR;protein_id=YP_005224372.1;transl_table=11
    [...]