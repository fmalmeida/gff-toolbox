## Def help message
usage_plot="""
gff-toolbox:

            Plot

Enabling the visualization of a genomic region from a GFF using the DNA features python package

usage:
    gff-toolbox plot -h|--help
    gff-toolbox plot check-gff [ --input <gff> ]
    gff-toolbox plot --plot features ( --contig <contig_name> ) [ --input <gff> | --fofn <file> ] [ --start <start_base> --end <end_base> --feature <feature_type> --identification <id> --title <title> --label <label> --color <color> --output <png_out> --width <width> --height <height> ]
    gff-toolbox plot --plot ideogram ( --ref_fasta <ref_fasta> ) [ --input <gff> --feature <feature_type> --chr_minsize <int> --chr_maxsize <int> --title <title> --width <width> --height <height> ]

options:

                                ###############
                                ### General ###
                                ###############

    -h --help                               Show this screen.

    -v --version                            Show version information.

    check-gff                               Does a simple parsing of the GFF file so the user knows the 
                                            available qualifiers that can be used as gene identifiers.
                                            GFF qualifiers are retrieved from the 9th column.
                                            Same as gff-toolbox overview command.
                                
                                ##############
                                ### Inputs ###
                                ##############

    -i, --input=<gff>                       Used to plot dna features from a single GFF file [Default: stdin].

    -f, --fofn=<file>                       Used to plot dna multiple features from multiple GFF files. 
                                            Contents must be in csv format with 3 columns:
                                            gff,custom_label,color (HEX format). Features from each GFF will have the color set in the 3rd column.
                                            Useful to compare annotations and to plot features with different colors if you separate them into multiple
                                            gff files each containing one type of feature.

                                ################################
                                ### Definitions for any plot ###
                                ################################
                            
    -p, --plot=<plot_type>                  Which plot type to draw? Options: features;ideogram. [Default: features].

    -t, --title=<title>                     Legend/plot title. [Default: Gene Plot].

    -l, --label=<label>                     Custom label for plotting features. This is the string that 
                                            appears in the legend. [Default: Gene].

    --width=<width>                         Plot width ratio. [Default: 20].

    --height=<height>                       Plot height ratio. [Default: 5].

    -o, --output=<png_out>                  Output SVG/PNG filename. [Default: ./out.png].
                                            You can output SVG with: "-o out.svg"
    
    --color=<color>                         HEX entry for desired plotting color. [Default: #ccccff].
    
    --feature=<feature_type>                Type of the GFF feature (3rd column) which you want to plot.
                                            It is possible to set more than one feature to be
                                            plotted by giving it separated by commas, eg. CDS,rRNA.

                                #####################################
                                ### Definitions for features plot ###
                                #####################################

    --start=<start_base>                    Starting position for plotting. [Default: 1].

    --end=<end_base>                        Ending position for plotting. [Default: 500].

    --contig=<contig_name>                  Name of the contig that you want to plot. Required for the features plot.

    --identification=<id>                   Which GFF qualifier must be used as gene identification to write in the plot?
                                            Please check for available qualifiers with 'check-gff'. Must be the exact name
                                            of the key found in the attributes columns. [Default: ID].

                                #####################################
                                ### Definitions for ideogram plot ###
                                #####################################
    
    --ref_fasta=<ref_fasta>                 FASTA file to be used as reference when plotting the karyotypes. Sequence chrs in FASTA and 
                                            the GFF must be related (same name).
    
    --chr_minsize=<int>                     Min. size of contigs to plot. [Default: 1]

    --chr_maxsize=<int>                     Max. size of contigs to plot.

example:

    ## Plotting CDS and/or rRNA features found inside the region between base 1 and 10000 from
    ## contig_1_segment0 sequence. Without much customization. Giving a custom label for the
    ## genes to appear in the legend and giving a different legend title.

$ gff-toolbox plot -i Kp_ref.gff --contig NC_016845.1 --feature CDS,rRNA --start 10000 --end 20000 -l "Generic features (CDS and rRNAs)" -t "Kp annotation"

    ## Plotting CDS, rRNA and tRNA features with different colors. Setting the genomic region
    ## in NC_016845.1:10000-20000. Checkout the example fofn file to understand it better.

$ gff-toolbox plot -f kp_gffs.fofn --start 10000 --end 20000 --contig NC_016845.1 -t "Kp annotation" --feature CDS,rRNA,tRNA

    ## Same as above.
    ## This time, instead of plotting gene names we plot the gene products by setting the
    ## parameter --identification to the exact name of the key in the attributes column.

$ gff-toolbox plot -f kp_gffs.fofn --start 10000 --end 20000 --contig NC_016845.1 -t "Kp annotation" --feature CDS,rRNA,tRNA --identification product

    ## Generating an ideogram plot using karyoploteR package
    ## plotting only the rRNA features

$ gff-toolbox plot -p ideogram --ref_fasta Kpneumoniae_genome.fasta --input Kp_ref.gff --feature rRNA
"""

##################################
### Loading Necessary Packages ###
##################################
from dna_features_viewer import *
from pprintpp import pprint
from BCBio import GFF
from BCBio.GFF import GFFExaminer
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import binascii
import gzip
import tempfile
import sys
from pyfaidx import Fasta
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.collections import BrokenBarHCollection

###################
### Stdin Check ###
###################
def stdin_checker(input, sys_contents):
    # Checking for stdin
    if input == "stdin":
        tmp = tempfile.NamedTemporaryFile(mode = "w+t", delete = False) # Create tmp file to work as input
        temp_file = open(tmp.name, 'w')
        for line in sys_contents:
            temp_file.writelines(f"{line}")

        temp_file.seek(0)
        input = tmp.name

    else:
        pass

    return input

###################
### Gzip opener ###
###################
def gzip_opener(input, mode_in):
    if  binascii.hexlify(open(input, 'rb').read(2)) == b"1f8b" or input.endswith(".gz"):
        return gzip.open(input, mode=mode_in)
    else:
        return open(input, mode=mode_in)

##################################################
### Function for checking available qualifiers ###
##################################################
def check_gff(infile):

    # Stdin check
    sys_contents = ""
    if infile == "stdin":
        sys_contents = sys.stdin.readlines()
    else:
        pass

    # GFF overview
    print("GFF overview:\n")
    examiner = GFFExaminer()
    pprint(examiner.available_limits(gzip_opener(stdin_checker(infile, sys_contents), "rt")))
    print("")

    # Check qualifiers
    for rec in GFF.parse(gzip_opener(stdin_checker(infile, sys_contents), "rt")):
        print("Example of the GFF's first line available qualifiers from the 9th column:\n")
        print(rec.features[0])
        print("\nPlease select only one of the available qualifiers to be used as gene identification!")
        exit()

####################################
### Functions for features plots ###
####################################

# Single GFF
def features_single_gff(infile, start, end, contig, feature, qualifier, coloring, custom_label, outfile, plot_title, plot_width, plot_height):

    # Stdin check
    sys_contents = ""
    if infile == "stdin":
        sys_contents = sys.stdin.readlines()
    else:
        pass

    # Subset GFF based on chr and feature type
    limit_info = dict()
    limit_info["gff_id"] = [contig]
    if feature != None:
        limit_info['gff_type'] = list(feature.split(','))

    # Create empty features and legened list
    features = []

    ## Populate features list
    ## Filtering by location
    start_nt = int(start)
    end_nt   = int(end)
    length   = end_nt - start_nt

    for rec in GFF.parse(gzip_opener(stdin_checker(infile, sys_contents), "rt"), limit_info=limit_info):
        for i in range(0, len(rec.features)):
            if ( int(rec.features[i].location.start) >= int(start_nt) and int(rec.features[i].location.end) <= int(end_nt) ):

                if (str(rec.features[i].location.strand) == "+"):
                    strand=+1
                else:
                    strand=-1

                ## Label not in the gene plot
                    if (qualifier in rec.features[i].qualifiers):
                        if (rec.features[i].qualifiers[qualifier][0] == "true"):
                            input= GraphicFeature(start=int(rec.features[i].location.start), end=int(rec.features[i].location.end),
                                                  strand=int(strand), color=coloring)
                        else: 
                            input = GraphicFeature(start=int(rec.features[i].location.start), end=int(rec.features[i].location.end),
                                            strand=int(strand), label=str(rec.features[i].qualifiers[qualifier][0]), color=coloring)
                    else:
                        input= GraphicFeature(start=int(rec.features[i].location.start), end=int(rec.features[i].location.end),
                                          strand=int(strand), color=coloring)

                # Append
                features.append(input)


    # Draw plot
    record = GraphicRecord(sequence_length=length, features=features, first_index=start_nt)
    ax, _ = record.plot(figure_width=int(plot_width), figure_height=int(plot_height))
    ## Label not in the gene plot (using separate legend box)
    legend = ax.legend(handles=[mpatches.Patch(facecolor=coloring, label="{0}".format(custom_label), linewidth = 0.5, edgecolor = 'black')],
                       loc = 1, title=plot_title, fontsize = 'medium', fancybox = True)
    ax.figure.savefig(outfile, bbox_inches='tight')

# Multiple GFFs
def features_multiple_gff(input_fofn, start, end, contig, qualifier, feature, outfile, plot_title, plot_width, plot_height):

    # Open list of filenames containing GFFs
    file = open(input_fofn, 'r')
    content = file.readlines()

    # Create empty features and legened list
    features = []
    legend_entries = []

    # Begin Parsing
    for line in content:
        data = line.strip().split(",", 3)
        infile    = data[0]
        labeling  = data[1]
        coloring  = data[2]

        # Subset GFF based on chr and feature type
        limit_info = dict()
        limit_info["gff_id"] = [contig]
        if feature != None:
            limit_info['gff_type'] = list(feature.split(','))

        # Load GFF and its sequences
        gff = GFF.parse(gzip_opener(infile, "rt"), limit_info=limit_info)

        ## Populate features list
        ## Filtering by location
        start_nt = int(start)
        end_nt   = int(end)
        length   = end_nt - start_nt

        for rec in gff:
            for i in range(0, len(rec.features)):
                if ( int(rec.features[i].location.start) >= int(start_nt) and int(rec.features[i].location.end) <= int(end_nt) ):

                    if (str(rec.features[i].location.strand) == "+"):
                        strand=+1
                    else:
                        strand=-1

                    ## Label not in the gene plot
                    if (qualifier in rec.features[i].qualifiers):
                        if (rec.features[i].qualifiers[qualifier][0] == "true"):
                            input= GraphicFeature(start=int(rec.features[i].location.start), end=int(rec.features[i].location.end),
                                                  strand=int(strand), color=coloring)
                        else: 
                            input = GraphicFeature(start=int(rec.features[i].location.start), end=int(rec.features[i].location.end),
                                            strand=int(strand), label=str(rec.features[i].qualifiers[qualifier][0]), color=coloring)
                    else:
                        input= GraphicFeature(start=int(rec.features[i].location.start), end=int(rec.features[i].location.end),
                                          strand=int(strand), color=coloring)

                    # Append DNA features plot
                    features.append(input)

        # Append to legend
        legend_entries.append(
            mpatches.Patch(facecolor=coloring, label="{0}".format(labeling), linewidth = 0.5, edgecolor = 'black')
        )


    # Draw plot
    record = GraphicRecord(sequence_length=length, features=features, first_index=start_nt)
    ax, _ = record.plot(figure_width=int(plot_width), figure_height=int(plot_height))
    ## Label not in the gene plot (using separate legend box)
    legend = ax.legend(handles=legend_entries, loc = 1, title=plot_title, fontsize = 'medium', fancybox = True)
    ax.figure.savefig(outfile, bbox_inches='tight')

####################################
### Functions for ideogram plots ###
####################################
def generate_bed(ref_fasta):
    rows = []
    for rec in Fasta(ref_fasta):
        rows.append([ rec.name, 1, len(rec), "whitesmoke" ])
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "colors"])

# the following functions are heavily inspired by the follwoing gist:
# https://gist.github.com/daler/c98fc410282d7570efc3
def chromosome_collections(df, y_positions, height, **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'])
    if del_width:
        del df['width']

def ideogram_plot(ref_fasta):

    # call chromsizes
    chromsizes = generate_bed(ref_fasta)

    # initialize a figure/
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    # add regions to the axis
    chrom_height = 1  # Height of each ideogram
    chrom_spacing = 1  # Spacing between consecutive ideograms

    gene_height = 0.4  # Height of the gene track < `chrom_spacing`
    gene_padding = 0.1  # Padding between the top of a gene track ideogram

    # Decide which chromosomes to use
    chromosome_list = list(chromsizes.chrom)

    # Keep track of the y positions genes for each chromosome
    # and center of each ideogram (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase = {}
    gene_ybase = {}
    chrom_centers = {}

    # Items in the beginning of `chromosome_list` will
    # appear at the top of the plot
    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.
        gene_ybase[chrom] = ybase - gene_height - gene_padding
        ybase += chrom_height + chrom_spacing

    # add chromsizes to ax
    ideo = chromsizes

    # Add a new column for width
    ideo['width'] = ideo.end - ideo.start

    # Add a collection of the query chromsizes
    for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, edgecolor='k'):
        ax.add_collection(collection)
    
    # legend
    leg = []
    leg_lab = []

    # add to ax
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)

    ax.legend(leg, leg_lab, loc=4)
    ax.axis("tight")

    # xtick formatting: position in Mbps
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, n: int(x/1e6)))

    # labels
    ax.set_title(f"Karyopype", fontsize=14)
    plt.xlabel('Chromosome Position (Mbp)', fontsize=14)
    plt.ylabel(f'Chromosome', fontsize=14)
    plt.savefig(f'_karyopype.png')

    return(plt)