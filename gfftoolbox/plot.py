## Def help message
usage_plot="""
gff-toolbox:

            Plot

Enabling the visualization of a genomic region from a GFF using the DNA features python package

usage:
    gff-toolbox plot -h|--help
    gff-toolbox plot check-gff ( --input <gff> )
    gff-toolbox plot ( --input <gff> | --fofn <file> ) ( --contig <contig_name> ) [ --start <start_base> --end <end_base> --feature <feature_type> --identification <id> --title <title> --label <label> --color <color> --output <png_out> --width <width> --height <height> ]

options:
    -h --help                               Show this screen.

    -v --version                            Show version information.

    check-gff                               Does a simple parsing of the GFF file so the user knows the available qualifiers that
                                            can be used as gene identifiers. GFF qualifiers are retrieved from the 9th column.
                                            Same as gff-toolbox overview command.

    -i, --input=<gff>                       Used to plot dna features from a single GFF file.

    -f, --fofn=<file>                       Used to plot dna multiple features from multiple GFF files. Contents must be in csv format with 3 columns:
                                            gff,custom_label,color (HEX format). Features from each GFF will have the color set in the 3rd column.
                                            Useful to compare annotations and to plot features with different colors if you separate them into multiple
                                            gff files each containing one type of feature.

    --start=<start_base>                    Starting position for plotting. [Default: 1].

    --end=<end_base>                        Ending position for plotting. [Default: 500].

    --contig=<contig_name>                  Name of the contig which you want to plot.

    --identification=<id>                   Which GFF qualifier must be used as gene identification to write in the plot?
                                            Please check for available qualifiers with 'check-gff'. [Default: ID].

    -t, --title=<title>                     Legend/plot title. [Default: Gene Plot].

    -l, --label=<label>                     Custom label for plotting features. This is the string that appears in the legend. [Default: Gene].

    --feature=<feature_type>                Type of the GFF feature (3rd column) which you want to plot. It is possible to set more than one feature to be
                                            plotted by giving it separated by commas, eg. CDS,rRNA. [Default: gene].

    --color=<color>                         HEX entry for desired plotting color. [Default: #ccccff].

    --width=<width>                         Plot width ratio. [Default: 20].

    --height=<height>                       Plot height ratio. [Default: 5].

    -o, --output=<png_out>                  Output PNG filename. [Default: ./out.png].

example:

    ## Plotting CDS and/or rRNA features found inside the region between base 1 and 1000 from
    ## contig_1_segment0 sequence. Without much customization. Giving a custom label for the
    ## genes to appear in the legend and giving a different legend title.

$ gff-toolbox plot -i test/input.gff --contig contig_1_segment0 --feature CDS,rRNA --start 1 --end 10000 -l "Generic features (CDS and rRNAs)" -t "Kp annotation"
"""

##################################
### Loading Necessary Packages ###
##################################
from dna_features_viewer import *
import Bio.SeqIO
from pprintpp import pprint
from BCBio import GFF
from BCBio.GFF import GFFExaminer
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


##################################################
### Function for checking available qualifiers ###
##################################################
def check_gff(infile):

    # GFF overview
    print("GFF overview:\n")
    examiner = GFFExaminer()
    in_handle = open(infile)
    pprint(examiner.available_limits(in_handle))
    in_handle.close()
    print("")

    # Load GFF and its sequences
    gff = GFF.parse(infile)

    # Check qualifiers
    for rec in gff:
        print("Example of the GFF's first line available qualifiers from the 9th column:\n")
        print(rec.features[0])
        print("\nPlease select only one of the available qualifiers to be used as gene identification!")
        exit()

######################################################
### Function for execution with a single GFF input ###
######################################################

def single_gff(infile, start, end, contig, feature, qualifier, coloring, custom_label, outfile, plot_title, plot_width, plot_height):

    # Subset GFF based on chr and feature type
    limit_info = dict(
            gff_id   = [contig],
            gff_type = list(feature.split(','))
    )

    # Load GFF and its sequences
    gff = GFF.parse(infile, limit_info=limit_info)

    # Create empty features and legened list
    features = []

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
                input= GraphicFeature(start=int(rec.features[i].location.start), end=int(rec.features[i].location.end),
                                      strand=int(strand), label=str(rec.features[i].qualifiers[qualifier][0]), color=coloring)

                # Append
                features.append(input)


    # Draw plot
    record = GraphicRecord(sequence_length=length, features=features, first_index=start_nt)
    ax, _ = record.plot(figure_width=int(plot_width), figure_height=int(plot_height))
    ## Label not in the gene plot (using separate legend box)
    legend = ax.legend(handles=[mpatches.Patch(facecolor=coloring, label="{0}".format(custom_label), linewidth = 0.5, edgecolor = 'black')],
                       loc = 1, title=plot_title, fontsize = 'medium', fancybox = True)
    ax.figure.savefig(outfile, bbox_inches='tight')



#######################################################
### Function for execution with multiple GFF inputs ###
#######################################################

def multiple_gff(input_fofn, start, end, contig, qualifier, feature, outfile, plot_title, plot_width, plot_height):

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
        limit_info = dict(
                gff_id   = [contig],
                gff_type = list(feature.split(','))
        )

        # Load GFF and its sequences
        gff = GFF.parse(infile, limit_info=limit_info)

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
                    input= GraphicFeature(start=int(rec.features[i].location.start), end=int(rec.features[i].location.end),
                                          strand=int(strand), label=str(rec.features[i].qualifiers[qualifier][0]), color=coloring)

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
