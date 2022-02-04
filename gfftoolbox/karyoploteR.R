#!/usr/bin/env Rscript

################
### Template ###
################
# Felipe Almeida <marques.felipe@aluno.unb.br>

#################
### Tutorials ###
#################
# https://bernatgel.github.io/karyoploter_tutorial/#Tutorial
# https://www.bioconductor.org/packages/release/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html
# https://www.rapidtables.com/web/color/RGB_Color.html

#####################
### Help function ###
#####################
suppressMessages(library(docopt))
'usage: karyoploteR.R [ -h|--help ] ( --ref=<file> --features=<file> --yaml=<file>)
options:
  -h, --help                      Show help.
  -r, --ref=<file>                BED file to be used as input for karyotypes.
  -f, --features=<file>           BED file containing the features to be plotted.
  -y, --yaml=<file>               YAML file containing plot configuration' -> doc
opt <- docopt(doc)

######################
### Load libraries ###
######################
print(" > Loading libraries!")
suppressMessages(library(karyoploteR))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))
suppressMessages(library(plyranges))
suppressMessages(library(dplyr))
suppressMessages(library(yaml))

#############################
### Load params from YAML ###
#############################
print(" > Loading YAML file!")
config = yaml.load_file(opt$yaml)

#######################################
### Function to load BED as GRanges ###
#######################################
grFromBed <- function(bed, minSize=1, maxSize=NULL) {
  
  if (!is.null(maxSize)) {
    # Open bed
    genome.bed <- read.csv(sep = "\t", header = FALSE,
                           stringsAsFactors = FALSE,
                           file = bed,
                           col.names = c('Chr', 'Start', 'End')) %>%
      dplyr::filter((End - Start) >= minSize & (End - Start) <= maxSize)
    
    # Set BED start to one if zero
    genome.bed[, 1][genome.bed[, 1] == 0] <- 1
    
    # Return GRanges
    return(
      with(genome.bed, GRanges(Chr, IRanges(Start, End)))
    )
    
  } else {
    
    # Open bed
    genome.bed <- read.csv(sep = "\t", header = FALSE,
                           stringsAsFactors = FALSE,
                           file = bed,
                           col.names = c('Chr', 'Start', 'End')) %>%
      dplyr::filter((End - Start) >= minSize)
    
    # Set BED start to one if zero
    genome.bed[, 1][genome.bed[, 1] == 0] <- 1
    
    # Return GRanges
    return(
      with(genome.bed, GRanges(Chr, IRanges(Start, End)))
    )
    
  }
}

####################
### Loading data ###
####################
print(" > Loading BED files!")
# Import karyotypes BED as GRanges
if (config$chr_maxsize != "ALL") {
  genome.gr <- grFromBed(bed=opt$ref, minSize = as.integer(config$chr_minsize), maxSize = as.integer(config$chr_maxsize))
} else {
  genome.gr <- grFromBed(bed=opt$ref, minSize = as.integer(config$chr_minsize))
}

# Import features as GRanges
if (config$feat_maxsize != "ALL") {
  features.gr <- grFromBed(bed=opt$features, minSize = as.integer(config$feat_minsize), maxSize = as.integer(config$feat_maxsize))
} else {
  features.gr <- grFromBed(bed=opt$features, minSize = as.integer(config$feat_minsize))
}

#####################################
### Definition of plot parameters ###
#####################################
print(" > Plotting!")
## Initiate an instance for the plot
svg(config$output_filename, width = as.integer(config$width), height = as.integer(config$height))

## In this area, we use the function getDefaultPlotParams so we have a 
## small list which contains all the necessary parameters used
## when plotting a karyoploteR ideogram. Whith this list we can
## change some parameters and customize the plot a little bit
##
## Read: https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotParams/PlotParams.html
pp <- getDefaultPlotParams(plot.type = 1)
pp$bottommargin   <- config$bottommargin
pp$topmargin      <- config$topmargin
pp$ideogramheight <- config$ideogramheight
pp$data1outmargin <- config$data1outmargin
pp$data1inmargin  <- config$data1inmargin
pp$data1height    <- config$data1height


############################
### Karyotype definition ###
############################

## In this area, we use the function plotKaryotype do draw the ideogram
## blocks based on a given GRanges object with the chromosome sizes
##
## kpDataBackground changes the color of the data panel background
## kpAddBaseNumbers is used to draw the position ticks which enable
## the visualization of the genome position (base ticks/markers)

# definition of karyotypes
kp <- plotKaryotype(
  plot.type = 1,
  main = as.character(config$plot_title),
  cex = config$title_font_size, 
  genome = genome.gr,
  plot.params = pp,
  labels.plotter = NULL
)

# The options yoffset and xoffset control the position where the names will appear
kpAddChromosomeNames(
  kp,
  cex = config$chr_name_size,
  yoffset = config$chr_name_distance
)

# background color of karyotypes
kpDataBackground(
  kp,
  color = config$karyo_color
)

# add ticks to karyotype
kpAddBaseNumbers(
  kp,
  tick.dist = config$ticks_window,
  tick.len = config$ticks_length,
  cex = config$ticks_font_size,
  add.units = TRUE
)

#######################################
### Plot features as "gene density" ###
#######################################

## In this area, we use the function kpPlotDensity to create a histogram
## of the density of genes found in a given genome window size
## data.panel = "ideogram" Is used to plot inside the ideograms

# plot density
kpPlotDensity(
  kp,
  features.gr,
  window.size = config$feat_density_window,
  data.panel = 1,
  col=config$feat_density_color,
  border=config$feat_density_color,
  r0 = 0,       # begin from where in data panel?
  r1 = 0.5      # finish where in data panel?
)

#####################################
### Plot features as "tick marks" ###
#####################################

## In this are, we use the function kpPlotRegions to draw little blocks
## that represents the genomic regions of all the features inside a given
## GRanges object
## kpAddLabels give us a little description of the data

# plot region
kpPlotRegions(
  kp,
  data=features.gr,
  col=config$feat_bars_color,
  layer.margin = 0.05, 
  border=config$feat_bars_color,
  r0=0.55,      # begin from where in data panel?
  r1=1          # finish where in data panel?
)

## Add legend
legend(
  x = "bottomright", 
  legend = c(
    config$feat_plot_label,
    config$feat_density_plot_label
  ),
  fill = c(
    config$feat_bars_color,
    config$feat_density_color
  ),
  border=NA, 
  cex = config$feat_label_font_size
)

# Finish the plot instance
dev.off()
