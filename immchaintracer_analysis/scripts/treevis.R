# ----------------------------------------------------------------------------
# treevis.R
# version: 1.25 (2021-06-21)
# author : Peter Blazso
#
# Load and visualize lineage trees
# -----------------------------------------------------------------------------

library(igraph)
library(colorspace)
library(alakazam)

divergingValueToColor <- function( x, lim=3, gran=10  )
{
  origo = gran+1
  # generate pre-tested Blue-Gray-Red palette
  colcount = gran*2+1
  pal <- diverging_hcl(n=colcount,h1=255,h2=12,c1=130,
                       cmax=80,l1=31,l2=78,p1=1,p2=1.3)
  
  # set extreme thresholds (values must fall between these two)
  x[ x > lim ] = lim
  x[ x < (-lim) ] = -lim
  
  # calculate stepping
  step <- lim/gran
  
  # lookup colors from paletteand return this vector
  return( pal[ origo+round(x/step) ] )
}

# MAIN CODE -------------------------------------------------------------------

# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
grfile <- args[1]          # input GML graph data filename - 1st arg
trfile <- args[2]          # output SVG (tree view) - 2nd arg

grfile_parts <- strsplit(basename(gsub(".gml","",grfile)), "_")[[1]]
subj    <- grfile_parts[1]    # subject ID
lineage <- grfile_parts[2]    # cellular compartment
title   <- paste0(subj," - ",lineage)


# read graph from gml file
graph <- read.graph(grfile, format="gml")

# ---- prepare graph for plot -------------------------------------------------
# correct for negative values in case of putative "Inferred" sequences
V(graph)$duplicatecount[ V(graph)$duplicatecount < 0 ] = 1

# specific vertex colors represent BASELINe sigma
V(graph)$color = divergingValueToColor( as.numeric(V(graph)$baselinesigma) )

# vertex frame colors represent compartmental affiliation
V(graph)$frame.color = "black"
V(graph)$frame.color[ V(graph)$repertoireid == "MN" ] = "#0000FF"
V(graph)$frame.color[ V(graph)$repertoireid == "AN" ] = "#FF0000"
V(graph)$frame.color[ V(graph)$repertoireid == "MB" ] = "#007D00"
#V(graph)$frame.color[ V(graph)$repertoireid == "CL" ] = "#110022"
V(graph)$frame.color[ grepl( "Inferred", V(graph)$name ) ] = "darkgrey"

# vertex size -> categorize
V(graph)$duplicatecount[ V(graph)$duplicatecount < 10  ] = 5   # small
V(graph)$duplicatecount[ V(graph)$duplicatecount >= 10 &
                         V(graph)$duplicatecount < 50  ] = 20  # medium
V(graph)$duplicatecount[ V(graph)$duplicatecount >= 50 ] = 40  # large

# constant sizes of specific vertices (putative germline, cloned, inferred seqs)
V(graph)$size = V(graph)$duplicatecount
V(graph)$size[ V(graph)$name == "1_Germline" ] = 5
V(graph)$size[ V(graph)$repertoireid == "CL" ] = 5
V(graph)$size[ grepl( "Inferred", V(graph)$name ) ] = 2

# vertex shapes
V(graph)$shape = "circle"
V(graph)$shape[ V(graph)$name == "1_Germline" ] = "square"

# edge properties
E(graph)$width = 2
E(graph)$color = "darkgrey"

# if it is a combined tree (with "order" of node groups)
# set properties to distinguish between different "orders"
if( !is.null(V(graph)$order) )
{
  # all neighbors of putative germline ("root") -> order=1
  V(graph)[ neighbors(graph, V(graph)[ "1_Germline" ]) ]$order = 1
  
  # change properties of edges based on order
  destgr <- V(graph)[ igraph::as_data_frame(graph,"edges")$to ]
#  E(graph)$width = ifelse( destgr$order == 1, 2, 1 )
  E(graph)$color[ destgr$order == 2 ] = "#bb7100"
} 
  

# create SVG file for compartmental view
svg( filename=trfile )

# if there are more than one vertex - there is a reason to plot
if( vcount(graph) > 1 )
{
    # plot graph with tree-layout
    lo <- layout_as_tree(graph, root="1_Germline")
    # make the plot (annotated)
    plot(as.undirected(graph, mode="each"),
         vertex.label.dist = 1,
         vertex.label.family="sans",
         edge.label = E(graph)$weight,
         edge.label.family="sans",
         layout=lo, main=title )

} else {
    # draw an empty canvas
    plot.new()
}

# save plot
dev.off()
