# -----------------------------------------------------------------------------
# treemap.R
# version: 1.0 (2020-09-28)
# author : Peter Blazso, Krisztian Csomos
#
# Generates uniformly color-coded treemap of V-gene usage/abundance
# -----------------------------------------------------------------------------

library( readr )
library( dplyr )
library( treemap )
library( colorspace )
library( alakazam )

# size of random downsampling of data
sampsize = 2500

# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
ccfile <- args[1]          # color-code file - 1st argument
dbfile <- args[2]          # AIRR data file  - 2nd argument
tmfile <- args[3]          # treemap file    - 3rd argument

# dissect input database filename
dbfile_parts <- strsplit(basename(dbfile), "_")[[1]]
subj <- dbfile_parts[1]    # subject ID
comp <- dbfile_parts[2]    # B-cell compartment

# read uniform color-code table
cc_df   <- read_tsv( ccfile )

# read AIRR data
airr_df <- read_tsv( dbfile )

# downsample AIRR data if original sample size is greater than limit
if( nrow (airr_df) > sampsize )
  airr_df <- airr_df %>% sample_n( sampsize )

# assign uniform hue values to different V-genes
airr_df$v_gene <- getGene( airr_df$v_call )
airr_df   <- airr_df %>% left_join( cc_df, by=c("v_gene"="gene") )

# mark distinct sequences in V-gene groups by (brightness) value levels
vals <-  seq( 1, 0.2, by=-0.01 )
airr_df$color <- hex(HSV(airr_df$hue, 0.8, vals))

# plot treemap
tiff( filename=tmfile, width=3600, height=2880, res=600, compression="zip" )
treemap( airr_df, index=c("v_call","sequence_id"),
         vSize="duplicate_count", vColor="color",
         type="color",
         title="", # title=paste( "AIRR TreeMap","-",paste0( subj,"_",comp ) ),
         algorithm="squarified", fontsize.labels = 0,
         border.lwds=c( 1,0.25 ),
         border.col = c("black","black") )
dev.off()
