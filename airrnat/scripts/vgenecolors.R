# -----------------------------------------------------------------------------
# vgenecolors.R
# version: 1.0 (2020-09-21)
# author : Peter Blazso
#
# Collects V-genes from specified repertoires and assigns unified color codes
# -----------------------------------------------------------------------------

library( openxlsx )
library( dplyr )
library( readr )
library( colorspace )

# arguments
args    <- commandArgs(TRUE)
outfile <- args[1]  # 1st argument is destination datafile

# all other arguments are input gene usage files
# read in and merge all data into one tibble
df <- tibble()
for( arg in args[2:length(args)] )
{
  # read sheet #1 (with V-gene usage list) from next XLSX file
  df <- bind_rows( df, read.xlsx( arg, sheet=1 ) )
}

# make a unique list of V-genes
df <- df %>% select( gene ) %>% group_by( gene ) %>% summarize() %>% ungroup()

# generate additional id to sort the genes correctly
df$id <- sprintf( "%s%03d",
                  unlist( lapply(strsplit(df$gene,"-"),getElement,1)),
                  as.numeric( 
                    unlist( lapply(strsplit(df$gene,"-"),getElement,2) ))
                )
# sort gene names
df <- df %>% arrange( id ) %>% select( -id )

# assign hue and sample color values to this sorted list
df$hue   <- round(seq(1, 360, length.out = nrow(df)))
df$color <- hex(HSV(df$hue, 0.8, 1))

# write out color-coded list
write_tsv( df, outfile )
