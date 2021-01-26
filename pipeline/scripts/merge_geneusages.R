# -----------------------------------------------------------------------------
# merge_geneusages.R
# version: 1.0 (2020-09-28)
# author : Peter Blazso
#
# Collects and merges gene usage data from specified merged VDJ counts
# -----------------------------------------------------------------------------

library( dplyr )
library( readr )

# arguments
args    <- commandArgs(TRUE)
outfile <- args[1]  # 1st argument is destination datafile

# all other arguments are input gene usage files
# read in and merge all data into one tibble
df <- tibble()
for( arg in args[2:length(args)] )
{
  df <- bind_rows( df, read_tsv( arg ) )
}

# write out list
write_tsv( df, outfile )
