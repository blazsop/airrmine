# ----------------------------------------------------------------------------
# geneusage_m.R
# version: 1.00 (2020-09-28)
# author : Peter Blazso
#
# Calculate V(D)J gene usages in repertoire, merge them and save as one line
# -----------------------------------------------------------------------------

library(alakazam)
library(readr)
library(dplyr)
library(tidyr)

# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
dbfile <- args[1]          # input data filename - 1st arg
gufile <- args[2]          # output gene usage .tsv filename - 2nd arg

# parse filename to extract repertoire info
dbfile_parts <- strsplit(basename(dbfile), "_")[[1]]

# read AIRR ".tsv" datafile based on the given input argument
seqdb <- read_tsv(dbfile)

# quantify gene usages and merge the results into one dataframe
v_gu <- countGenes( seqdb,gene="v_call", mode="gene",copy="duplicate_count" )
d_gu <- countGenes( seqdb,gene="d_call", mode="gene",copy="duplicate_count" )
j_gu <- countGenes( seqdb,gene="j_call", mode="gene",copy="duplicate_count" )
c_gu <- countGenes( seqdb,gene="isotype", mode="asis",copy="duplicate_count" )
vdjc_gu <- select( bind_rows( v_gu, d_gu, j_gu, c_gu ), c("seq_count","gene") )

gu <- pivot_wider( vdjc_gu, names_from = "gene", values_from = "seq_count" )

# add repertoire (sample) name column
dbfile_parts <- strsplit(basename(dbfile), "_")[[1]]
gu$repertoire <- paste0( dbfile_parts[1],"_",dbfile_parts[2] )

# write out results
write_tsv( gu, gufile )
