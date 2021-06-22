# ----------------------------------------------------------------------------
# clstats.R
# version: 1.00 (2020-10-09)
# author : Peter Blazso
#
# Generate lineage statistics
# -----------------------------------------------------------------------------

library(readr)
library(dplyr)

# input variables (mostly generated from parsed input filename)
args    <- commandArgs(TRUE)  # read arguments into "args" vector
dbfile  <- args[1]            # data filename should be the 1st argument
outfile <- args[2]            # output data filename should be the 2nd argument

# file, repertoire names, directories, etc.
dbfile_parts <- strsplit(basename(dbfile), "_")[[1]]
subj <- dbfile_parts[1]    # subject ID
line <- dbfile_parts[2]    # B-cell lineage

tsv_lin <- read_tsv( dbfile )

# group clones and count them
lincounts <- tsv_lin %>% group_by( repertoire_id ) %>% 
                         summarize( count=n(), copy=sum(duplicate_count) )

# output data into a tab separated datafile
outdata <- data.frame( Lineage=paste0(subj,"_",line),
                       MN_seq=lincounts$count[ lincounts$repertoire_id == "MN" ][1],
                       MN_copy=lincounts$copy[ lincounts$repertoire_id == "MN" ][1],
                       AN_seq=lincounts$count[ lincounts$repertoire_id == "AN" ][1],
                       AN_copy=lincounts$copy[ lincounts$repertoire_id == "AN" ][1], 
                       MB_seq=lincounts$count[ lincounts$repertoire_id == "MB" ][1],
                       MB_copy=lincounts$copy[ lincounts$repertoire_id == "MB" ][1] )
write_tsv( outdata, outfile )
