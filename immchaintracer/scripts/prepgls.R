# -----------------------------------------------------------------------------
# prepgls.R
# version: 1.10 (2021-06-21)
# author : Peter Blazso
#
# Prepare raw AIRR germline list for subsequent analyses
# -----------------------------------------------------------------------------

library(readr)
library(dplyr)
library(alakazam)

# columns that need to be present and in the same exact order 
match_cols <- c( "sequence_id", "sequence", "productive", "junction",
                 "locus", "isotype", "v_call", "j_call",
                 "junction_length", "np1_length", "np2_length",
                 "sequence_alignment", "d_call", 
                 "v_sequence_start", "v_sequence_end", 
                 "d_sequence_start", "d_sequence_end",
                 "j_sequence_start", "j_sequence_end",
                 "v_germline_start", "v_germline_end", 
                 "d_germline_start", "d_germline_end",
                 "j_germline_start", "j_germline_end",
                 "repertoire_id", "duplicate_count"  )

# -----------------------------------------------------------------------------

# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
dbfile  <- args[1]     # input TSV germline database - 1st argument
outfile <- args[2]     # PREPARED TSV database - 2nd argument

# read data files into tibbles
rep <- read_tsv(dbfile)

# specify "isotype", "duplicate_count" and "repertoire_id" directly
rep$isotype <- "hIGHM"
rep$duplicate_count <- 1
rep$repertoire_id <- "CL"

# concatenate CDR1+FWR2+CDR2+FWR3+CDR3+FWR4 = "trim" (or overwrite) FWR1
rep$sequence_alignment <- paste0( gsub(".", "\\.", rep$fwr1),
                                  rep$cdr1, rep$fwr2,
                                  rep$cdr2, rep$fwr3,
                                  rep$cdr3, rep$fwr4 )

# clear placeholders (".") to obtain the raw, concatenated sequence
rep$sequence <- gsub( ".", "", rep$sequence_alignment, fixed=T )


# select only a subset of columns to get in sync with repertoire datasets
rep <- rep %>% select( all_of(match_cols) )

# write out the results into a new database
write_tsv(rep, path=outfile, quote_escape=F)
