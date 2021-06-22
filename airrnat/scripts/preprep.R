# -----------------------------------------------------------------------------
# preprep.R
# version: 1.00 (2020-09-21)
# author : Peter Blazso, Krisztian Csomos
#
# Prepare aligned/annotated AIRR - add isotypes, counts from an external source
# -----------------------------------------------------------------------------

library(readr)
library(dplyr)

# -----------------------------------------------------------------------------
# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
tab_file    <- args[1]     # TAB database filename should be the 1st argument
irep_file   <- args[2]     # IREP database filename should be the 2nd argument
newtab_file <- args[3]     # PREPARED TAB database filename should be the 3rd

# read data files into tibbles
rep.tsv <- read_tsv(tab_file)
rep.csv <- read_csv(irep_file)

# keep important columns and convert all sequences to fully UPPERCASE
rep.csv <- select( rep.csv, joinedSeq, isotype=C, duplicate_count=copy )
rep.csv$joinedSeq <- toupper( rep.csv$joinedSeq )

# add "isotype" and "duplicate_count" from an external (IRepertoire) source
rep <- inner_join( rep.tsv, rep.csv, by=c("sequence"="joinedSeq") )
rm( rep.tsv, rep.csv )

# add "repertoire_id" from filename
dbfile_parts <- strsplit(basename(tab_file), "_")[[1]]
rep$repertoire_id = dbfile_parts[2]    # B-cell compartment

# write out the results into a new database
write_tsv(rep, path=newtab_file, quote_escape=F)
