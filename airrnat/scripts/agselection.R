# -----------------------------------------------------------------------------
# agselection.R
# version: 1.0 (2020-09-25)
# author : Peter Blazso
#
# Quantification of antigene selection pressure using BASELINe method
# ref.: https://shazam.readthedocs.io/en/stable/vignettes/Baseline-Vignette/
# -----------------------------------------------------------------------------

library( readr )
library( dplyr )
library( shazam )

# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
dbfile  <- args[1]         # AIRR data file  - 1st argument
outfile <- args[2]         # BASELINe extended file    - 2nd argument

# read database file
indb <- read_tsv( dbfile )

# calculate the posterior probability density function (PDF) for each sequence
db <- calcBaseline(indb,
                   sequenceColumn = "sequence_alignment",
                   germlineColumn = "germline_alignment_d_mask",
                   testStatistic  = "focused")

# calculate BASELINe summary statistics for each sequence
db <- summarizeBaseline(db, returnType = "df")
db$baseline_sigma     <- ifelse(is.na(db$baseline_sigma),0,db$baseline_sigma)
db$baseline_ci_lower  <- ifelse(is.na(db$baseline_ci_lower),0,db$baseline_ci_lower)
db$baseline_ci_upper  <- ifelse(is.na(db$baseline_ci_upper),0,db$baseline_ci_upper)
db$baseline_ci_pvalue <- ifelse(is.na(db$baseline_ci_pvalue),0,db$baseline_ci_pvalue)


# merge input AIRR database with BASELINe calculations
outdb <- indb %>% left_join( db ) %>% select( -region )

# write out data
write_tsv( outdb, outfile )