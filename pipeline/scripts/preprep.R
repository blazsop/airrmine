# -----------------------------------------------------------------------------
# preprep.R
# version: 1.00 (2020-09-21)
# author : Peter Blazso, Krisztian Csomos
#
# Prepare raw AIRR for subsequent analyses
#   1. merge datasets (with externals like isotypes, counts)
#   2. pre-collapse (collapse identical sequences)
#   3. collapse (collapse genetically close sequences)
#   4. subsample (randomly choose specified number of prepared sequences)
# -----------------------------------------------------------------------------

library(readr)
library(dplyr)
library(alakazam)

# number of randomly sampled productive sequences
sample_no = 20000

# core proterties of sequences that must be matched for either collapsing
# or pre-collapsing
match_core_props <- c("locus", "isotype", "v_call", "j_call",
                      "junction_length", "np1_length", "np2_length" )

match_extr_props <- c( "sequence", "sequence_alignment",
                       "productive", "junction", "d_call",
                       "v_sequence_start", "v_sequence_end",
                       "d_sequence_start", "d_sequence_end",
                       "j_sequence_start", "j_sequence_end",
                       "v_germline_start", "v_germline_end",
                       "d_germline_start", "d_germline_end",
                       "j_germline_start", "j_germline_end",
                       "repertoire_id" )


# -----------------------------------------------------------------------------
isCloseSequence <- function( seq1, seq2, dist=1 )
{
  if( nchar(seq1) == nchar(seq2) ) return (seqDist(seq1, seq2) <= dist)
  else                             return (FALSE)
}

# -----------------------------------------------------------------------------
findClosestSequence <- function( seq_tbl, rep_tbl, match_vars, dist=1,
                                 id="sequence_id", seq="sequence", count="duplicate_count" )
{
  # "seq_tbl" is a one-row tibble of a single, chosen sequence
  # extract close-to-"seq_tbl" group of sequences from "rep_tbl" where
  # all "match_vars" variables match
  extr_tbl <- semi_join( rep_tbl, seq_tbl, match_vars )

  # if there is no match, just return with NA
  if( dim(extr_tbl)[1] == 0 ) return (NA)

  # pick the subset of the matched group of sequences that
  # fall within "dist" Hamming distance
  extr_tbl <- extr_tbl[ sapply( pull(extr_tbl, seq),
                                isCloseSequence,
                                seq_tbl[[1,seq]], dist ), ]

  # if there is no match, just return with NA
  if( dim(extr_tbl)[1] == 0 ) return (NA)

  # if there is more than one sequence in this subset
  # choose the one with the highest count
  extr_tbl <- top_n( extr_tbl, 1, count )

  # if there is no match, just return with NA, otherwise return with sequence ID
  if( dim(extr_tbl)[1] == 0 ) return (NA)
  else return( extr_tbl[[1, id]] )
}

# -----------------------------------------------------------------------------
collapseCloseSequences <- function( rep, dist=1, thr_count=1, discard_low_copy=TRUE )
{
  #  pre-collapse identical sequences in the whole repertoire
  #  with matching specified properties
  rep.pc <- rep %>% group_by_at( .vars=c( "sequence", "productive", "junction",
                                 match_core_props, match_extr_props ) ) %>%
                    summarize( duplicate_count = sum(duplicate_count) )

  # add a unique ID for each sequence
  rep.pc$sequence_id = 1:nrow(rep.pc)

  # split pre-collapsed repertoire into low- and high-copy sequences
  rep.lcopy <- ungroup( rep.pc[ rep.pc$duplicate_count <= thr_count, ] )
  rep.hcopy <- ungroup( rep.pc[ rep.pc$duplicate_count >  thr_count, ] )
  rm(rep.pc)

  # find the closest sequence from high-copy sequence pool to all low-copy
  # sequences and if available, add its external ID ("id_close") for collapsing
  # CAUTION: this is a time consuming process! (60k sequence ~ 13 minutes)
  rep.lcopy$id_close <- as.vector( by( rep.lcopy, rep.lcopy$sequence_id,
                                       findClosestSequence,
                                       rep.hcopy, match_core_props ) )

  # collect and summarize counts of close low-copy sequences matched to
  # high copy sequences
  rep.dest <- ungroup(rep.lcopy) %>%
              filter( !is.na(id_close) ) %>%
              select( id_close, duplicate_count ) %>%
              group_by( id_close ) %>%
              summarize( duplicate_count_s=sum(duplicate_count)) %>%
              right_join( ungroup(rep.hcopy), by=c("id_close"="sequence_id") ) %>%
              rename( sequence_id=id_close ) %>%
              mutate( duplicate_count=duplicate_count+
                                      ifelse(is.na(duplicate_count_s),
                                             0,duplicate_count_s) )%>%
              select( -duplicate_count_s )

  # add unmatched low-copy sequences
  if( !discard_low_copy )
    rep.dest <- rep.lcopy %>%
                filter( is.na(id_close) ) %>%
                select( -id_close ) %>%
                bind_rows( rep.dest )

  return (rep.dest)
}


# -----------------------------------------------------------------------------

# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
tab_file    <- args[1]     # TAB database filename should be the 1st argument
irep_file   <- args[2]     # IREP database filename should be the 2nd argument
newtab_file <- args[3]     # COLLAPSED TAB database filename should be the 3rd

# read data files into tibbles
rep.tsv <- read_tsv(tab_file)
rep.csv <- read_csv(irep_file)

# convert all sequences to fully UPPERCASE
rep.csv$joinedSeq <- toupper( rep.csv$joinedSeq )

# add "isotype" and "duplicate_count" from an external (IRepertoire) source
rep <- select( rep.csv, joinedSeq, isotype=C, duplicate_count=copy ) %>%
       inner_join( rep.tsv, by=c("joinedSeq"="sequence") ) %>%
       select( -joinedSeq )
rm( rep.tsv, rep.csv )

# add "repertoire_id" from filename
dbfile_parts <- strsplit(basename(tab_file), "_")[[1]]
rep$repertoire_id = dbfile_parts[2]    # B-cell compartment

# concatenate CDR1+FWR2+CDR2+FWR3+CDR3+FWR4 = "trim" (or overwrite) FWR1
rep$sequence_alignment <- paste0( gsub(".", "\\.", rep$fwr1),
                                  rep$cdr1, rep$fwr2,
                                  rep$cdr2, rep$fwr3,
                                  rep$cdr3, rep$fwr4 )
# clear placeholders (".") to obtain the raw, concatenated sequence
rep$sequence <- gsub( ".", "", rep$sequence_alignment, fixed=T )

# collapse genetically very close sequences with default thresholds
# to eliminate, compensate for the majority of PCR errors
rep.collapsed <- collapseCloseSequences( rep, thr_count=4, discard_low_copy=T ) %>%
                 filter( productive == TRUE )

# write out the results into a new database
write_tsv(rep.collapsed, path=newtab_file, quote_escape=F)
