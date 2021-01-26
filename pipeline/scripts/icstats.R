# ----------------------------------------------------------------------------
# icstats.R
# version: 1.10 (2020-10-30)
# author : Peter Blazso
#
# Generate inter-compartmental statistics
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
comp <- dbfile_parts[2]    # B-cell compartment

tsv_clones <- read_tsv( dbfile )

# group clones and count them
clones <- tsv_clones %>% group_by( clone_id ) %>% 
          summarize( count=n(),
                     mn_an=ifelse(("MN" %in% repertoire_id) && ("AN" %in% repertoire_id), 1,0),
                     an_mb=ifelse(("AN" %in% repertoire_id) && ("MB" %in% repertoire_id), 1,0),
                     mb_mn=ifelse(("MB" %in% repertoire_id) && ("MN" %in% repertoire_id), 1,0),
                     mn_an_mb=ifelse(("MN" %in% repertoire_id) &&
                                     ("AN" %in% repertoire_id) &&
                                     ("MB" %in% repertoire_id), 1,0) )

# output data into a tab separated datafile
outdata <- data.frame( Repertoire=paste0(subj,"_",comp),
                       AvgClonCount=mean(clones$count),
                       MN_AN=sum(clones$mn_an)/nrow(clones),
                       AN_MB=sum(clones$an_mb)/nrow(clones),
                       MB_MN=sum(clones$mb_mn)/nrow(clones),
                       MN_AN_MB=sum(clones$mn_an_mb)/nrow(clones) )
write_tsv( outdata, outfile )
