# ----------------------------------------------------------------------------
# R_S_mutations.R
# version: 1.10 (2020-10-30)
# author : Boglarka Ujhazi, Peter Blazso
#
# Calculates R and S mutation rate/frequency
# -----------------------------------------------------------------------------

library(alakazam)
library(shazam)
library(readr)
library(dplyr)
library(openxlsx)

# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
dbfile  <- args[1]         # input TSV filename (with germlines) - 1st argument
outfile <- args[2]         # output XLSX - 2nd argument

# load data
db <- read_tsv( dbfile )

# calculate R and S mutation counts and sort by isotype
#dbcnt <- observedMutations( db,
#                            sequenceColumn="sequence_alignment",
#                            germlineColumn="germline_alignment_d_mask",
#                            regionDefinition=NULL,
#                            frequency=FALSE,
#                            nproc=1 ) %>%
#         arrange( isotype )

# calculate R and S mutation frequencies and sort by isotype
# dbfrq <- observedMutations( db,
#                             sequenceColumn="sequence_alignment",
#                             germlineColumn="germline_alignment_d_mask",
#                             regionDefinition=NULL,
#                             frequency=TRUE,
#                             nproc=1 ) %>%
#          arrange( isotype )

# calculate combined R and S mutation frequencies and sort by isotype
dbcmb <- observedMutations( db,
                            sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            frequency=FALSE,
                            combine=TRUE,
                            nproc=1) %>%
         arrange( isotype )

avgtab=data.frame( Average=mean(dbcmb$mu_count) )


# create Excel workbook
wb <- createWorkbook()
# addWorksheet(wb, "Count")
# addWorksheet(wb, "Frequency")
addWorksheet(wb, "ALL_combined")
addWorksheet(wb, "Average_combined")


# write out data into Excel file
# writeDataTable( wb, "Count",     dbcnt )
# writeDataTable( wb, "Frequency", dbfrq )
writeDataTable( wb, "ALL_combined",  dbcmb )
writeDataTable( wb, "Average_combined", avgtab )

# write out results
saveWorkbook( wb, outfile, overwrite=T )


#Plot
# g2 <- ggplot(db_obs_com, aes(x=isotype, y=mu_freq, fill=isotype)) + theme_bw() + ggtitle("Total mutations") +
# xlab("Isotype") + ylab("Mutation frequency") + scale_fill_manual(name="Isotype", values=IG_COLORS) + geom_boxplot()
# plot(g1)



