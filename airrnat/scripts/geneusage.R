# ----------------------------------------------------------------------------
# geneusage.R
# version: 1.00 (2020-09-28)
# author : Peter Blazso
#
# Calculates V(D)J gene and isotype usages in repertoire
# -----------------------------------------------------------------------------

library(alakazam)
library(openxlsx)

# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
dbfile <- args[1]          # input data filename - 1st arg
gufile <- args[2]          # output gene usage XLSX filename - 2nd arg

# read AIRR ".tsv" datafile based on the given input argument
seqdb <- read.delim(dbfile, stringsAsFactors=F)

# create Excel workbook
wb <- createWorkbook()
addWorksheet(wb, "V genes")
addWorksheet(wb, "D genes")
addWorksheet(wb, "J genes")

# quantify gene usages and put the results into separate worksheets
writeDataTable( wb, "V genes", countGenes(seqdb,gene="v_call",
                                          mode="gene",copy="duplicate_count") )
writeDataTable( wb, "D genes", countGenes(seqdb,gene="d_call",
                                          mode="gene",copy="duplicate_count") )
writeDataTable( wb, "J genes", countGenes(seqdb,gene="j_call",
                                          mode="gene",copy="duplicate_count") )

# write out results
saveWorkbook( wb, gufile, overwrite=T )
