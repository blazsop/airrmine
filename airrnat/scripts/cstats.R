# ----------------------------------------------------------------------------
# cstats.R
# version: 1.00 (2020-10-09)
# author : Peter Blazso
#
# Generate (intra-)compartmental statistics (counts, D50, Shannon, etc.)
# -----------------------------------------------------------------------------

library(readr)
library(dplyr)
library(shazam)

# input variables (mostly generated from parsed input filename)
args    <- commandArgs(TRUE)  # read arguments into "args" vector
dbfile  <- args[1]            # data filename should be the 1st argument
outfile <- args[2]            # output data filename should be the 2nd argument

# file, repertoire names, directories, etc.
dbfile_parts <- strsplit(basename(dbfile), "_")[[1]]
subj <- dbfile_parts[1]    # subject ID
comp <- dbfile_parts[2]    # immune cell compartment

# load clonally-clustered (!) repertoire data
db <- read_tsv( dbfile )

# unique/total index
ut <- nrow(db) / sum(db$duplicate_count)

# sort sequence counts by abundance, most abundant goes first
ser = sort( db$duplicate_count, decreasing = T )

# total sequence number
tsno   <- nrow( db )

# total copy number
tcno   <- sum( db$duplicate_count )

# average clone size = group clones, count the sequences and get their mean
clones <- db %>% group_by( clone_id ) %>% summarize( count=n() )
acs    <- mean(clones$count)

# D50 index
sumi = 0
dc50 = round( sum( ser ) / 2 )
for( i in 1:length(ser) )
{
  sumi = sumi + ser[i]
  if( sumi >= dc50 ) break
}
d50 <- i/length(ser) * 100

limit = 1000
cap <- ifelse( length(ser) >= limit, limit, length(ser) )
ser1k  <- ser[ 1:cap ]

# Shannon-Wiener Diversity Index
swdi <- -sum(      db$duplicate_count/sum(db$duplicate_count) *
              log( db$duplicate_count/sum(db$duplicate_count), base=exp(1) ) )

swdi1k <- -sum( ser1k/sum(ser1k) *
                log( ser1k/sum(ser1k), base=exp(1) ) )

# Simpson's Diversity Index
sdi <- 1 - (sum( db$duplicate_count * (db$duplicate_count - 1) )/
           (sum( db$duplicate_count ) * ( sum(db$duplicate_count)-1) ))

sdi1k <- 1 - ( sum(ser1k * (ser1k - 1))/
               (sum(ser1k) * (sum(ser1k)-1)) )

# frequency of top 100
ft100 <- sum( ser[1:100] ) / sum( db$duplicate_count )

# frequency of top 10
ft10 <- sum( ser[1:10] ) / sum( db$duplicate_count )

# distinct frequencies of top 10
top1  <- ser[1] / sum( db$duplicate_count )
top2  <- ser[2] / sum( db$duplicate_count )
top3  <- ser[3] / sum( db$duplicate_count )
top4  <- ser[4] / sum( db$duplicate_count )
top5  <- ser[5] / sum( db$duplicate_count )
top6  <- ser[6] / sum( db$duplicate_count )
top7  <- ser[7] / sum( db$duplicate_count )
top8  <- ser[8] / sum( db$duplicate_count )
top9  <- ser[9] / sum( db$duplicate_count )
top10 <- ser[10] / sum( db$duplicate_count )

outdata <- data.frame( Repertoire=paste0(subj,"_",comp),
                       TotalSeq=tsno, TotalCopy=tcno,
                       AvgClonCount=acs, U_T=ut, D50=d50,
                       Shannon=swdi, Shannon_top1k=swdi1k,
                       Simpson=sdi, Simpson_top1k=sdi1k,
                       Top100all=ft100, Top10all=ft10,
                       T1=top1, T2=top2, T3=top3, T4=top4, T5=top5,
                       T6=top6, T7=top7, T8=top8, T9=top9, T10=top10 )

# output data into a tab separated datafile
write_tsv( outdata, outfile )