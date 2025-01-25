# -----------------------------------------------------------------------------
# findthreshold.R
# version: 1.3 (2025-01-23)
# author : Peter Blazso
#
# Finds distance threshold that separates clonally related from unrelated
# -----------------------------------------------------------------------------

# load libraries
library(readr, quietly=T)
library(shazam, quietly=T)
library(dplyr, quietly=T)

# insignificant peak of density
peak_thr = 0.05  # below 5% height of mean peak height is noise

# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
dbfile <- args[1]          # data filename should be the 1st argument
thfile <- args[2]          # threshold data filename should be the 2nd argument

dbfile_parts <- strsplit(basename(dbfile), "_")[[1]]
subj <- dbfile_parts[1]    # subject ID
comp <- dbfile_parts[2]    # B-cell compartment

outdir   <- dirname( thfile )
outbname <- paste( outdir, paste(subj,comp,sep="_"), sep="/")

# ---- read database file & calculate length-corrected Hamming distance matrix -

# read IGMT ".tsv" datafile based on the given input argument
# and only keep a small sample fraction of it
seqdb <- read_tsv( dbfile, show_col_types=F )

if( nrow(seqdb) >= 100000 )
  seqdb <- sample_n( seqdb, 100000, replace=F )

# Use genotyped V assignments, Hamming distance, and normalize by junction length
dist <- distToNearest(seqdb, model="ham", first=FALSE, normalize="len", progress=TRUE )

# ---- determine distance threshold --------------------------------------------

# smooth noise in order to help local max determination
d <- na.omit(dist$dist_nearest)
dens <- density(d, bw=bw.bcv(d) )

# flatten insignificant noise below mean
ydif <- ifelse( dens$y < mean(dens$y)*peak_thr, 0, dens$y )

# determine position(s) of all peaks
peak_pos <- which(diff(diff(ydif)>=0)<0) + 1

# sort out significant peaks
mean_peaks <- mean(dens$y[peak_pos])
sigpeak_pos <- peak_pos[ dens$y[peak_pos] >= (mean_peaks * peak_thr) ]

# determine the location of the end of peak = threshold position
thr_pos = sigpeak_pos
thold = dens$x[ length(dens$y) ]

# in v1.3: treat the situation when there is only one significant peak
# are there more than one peaks? 
if( length(sigpeak_pos) > 1 )
{
  # yes: find threshold position -> right after the first significant peak
  thr_pos <- which(diff(dens$y[ sigpeak_pos[1]:length(dens$y) ]) >= 0 )[[1]] + sigpeak_pos[1]
  # determine threshold
  thold <- dens$x[ thr_pos ]
}
# no: keep the default value = end position of the range


# draw plot and put this into a tiff file
tiff( filename=paste(outbname,"threshold_histo.tiff", sep="_"),
      width=2400, height=1200, res=300, compression="zip" )
plot(dens);
abline( v=thold, col="blue", lty="dashed" )
abline( h=mean_peaks * peak_thr, col="brown", lty="dotted")
dev.off()

# output threshold value into a textfile
thfile <- file( thfile )
writeLines( toString(thold), con=thfile )
close( thfile )
