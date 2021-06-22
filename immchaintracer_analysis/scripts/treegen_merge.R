# ----------------------------------------------------------------------------
# treegen_merge.R
# version: 1.00 (2020-11-14)
# author : Peter Blazso
#
# Merge related clonotypes and generate a combined lineage tree as graph
# -----------------------------------------------------------------------------

library(parallel)
library(alakazam)
library(igraph)
library(dplyr)
library(readr)

# -----------------------------------------------------------------------------
# GLOBAL CONSTANTS
# -----------------------------------------------------------------------------
cores_to_use = 16
dnapars_bin  = "/usr/lib/phylip/bin/dnapars"
# dnapars_bin  = "~/software/bin/dnapars"
logging      = FALSE


# -----------------------------------------------------------------------------
# PREDEFINED FUNCTIONS
# -----------------------------------------------------------------------------

# wrap a rename mechanism for "buildPhylipLineage" to make unique vertex names
# and generate single-vertex graphs, as well
buildUPhylipLineage <- function( clone_in, dnapars_bin,
                                 dist_mat=getDNAMatrix(gap=0), 
                                 rm_temp=FALSE, verbose=FALSE)
{
  # there is more than one sequences -> generate graph
  if( nrow(clone_in@data) > 1 )
  {
    graph <- buildPhylipLineage(clone=clone_in, phylip_exec=dnapars_bin,
                                dist_mat=dist_mat,
                                rm_temp=rm_temp, verbose=verbose  )

    # make vertex names/labels unique ("Inferred..." -> "cloneid_Inferred...")
    V(graph)$name  <- paste0( clone_in@clone, "_", V(graph)$name  )
    V(graph)$label <- paste0( clone_in@clone, "_", V(graph)$label )
    V(graph)$clone_id <- clone_in@clone

    # if vertices are inferred, BASELINe parameters are NA, convert to 0
    V(graph)$baseline_sigma = ifelse( is.na(V(graph)$baseline_sigma),0,
                                            V(graph)$baseline_sigma ) 
  }
  # there is only one seqence ->
  # create an single-vertex igraph object without edges
  else
  {
    # assemble vertex name/label
    uid <- paste0( clone_in@clone,"_",clone_in@data$sequence_id[1] )

    # generate a single-vertex "graph" out of a single sequence clone
    graph <- make_empty_graph( directed=F ) %>%
             add_vertices(1,
                          clone     = clone_in@data$clone[1],
                          v_gene    = clone_in@data$v_gene[1],
                          j_gene    = clone_in@data$j_gene[1],
                          junc_len  = clone_in@data$junc_len[1],
                          order     = clone_in@data$order[1],
                          name      = uid,
                          sequence  = clone_in@data$sequence[1],
                          collapse_count  = clone_in@data$collapse_count[1],
                          duplicate_count = clone_in@data$duplicate_count[1],
                          isotype   = clone_in@data$isotype[1],
                          repertoire_id = clone_in@data$repertoire_id[1],
                          baseline_sigma = clone_in@data$baseline_sigma[1],
                          baseline_ci_pvalue = clone_in@data$baseline_ci_pvalue[1],
                          clone_id  = clone_in@clone,
                          label     = uid )
  }

  return( graph )
}

# -----------------------------------------------------------------------------
# MAIN PROGRAM CODE
# -----------------------------------------------------------------------------
# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
db1file <- args[1]         # 1st data file should be the 1st argument (earlier clonotype)
db2file <- args[2]         # 2nd data file should be the 2nd argument (later clonotype)
gmlfile <- args[3]         # output GML filename should be the 3rd argument

dbfile_parts <- strsplit(basename(db1file), "_")[[1]]
subj <- dbfile_parts[1]    # subject ID
comp <- dbfile_parts[2]    # cellular compartment
graph_fbase <- paste(subj,comp,sep="_")

# ---- read and merge database files
# read IGMT ".tab" datafiles based on the given input arguments and take cloned
# sequence ("CL") only once
seq1db <- read_tsv(db1file)
seq1db$order = 1
seq2db <- read_tsv(db2file) %>% dplyr::filter( repertoire_id != "CL" )
seq2db$order = 2
seqdb <- bind_rows( seq1db, seq2db )
rm( seq1db, seq2db )

# ---- make seq IDs unique ("index")
seqdb$sequence_id = paste0( "ID_", seqdb$sequence_id )

# change "." placeholders to "N"s (for dnapars)
seqdb$sequence <- gsub( "\\.", "N", seqdb$sequence_alignment )
clone <- new("ChangeoClone",
             data=seqdb,
             clone="1",
             germline=gsub( "\\.","N",seqdb$germline_alignment[1] ),
             v_gene=getGene(seqdb$v_call[1]),
             j_gene=getGene(seqdb$j_call[1]),
             junc_len=seqdb$junction_length[1] )

# build tree: call dnapars
graph <- buildUPhylipLineage( clone, dnapars_bin, rm_temp=T, verbose=logging )

# ---- save raw graph data in a file
write.graph(graph, file=gmlfile, format="gml")

