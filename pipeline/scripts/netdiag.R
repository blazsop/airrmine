# ----------------------------------------------------------------------------
# netdiag.R
# version: 1.00 (2020-09-08)
# author : Peter Blazso, Krisztian Csomos
#
# Generate clonal network diagrams (plotted as multiple graphs)
# -----------------------------------------------------------------------------

library(parallel)
library(alakazam)
library(igraph)
library(dplyr)

# -----------------------------------------------------------------------------
# GLOBAL CONSTANTS
# -----------------------------------------------------------------------------
cores_to_use = 16
dnapars_bin  = "/usr/lib/phylip/bin/dnapars"
#dnapars_bin  = "~/software/bin/dnapars"
logging      = TRUE

sampdepth    = 2500  # subsampling "depth"


# -----------------------------------------------------------------------------
# PREDEFINED FUNCTIONS
# -----------------------------------------------------------------------------

# wrap a rename mechanism for "buildPhylipLineage" to make unique vertex names
# and generate single-vertex graphs, as well
buildUPhylipLineage <- function( clone_in, dnapars_bin,
                                 dist_mat=getDNAMatrix(gap=0),
                                 rm_temp=FALSE, verbose=FALSE)
{
  # there are more than one sequence -> generate graph
  if( nrow(clone_in@data) > 1 )
  {
    graph <- buildPhylipLineage(clone=clone_in, phylip_exec=dnapars_bin,
                                dist_mat=dist_mat,
                                rm_temp=rm_temp, verbose=verbose  )

    # make vertex names/labels unique ("Inferred..." -> "cloneid_Inferred...")
    V(graph)$name  <- paste0( clone_in@clone, "_", V(graph)$name  )
    V(graph)$label <- paste0( clone_in@clone, "_", V(graph)$label )
    V(graph)$clone_id <- clone_in@clone
  }
  # there is only one seqence ->
  # create a single-vertex igraph object without edges
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
                          name      = uid,
                          sequence  = clone_in@data$sequence[1],
                          collapse_count  = clone_in@data$collapse_count[1],
                          duplicate_count = clone_in@data$duplicate_count[1],
                          isotype   = clone_in@data$isotype[1],
                          repertoire_id = clone_in@data$repertoire_id[1],
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
dbfile  <- args[1]         # data filename should be the 1st argument
gmlfile <- args[2]         # output GML filename should be the 2nd argument

dbfile_parts <- strsplit(basename(dbfile), "_")[[1]]
subj <- dbfile_parts[1]    # subject ID
comp <- dbfile_parts[2]    # cellular compartment
graph_fbase <- paste(subj,comp,sep="_")

# ---- read database file
# read AIRR datafile based on the given input argument
# and only keep a small sample fraction of it
seqdb <- read.delim(dbfile, stringsAsFactors=F)

# ---- subsample data
sampdepth <- ifelse( nrow(seqdb) > sampdepth, sampdepth, nrow(seqdb) )
seqdb <- sample_n( seqdb, sampdepth )


# ---- make seq IDs unique ("index")
seqdb$sequence_id = paste0( "ID_", seqdb$sequence_id )

# ---- group by clones
clones <- seqdb %>% group_by(clone_id) %>%
          do(CHANGEO=makeChangeoClone(., num_fields="duplicate_count",
                                         seq_fields=c("isotype", "repertoire_id"),
                                      germ="germline_alignment_d_mask",
                                      verbose=logging))

# ---- build unique lineages of grouped clones
graphs <- mclapply(clones$CHANGEO, buildUPhylipLineage,
                   dnapars_bin, rm_temp=T,
                   verbose=logging, mc.cores = cores_to_use)

# ---- merge vector of igraph objects into one combined object for plotting
sumgraph    <- make_empty_graph( directed=F )
sumvertices <- igraph::as_data_frame(sumgraph, what="vertices")
sumedges    <- igraph::as_data_frame(sumgraph, what="edges")
for( graph in graphs )
{
  sumvertices <- rbind( sumvertices, igraph::as_data_frame(graph, what="vertices") )
  sumedges    <- rbind( sumedges,    igraph::as_data_frame(graph, what="edges") )
}
sumgraph <- graph_from_data_frame( sumedges, directed=F, vertices=sumvertices )

# remove putative (fictive) germline vertices/edges to represent reality better
graph <- delete_vertices(sumgraph,
                         V(sumgraph)$name[ grepl("Germline",V(sumgraph)$name) ])

# ---- save raw graph data in the output file
write.graph(graph, file=gmlfile, format="gml")
