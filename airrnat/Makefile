# ------------------------------------------------------------------------------
# AIRRNAT analytic pipeline system
#
# by Peter Blazso
# 2021-06-20
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# GLOBAL SETTINGS
# ------------------------------------------------------------------------------
SHELL     := /bin/bash   # DO NOT EDIT! Using "bash" is mandatory!

# immune receptor chain type
LOCUS     := Ig
SUBLOCUS  := H
CHAIN     := $(LOCUS)$(SUBLOCUS)

# species
SPECIES    := human
SPECIES_L  := Homo_sapiens

# directories containg reference data and code
# IMGT and IgBLAST's reference database files
IGDATADIR := ../igdata
# data processing scripts
SCRIPTDIR := scripts

# source/destination data directories 
# IMGT annotated .txz files
IMGTDIR   := imgt
# IRepertoire dataset of annotasted .csv files
IREPDIR   := irep
# FASTA extraction
FASTADIR  := fasta
# raw, annotated AIRR dataset (without external data)
TSVDIR    := tsv
# annotated and extended (merged) AIRR dataset
PREPDIR   := prep
# collapsed AIRR dataset (source of all analytics)
COLLDIR   := coll
# destination for processed and supplementary data
DESTDIR   := dest
# clonally clustered AIRR dataset
CLONEDIR  := clone
# background data and figures of AIRR network diagrams
GRAPHDIR  := graph

# filename suffices
IMGT_SUFF := .txz
IREP_SUFF := _irep.csv
PREP_SUFF := _prepared.tsv
COLL_SUFF := _collapsed.tsv

# source (input) dataset files - where the pipeline begins
# Note: uncomment the line-pairs below which are applicable.
#IMGT      = $(wildcard $(IMGTDIR)/*$(IMGT_SUFF)) # individual IMGT output files
#REPS      = $(notdir $(subst $(IMGT_SUFF),,$(IMGT))) # raw repertoire names

IREP      = $(wildcard $(IREPDIR)/*$(IREP_SUFF)) # individual IREP output files
REPS      = $(notdir $(subst $(IREP_SUFF),,$(IREP))) # raw repertoire names

#PREP      = $(wildcard $(PREPDIR)/*$(PREP_SUFF)) # individual PREP output files
#REPS      = $(notdir $(subst $(PREP_SUFF),,$(PREP))) # raw repertoire names

#COLL      = $(wildcard $(COLLDIR)/*$(COLL_SUFF)) # individual COLL output files
#REPS      = $(notdir $(subst $(PREP_SUFF),,$(COLL))) # raw repertoire names

# destination datafile lists, single (aggregated) data files
#EXTRFAS   = $(addprefix $(FASTADIR)/, $(addsuffix .fasta, $(REPS)))
#ANNOTATE  = $(addprefix $(TSVDIR)/,   $(addsuffix .tsv, $(REPS)))
#EXTRTSV   = $(addprefix $(TSVDIR)/,   $(addsuffix _db-pass.tsv, $(REPS)))
PREPREP   = $(addprefix $(PREPDIR)/,  $(addsuffix $(PREP_SUFF), $(REPS)))
COLLREP   = $(addprefix $(COLLDIR)/,  $(addsuffix $(COLL_SUFF), $(REPS)))
GENEUSG   = $(addprefix $(DESTDIR)/,  $(addsuffix _genes.xlsx, $(REPS)))
COLCODES  = $(DESTDIR)/v_gene_colorcodes.tsv
DISTTH    = $(addprefix $(DESTDIR)/,  $(addsuffix _dist_thrs.data, $(REPS)))
TREEMAPS  = $(addprefix $(DESTDIR)/,  $(addsuffix _tmap.tiff, $(REPS)))
DEFCLONES = $(addprefix $(CLONEDIR)/, $(addsuffix _allclone.tsv, $(REPS)))
TMAPS_CL  = $(addprefix $(DESTDIR)/,  $(addsuffix _tmap_cl.tiff, $(REPS)))
CSTATS    = $(addprefix $(CLONEDIR)/, $(addsuffix _cstats.tsv, $(REPS)))
GENGLS    = $(addprefix $(CLONEDIR)/, $(addsuffix _withgl.tsv, $(REPS)))
RSMUT     = $(addprefix $(CLONEDIR)/, $(addsuffix _rsmut.xlsx, $(REPS)))
AGSEL     = $(addprefix $(CLONEDIR)/, $(addsuffix _agsel.tsv, $(REPS)))
CGRAPH    = $(addprefix $(GRAPHDIR)/, $(addsuffix .gml, $(REPS)))
VISND     = $(addprefix $(GRAPHDIR)/, $(addsuffix _AIRRND.svg, $(REPS)))

# ------------------------------------------------------------------------------
.PHONY: all clean

all: geneusg treemaps tmaps_cl cstats rsmut agsel visnd
	@echo done.

# ------------------------------------------------------------------------------
# CONVERT (EXTRACT) IMGT annotated DB files into TSV databases (optional)
# ------------------------------------------------------------------------------
#$(TSVDIR)/%_db-pass.tsv: $(IMGTDIR)/%.txz
#	@echo Converting IMGT database $< into TSV file...
#	@mkdir -p $(TSVDIR)
#	@MakeDb.py imgt -i $< --outdir $(TSVDIR) --extended 

# ------------------------------------------------------------------------------
# EXTRACT sequences from IRepertoire DB files into FASTA files
# ------------------------------------------------------------------------------
#$(EXTRFAS) : $(FASTADIR)/%.fasta : $(IREPDIR)/%$(IREP_SUFF)
#	@echo Extracting sequences from IRepertoire datafile $< into FASTA file...
#	@mkdir -p $(FASTADIR)
#	@gawk -f $(SCRIPTDIR)/irep2fasta.gawk $< > $@

# ------------------------------------------------------------------------------
# ANNOTATE extracted FASTA sequences from IRepertoire DB files with IgBLAST
# ------------------------------------------------------------------------------
#$(ANNOTATE) : $(TSVDIR)/%.tsv : $(FASTADIR)/%.fasta
#	@echo Annotate $(SPECIES) $(CHAIN) sequences from $< to AIRR dataset...
#	@mkdir -p $(TSVDIR)
#	@export IGDATA=$(IGDATADIR); \
#         igblastn -germline_db_V $(IGDATADIR)/imgt/$(SPECIES)/IGHV \
#                  -germline_db_D $(IGDATADIR)/imgt/$(SPECIES)/IGHD \
#                  -germline_db_J $(IGDATADIR)/imgt/$(SPECIES)/IGHJ \
#                  -auxiliary_data optional_file/$(SPECIES)_gl.aux \
#                  -domain_system imgt -ig_seqtype $(LOCUS) -organism $(SPECIES) \
#                  -outfmt 19 \
#                  -num_threads=3 \
#                  -query $< \
#                  -out $@

# ------------------------------------------------------------------------------
# PREPARE dataset: merge annotated sequences with external (IRepertoire) data
# ------------------------------------------------------------------------------
$(PREPDIR)/%_prepared.tsv : $(TSVDIR)/%_db-pass.tsv
	@echo "Preparing (merging IMGT and IRep) dataset $* ..."
	@mkdir -p $(PREPDIR)
	@Rscript --vanilla $(SCRIPTDIR)/preprep.R $< $(IREPDIR)/$*$(IREP_SUFF) $@

# ------------------------------------------------------------------------------
# COLLAPSE dataset: pre-collapse, collapse and discard unproductive sequences
# ------------------------------------------------------------------------------
$(COLLDIR)/%_collapsed.tsv : $(PREPDIR)/%_prepared.tsv
	@echo "Collapsing (pre-collapsing, collapsing and filtering) dataset $* ..."
	@mkdir -p $(COLLDIR)
	@Rscript --vanilla $(SCRIPTDIR)/collapserep.R $< $@

# ------------------------------------------------------------------------------
# CALCULATE gene (V,D,J) usages
# ------------------------------------------------------------------------------
$(GENEUSG) : $(DESTDIR)/%_genes.xlsx : $(COLLDIR)/%_collapsed.tsv
	@echo "Calculating V(D)J gene usages in repertoire $* ..."
	@mkdir -p $(DESTDIR)
	@Rscript --vanilla $(SCRIPTDIR)/geneusage.R $< $@
	@Rscript --vanilla $(SCRIPTDIR)/geneusage_m.R $< $(DESTDIR)/$*_genes.tsv

# ------------------------------------------------------------------------------
# GENERATE unified color codes for collected V-gene list
# ------------------------------------------------------------------------------
$(COLCODES) : $(GENEUSG)
	@echo "Collecting V-genes from dataset and assigning color codes..."
	@Rscript --vanilla $(SCRIPTDIR)/vgenecolors.R $(COLCODES) $(GENEUSG)

# ------------------------------------------------------------------------------
# GENERATE treemaps based on gene usage diversity and abundance
# ------------------------------------------------------------------------------
$(TREEMAPS) : $(DESTDIR)/%_tmap.tiff : $(COLLDIR)/%_collapsed.tsv $(COLCODES)
	@echo "Generating treemaps for repertoires..."
	@Rscript --vanilla $(SCRIPTDIR)/treemap.R $(COLCODES) $< $@

# ------------------------------------------------------------------------------
# CALCULATE Hamming distances, distance thresholds and generate histograms
# ------------------------------------------------------------------------------
$(DISTTH) : $(DESTDIR)/%_dist_thrs.data: $(COLLDIR)/%_collapsed.tsv
	@echo Calculating Hamming distance threshold for clonal clustering on $< ...
	@mkdir -p $(DESTDIR)
	@Rscript --vanilla $(SCRIPTDIR)/findthreshold.R $< $@

# ------------------------------------------------------------------------------
# DEFINE CLONES in prepared repertoires
# ------------------------------------------------------------------------------
$(DEFCLONES) : $(CLONEDIR)/%_allclone.tsv : $(COLLDIR)/%_collapsed.tsv $(DESTDIR)/%_dist_thrs.data
	@echo Defining clones in repertoire $* ...
	@mkdir -p $(CLONEDIR)
	@$(eval DIST=$(shell cat $(DESTDIR)/$*_dist_thrs.data))
	@DefineClones.py --nproc 6 --mode gene --act set --norm len --model ham \
	                 --dist $(DIST) -d $< -o $@

# ------------------------------------------------------------------------------
# GENERATE clonal treemaps based on gene usage diversity and abundance
# ------------------------------------------------------------------------------
$(TMAPS_CL) : $(DESTDIR)/%_tmap_cl.tiff : $(CLONEDIR)/%_allclone.tsv $(COLCODES)
	@echo "Generating clonal treemaps for repertoires..."
	@Rscript --vanilla $(SCRIPTDIR)/treemap_clonal.R $(COLCODES) $< $@

# ------------------------------------------------------------------------------
# CALCULATE (intra-)compartmental statistics and measures
# ------------------------------------------------------------------------------
$(CSTATS) : $(CLONEDIR)/%_cstats.tsv : $(CLONEDIR)/%_allclone.tsv
	@echo "Generating compartmental statistics/measures of $*"
	@Rscript --vanilla $(SCRIPTDIR)/cstats.R $< $@

# ------------------------------------------------------------------------------
# GENERATE GERMLINES in prepared repertoires
# ------------------------------------------------------------------------------
$(GENGLS) : $(CLONEDIR)/%_withgl.tsv : $(CLONEDIR)/%_allclone.tsv
	@echo Generating germline sequences for clonotypes in repertoire $* ...
	@CreateGermlines.py -d $< -g dmask --cloned \
	                    -r $(IGDATADIR)/imgt/$(SPECIES)/IGHV.fasta \
	                       $(IGDATADIR)/imgt/$(SPECIES)/IGHD.fasta \
	                       $(IGDATADIR)/imgt/$(SPECIES)/IGHJ.fasta \
	                    -o $@

# ------------------------------------------------------------------------------
# COUNT R/S mutation counts and frequencies in prepared repertoires
# ------------------------------------------------------------------------------
$(RSMUT) : $(CLONEDIR)/%_rsmut.xlsx : $(CLONEDIR)/%_withgl.tsv
	@echo Counting/calculating R/S mutation frequencies/rates in repertoire $* ...
	@Rscript --vanilla $(SCRIPTDIR)/R_S_mutations.R $< $@

# ------------------------------------------------------------------------------
# CALCULATE BASELINe (quantifying antigen selection pressure)
# ------------------------------------------------------------------------------
$(AGSEL) : $(CLONEDIR)/%_agsel.tsv : $(CLONEDIR)/%_withgl.tsv
	@echo Quantifying antigen selection pressure in repertoire $* ...
	@Rscript --vanilla $(SCRIPTDIR)/agselection.R $< $@

# ------------------------------------------------------------------------------
# GENERATE COMPARTMENTAL AIRR network diagrams (graphs)
# ------------------------------------------------------------------------------
$(CGRAPH) : $(GRAPHDIR)/%.gml : $(CLONEDIR)/%_agsel.tsv
	@echo Generating compartmental AIRR network diagram of repertoire $* ...
	@mkdir -p $(GRAPHDIR)
	@Rscript --vanilla $(SCRIPTDIR)/generate_nd.R $< $@

# ------------------------------------------------------------------------------
# PLOT AIRR network diagrams (graphs)
# ------------------------------------------------------------------------------
$(VISND) : $(GRAPHDIR)/%_AIRRND.svg : $(GRAPHDIR)/%.gml
	@echo Visualizing compartmental AIRR network diagram of repertoire $* ...
	@mkdir -p $(GRAPHDIR)
	@Rscript --vanilla $(SCRIPTDIR)/visualize_nd.R $< $@ \
                           $(GRAPHDIR)/$*_AIRRND_agsel.svg \
                           $(GRAPHDIR)/$*_AIRRND_itypes.svg \
                           $(GRAPHDIR)/$*_AIRRND_mutations.svg


#extrtsv:   $(EXTRTSV)
#extrfas:   $(EXTRFAS)
#annotate:  $(ANNOTATE)
preprep:   $(PREPREP)
collrep:   $(COLLREP)
geneusg:   $(GENEUSG)
colcodes:  $(COLCODES)
distth:    $(DISTTH)
treemaps:  $(TREEMAPS)
defclones: $(DEFCLONES)
tmaps_cl:  $(TMAPS_CL)
cstats:    $(CSTATS)
gengls:    $(GENGLS)
rsmut:     $(RSMUT)
agsel:     $(AGSEL)
cgraph:    $(CGRAPH)
visnd:     $(VISND)

clean:
	@rm -rf $(DESTDIR) $(CLONEDIR) $(GRAPHDIR)
	@echo clean.

debug:
	@echo "IMGT: $(IMGT)"
	@echo "IREP: $(IREP)"
	@echo "PREP: $(PREP)"
	@echo "COLL: $(COLL)"
	@echo "REPS: $(REPS)"
	@echo "EXTRTSV: $(EXTRTSV)"
	@echo "EXTRFAS: $(EXTRFAS)"
	@echo "ANNOTATE: $(ANNOTATE)"
	@echo "PREPREP: $(PREPREP)"
	@echo "COLLREP: $(COLLREP)"
	@echo "GENEUSG: $(GENEUSG)"
	@echo "COLCODES: $(COLCODES)"
	@echo "DISTTH: $(DISTTH)"
	@echo "TREEMAPS: $(TREEMAPS)"
	@echo "DEFCLONES: $(DEFCLONES)"
	@echo "TMAPS_CL: $(TMAPS_CL)"
	@echo "CSTATS: $(CSTATS)"
	@echo "GENGLS: $(GENGLS)"
	@echo "RSMUT: $(RSMUT)"
	@echo "AGSEL: $(AGSEL)"
	@echo "CGRAPH: $(CGRAPH)"
	@echo "VISND: $(VISND)"