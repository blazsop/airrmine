# ------------------------------------------------------------------------------
# ICT POST analytic pipeline  v1.00
#
# by Peter Blazso
# 2021-06-21
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

# source/target directories
IGDATADIR := ../igdata
SCRIPTDIR := scripts
LINDIR    := lin
GRAPHDIR  := graph

# filename suffices
TSV_SUFF := .tsv

# source/target files
LINFILES  = $(wildcard $(LINDIR)/*$(TSV_SUFF))  
LINS      = $(notdir $(subst $(TSV_SUFF),,$(LINFILES)))

# destination datafile lists, single (aggregated) data files
GENGLS    = $(addprefix $(GRAPHDIR)/, $(addsuffix _withgl.tsv, $(LINS)))
CLSTATS   = $(addprefix $(GRAPHDIR)/, $(addsuffix _clstats.tsv, $(LINS)))
AGSEL     = $(addprefix $(GRAPHDIR)/, $(addsuffix _agsel.tsv, $(LINS)))
TREEGEN   = $(addprefix $(GRAPHDIR)/, $(addsuffix .gml, $(LINS)))
TREEVIS   = $(addprefix $(GRAPHDIR)/, $(addsuffix _tree.svg, $(LINS)))

# ------------------------------------------------------------------------------
.PHONY: all clean

all: treevis
	@echo done.


# ------------------------------------------------------------------------------
# GENERATE GERMLINES in prepared lineages
# ------------------------------------------------------------------------------
$(GENGLS) : $(GRAPHDIR)/%_withgl.tsv : $(LINDIR)/%.tsv
	@echo Generating germline sequences for clonotypes in repertoire $* ...
	@mkdir -p $(GRAPHDIR)
	@CreateGermlines.py -d $< -g full --cloned \
	                    -r $(IGDATADIR)/imgt/$(SPECIES)/IGHV.fasta \
	                       $(IGDATADIR)/imgt/$(SPECIES)/IGHD.fasta \
	                       $(IGDATADIR)/imgt/$(SPECIES)/IGHJ.fasta \
	                    -o $@

# ------------------------------------------------------------------------------
# CALCULATE BASELINe (quantifying antigen selection pressure)
# ------------------------------------------------------------------------------
$(AGSEL) : $(GRAPHDIR)/%_agsel.tsv : $(GRAPHDIR)/%_withgl.tsv
	@echo Quantifying antigen selection pressure in lineage $* ...
	@mkdir -p $(GRAPHDIR)
	@Rscript --vanilla $(SCRIPTDIR)/agselection.R $< $@

# ------------------------------------------------------------------------------
# CALCULATE CLONOTYPE measures (e.g. clonotype contribution)
# ------------------------------------------------------------------------------
$(CLSTATS) : $(GRAPHDIR)/%_clstats.tsv : $(GRAPHDIR)/%_agsel.tsv
	@echo Calculating clonotype statistics in lineage $* ...
	@mkdir -p $(GRAPHDIR)
	@Rscript --vanilla $(SCRIPTDIR)/clstats.R $< $@

# ------------------------------------------------------------------------------
# GENERATE lineage tree
# ------------------------------------------------------------------------------
$(TREEGEN) : $(GRAPHDIR)/%.gml : $(GRAPHDIR)/%_agsel.tsv
	@echo Generating graph from lineage $* ...
	@mkdir -p $(GRAPHDIR)
	@Rscript --vanilla $(SCRIPTDIR)/treegen.R $< $@

# ------------------------------------------------------------------------------
# PLOT lineage tree
# ------------------------------------------------------------------------------
$(TREEVIS) : $(GRAPHDIR)/%_tree.svg : $(GRAPHDIR)/%.gml
	@echo Visualizing graph as tree of lineage $* ...
	@mkdir -p $(GRAPHDIR)
	@Rscript --vanilla $(SCRIPTDIR)/treevis.R $< $@ $(GRAPHDIR)/$*_tree_agsel.svg


gengls:    $(GENGLS)
agsel:     $(AGSEL)
clstats:   $(CLSTATS)
treegen:   $(TREEGEN)
treevis:   $(TREEVIS)

clean:
	@rm -rf $(GRAPHDIR)
	@echo clean.

debug:
	@echo "LINFILES: $(LINFILES)"
	@echo "LINS: $(LINS)"
	@echo "GENGLS: $(GENGLS)"
	@echo "AGSEL: $(AGSEL)"
	@echo "CLSTATS: $(CLSTATS)"
	@echo "TREEGEN: $(TREEGEN)"
	@echo "TREEVIS: $(TREEVIS)"
