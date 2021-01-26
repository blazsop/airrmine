# ------------------------------------------------------------------------------
# AIRR analytic pipeline system
#
# by Peter Blazso
# 2020-09-28
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
IGDATADIR := igdata
SCRIPTDIR := scripts
IREPDIR   := irep
FASTADIR  := fasta
TSVDIR    := tsv
PREPDIR   := prep
DESTDIR   := dest
CLONEDIR  := clone
GRAPHDIR  := graph

# filename suffices
IREP_SUFF := _irep.csv

# source/target files
IREP      = $(wildcard $(IREPDIR)/*$(IREP_SUFF)) # individual IREP output files
REPS      = $(notdir $(subst $(IREP_SUFF),,$(IREP))) # raw repertoire names

# destination datafile lists, single (aggregated) data files
#EXTRFAS   = $(addprefix $(FASTADIR)/, $(addsuffix .fasta, $(REPS)))
#ANNOTATE  = $(addprefix $(TSVDIR)/,   $(addsuffix .tsv, $(REPS)))
PREPREP   = $(addprefix $(PREPDIR)/,  $(addsuffix _prepared.tsv, $(REPS)))
GENEUSG   = $(addprefix $(DESTDIR)/,  $(addsuffix _genes.xlsx, $(REPS)))
COLCODES  = $(DESTDIR)/v_gene_colorcodes.tsv
TREEMAPS  = $(addprefix $(DESTDIR)/,  $(addsuffix _tmap.tiff, $(REPS)))
DISTTH    = $(addprefix $(DESTDIR)/,  $(addsuffix _dist_thrs.data, $(REPS)))
DEFCLONES = $(addprefix $(CLONEDIR)/, $(addsuffix _allclone.tsv, $(REPS)))
CSTATS    = $(addprefix $(CLONEDIR)/, $(addsuffix _cstats.tsv, $(REPS)))
ICSTATS   = $(addprefix $(CLONEDIR)/, $(addsuffix _icstats.tsv, $(REPS)))
GENGLS    = $(addprefix $(CLONEDIR)/, $(addsuffix _withgl.tsv, $(REPS)))
RSMUT     = $(addprefix $(CLONEDIR)/, $(addsuffix _rsmut.xlsx, $(REPS)))
AGSEL     = $(addprefix $(CLONEDIR)/, $(addsuffix _agsel.tsv, $(REPS)))
CGRAPH    = $(addprefix $(GRAPHDIR)/, $(addsuffix .gml, $(REPS)))
ICGRAPH   = $(addprefix $(GRAPHDIR)/, $(addsuffix _ic.gml, $(REPS)))
VISND     = $(addprefix $(GRAPHDIR)/, $(addsuffix _AIRRND.svg, $(REPS)))
ICVISND   = $(addprefix $(GRAPHDIR)/, $(addsuffix _AIRRND_ic.svg, $(REPS)))

# ------------------------------------------------------------------------------
.PHONY: all clean

all: cgraph
	@echo done.

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
# PREPARE (merge, collapse and subsample) dataset
# ------------------------------------------------------------------------------
$(PREPDIR)/%_prepared.tsv : $(TSVDIR)/%.tsv
	@echo "Preparing (merging, collapsing and subsampling) dataset $* ..."
	@mkdir -p $(PREPDIR)
	@Rscript --vanilla $(SCRIPTDIR)/preprep.R $< $(IREPDIR)/$*$(IREP_SUFF) $@

# ------------------------------------------------------------------------------
# CALCULATE gene (V,D,J) usages
# ------------------------------------------------------------------------------
$(GENEUSG) : $(DESTDIR)/%_genes.xlsx : $(PREPDIR)/%_prepared.tsv
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
$(TREEMAPS) : $(DESTDIR)/%_tmap.tiff : $(PREPDIR)/%_prepared.tsv $(COLCODES)
	@echo "Generating treemaps for repertoires..."
	@Rscript --vanilla $(SCRIPTDIR)/treemap.R $(COLCODES) $< $@


# ------------------------------------------------------------------------------
# CALCULATE Hamming distances, distance thresholds and generate histograms
# ------------------------------------------------------------------------------
$(DISTTH) : $(DESTDIR)/%_dist_thrs.data: $(PREPDIR)/%_prepared.tsv
	@echo Calculating Hamming distance threshold for clonal clustering on $< ...
	@mkdir -p $(DESTDIR)
	@Rscript --vanilla $(SCRIPTDIR)/findthreshold.R $< $@

# ------------------------------------------------------------------------------
# DEFINE CLONES in prepared repertoires
# ------------------------------------------------------------------------------
$(DEFCLONES) : $(CLONEDIR)/%_allclone.tsv : $(PREPDIR)/%_prepared.tsv $(DESTDIR)/%_dist_thrs.data
	@echo Defining clones in repertoire $* ...
	@mkdir -p $(CLONEDIR)
	@$(eval DIST=$(shell cat $(DESTDIR)/$*_dist_thrs.data))
	@DefineClones.py --nproc 6 --mode gene --act set --norm len --model ham \
	                 --dist $(DIST) -d $< -o $@

# ------------------------------------------------------------------------------
# CALCULATE (intra-)compartmental statistics and measures
# ------------------------------------------------------------------------------
$(CSTATS) : $(CLONEDIR)/%_cstats.tsv : $(CLONEDIR)/%_allclone.tsv
	@echo "Generating compartmental statistics/measures of $*"
	@Rscript --vanilla $(SCRIPTDIR)/cstats.R $< $@

# ------------------------------------------------------------------------------
# CALCULATE inter-compartmental connections related to lineage count
# ------------------------------------------------------------------------------
$(ICSTATS) : $(CLONEDIR)/%_icstats.tsv : $(CLONEDIR)/%_allclone.tsv
	@echo "Calculating inter-compartmental statistics in repertoire $* ..."
	@Rscript --vanilla $(SCRIPTDIR)/icstats.R $< $@

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
$(CGRAPH) : $(GRAPHDIR)/%.gml : $(CLONEDIR)/%_withgl.tsv
	@echo Generating compartmental AIRR network diagram of repertoire $* ...
	@mkdir -p $(GRAPHDIR)
	@Rscript --vanilla $(SCRIPTDIR)/netdiag.R $< $@

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

# ------------------------------------------------------------------------------
# GENERATE INTER-COMPARTMENTAL AIRR network diagrams (graphs)
# ------------------------------------------------------------------------------
$(ICGRAPH) : $(GRAPHDIR)/%_ic.gml : $(CLONEDIR)/%_agsel.tsv
	@echo Generating inter-compartmental AIRR network of repertoire $* ...
	@mkdir -p $(GRAPHDIR)
	@Rscript --vanilla $(SCRIPTDIR)/netdiag_ic.R $< $@

# ------------------------------------------------------------------------------
# PLOT AIRR network diagrams (graphs)
# ------------------------------------------------------------------------------
$(ICVISND) : $(GRAPHDIR)/%_AIRRND_ic.svg : $(GRAPHDIR)/%_ic.gml
	@echo Visualizing inter-compartmental AIRR network diagrams of repertoire $* ...
	@mkdir -p $(GRAPHDIR)
	@Rscript --vanilla $(SCRIPTDIR)/visualize_nd.R $< $@ \
                           $(GRAPHDIR)/$*_AIRRND_ic_agsel.svg \
                           $(GRAPHDIR)/$*_AIRRND_ic_itypes.svg \
                           $(GRAPHDIR)/$*_AIRRND_ic_mutations.svg


extrfas:   $(EXTRFAS)
annotate:  $(ANNOTATE)
preprep:   $(PREPREP)
geneusg:   $(GENEUSG)
colcodes:  $(COLCODES)
distth:    $(DISTTH)
treemaps:  $(TREEMAPS)
defclones: $(DEFCLONES)
cstats:    $(CSTATS)
icstats:   $(ICSTATS)
gengls:    $(GENGLS)
rsmut:     $(RSMUT)
agsel:     $(AGSEL)
icgraph:   $(ICGRAPH)
cgraph:    $(CGRAPH)
visnd:     $(VISND)
icvisnd:   $(ICVISND)

clean:
	@rm -rf $(DESTDIR) $(CLONEDIR) $(GRAPHDIR)
	@echo clean.

debug:
	@echo "IREP: $(IREP)"
	@echo "REPS: $(REPS)"
	@echo "EXTRFAS: $(EXTRFAS)"
	@echo "ANNOTATE: $(ANNOTATE)"
	@echo "PREPREP: $(PREPREP)"
	@echo "GENEUSG: $(GENEUSG)"
	@echo "DISTTH: $(DISTTH)"
	@echo "DEFCLONES: $(DEFCLONES)"
	@echo "CSTATS: $(CSTATS)"
	@echo "ICSTATS: $(ICSTATS)"
	@echo "GENGLS: $(GENGLS)"
	@echo "AGSEL: $(AGSEL)"
	@echo "CGRAPH: $(CGRAPH)"
	@echo "ICGRAPH: $(ICGRAPH)"
	@echo "VISND: $(VISND)"
	@echo "ICVISND: $(ICVISND)"
